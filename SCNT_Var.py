#!/usr/bin/env python

import sys, numpy, pyfaidx, cyvcf2
from argparse import RawTextHelpFormatter, ArgumentParser
from collections import defaultdict, Counter


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"


#need to modify this for submodules:
#need snp, indel (or SNPDEL), mei (melt), sv (lumpy/svtyper)
def get_args():
    parser = ArgumentParser()

    #add subparsers
    subparsers = parser.add_subparsers(title='modes', dest='mode')

    #snpdel specific
    parser_snp = subparsers.add_parser('SNP', help='SNP Calling Mode')
    parser_snp.set_defaults(func=snp)

    parser_fnr = subparsers.add_parser('FNR', help='SNP FNR Mode')
    parser_fnr.set_defaults(func=snp_fnr)

    #MEI mode
    parser_mei = subparsers.add_parser('MEI', help='MEI Calling Mode')
    parser_mei.set_defaults(func=mei)

    #SV mode
    parser_sv = subparsers.add_parser('SV', help='SV Calling Mode')
    parser_sv.set_defaults(func=sv)

    #list to hold required group for parser snp, mei, sv
    reqs = []
    for p in [parser_snp, parser_fnr, parser_mei, parser_sv]:
        req = p.add_argument_group('Required Options')
        p.add_argument('-i', metavar='input', 
                            help='Input VCF [stdin]')

        p.add_argument('-o', metavar='output', 
                            help='Output VCF [stdout]')

        req.add_argument('-m', metavar='map',
                            help='Sample Map Columns: \n(Animal, Sample, Source, Experiment, Case)',
                            required=True)
        reqs.append(req)

    #add REF to required opts for parser_snp function
    reqs[0].add_argument('-r', metavar='REF', 
                        help='Reference Genome (.fa, indexed)',
                        required=True)

    parser_snp.add_argument('--vaf', 
                        help='VAF Minimum for SNP in autosome and X of females [0.30]')

    parser_snp.add_argument('--mvaf', 
                        help='VAF Minimum for SNP in X/Y of males [0.95]')

    reqs[1].add_argument('-c', metavar='counts',
                        help='Counts output file',
                        required=True)

    # parse the arguments
    args = parser.parse_args()
    
    # bail if no input file or stream
    if args.i is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
        args.i = "-"
    
    # send back the user input
    return args

def unphred(array):
    return numpy.power(10., numpy.divide(array, -10.))

def load_sample_map(infile):

    infile = open(infile, 'r')
    sample_map = defaultdict(dict)
    animal_map = defaultdict(lambda: defaultdict(list))
    header = False

    for line in infile:

        line = line.strip().split("\t")

        if not header:
            assert line == ['Sample', 'Animal', 'Source', 'Experiment', 'Case', 'Sex'], \
                "Invalid sample map header, must match: \"Animal Sample Source Experiment Case Sex\""
            header = line
            continue

        sample, animal, source, experiment, case, sex = line

        assert case == 'SCNT' or case == 'Control', \
            "Invalid sample map. Case must be one of \"SCNT\" or \"Control\""

        assert sex == 'M' or sex == 'F', \
            "Invalid sample map. Sex must be one of \"F\" or \"M\""

        sample_map[sample]['Animal'] = animal
        sample_map[sample]['Source'] = source
        sample_map[sample]['Experiment'] = experiment
        sample_map[sample]['Case'] = case
        sample_map[sample]['Sex'] = sex

        animal_map[animal][case].append(sample)

    infile.close()
    return sample_map, animal_map



def same_animal(smap, samples, i , j):

    try:
        return smap[samples[j]]['Animal'] == smap[samples[i]]['Animal']

    except KeyError:
        sys.stderr.write("Error: IDs in sample map do not match input .vcf samples")
        exit(0)


# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():

    #parse arguments
    args = get_args()
    args.func(args)

def snp(args):

    genome = pyfaidx.Fasta(args.r)

    #should look specifically at X chrom VAFs in males
    VAF=0.30
    MVAF=0.95
    if args.vaf:    VAF=float(args.vaf)
    if args.mvaf:   MVAF=float(args.mvaf)

    sample_map, animal_map = load_sample_map(args.m)

    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)


    #get list of samples
    samples = reader.samples

    #add the two new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TYPE", "String", 1, "Varant Type (SNPS: TS/TV) (INDELS: INS/DEL)")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")
    reader.update("CONTEXT", "String", 1, "Trinucleotide Context")
    reader.update("ANIMAL", "String", 1, "ID of origin animal")
    reader.update("UDP", "Integer", 1, "Depth at uniq site")
    reader.update("AAGR", "Float", 1, "AAG/RR Ratio")
    reader.update("UGQ", "Float", 1, "Uniq Genotype Quality")

    reader.add_filter_to_header({"ID":"LowVAF", "Description":"Somatic VAF below threshold"})
    reader.add_filter_to_header({"ID":"MGP", "Description":"Variant present in MGP"})

    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    AAG_RR_MIN = numpy.power(10.,10.)
    RR_AAG_MIN = numpy.power(10.,5.)

    min_depth = 10
    max_depth = 250

    allosomes = set(["X", "Y"])

    #iterate over vars
    for var in reader:
        

        unique = True

        # set max alt VAF for snp or indel
        if var.is_snp:
            MAX_VAF = 0.05
        elif var.is_indel:
            MAX_VAF = 0.00
        else:
            sys.stderr.write("Skipping Variant: Not SNP/Indel")
            continue


        #get RR and AAG genotype likelihoods
        RR_PLs = unphred(var.gt_phred_ll_homref)
        AAG_PLs = unphred(var.gt_phred_ll_het)

        #get AAG/RR and RR/AAG likelihood ratios
        AAG_RR_ratios = numpy.true_divide(AAG_PLs, RR_PLs)
        RR_AAG_ratios = numpy.true_divide(RR_PLs, AAG_PLs)

        #get genotypes, depths, and alt allele depths
        GTs = var.gt_types
        DEPTHS = var.gt_depths
        ALT_DEPTHS = var.gt_alt_depths
        QUALS = var.gt_quals

        #get allele balances
        ABs = numpy.true_divide(ALT_DEPTHS, DEPTHS)

        for i in range(len(samples)):
            # dont waste time in uneeded loops
            if not unique:
                break

            AAG = 1
            MIN_VAF = VAF

            #change AAG to 1/1 and min vaf to male allosome min [0.95]
            if sample_map[samples[i]]['Sex']=="M" and var.CHROM in allosomes:
                AAG = 2
                MIN_VAF = MVAF

            VAF_FILT = False
            

            #criteria for presence in given sample
            if (max_depth > DEPTHS[i] > min_depth 
                and AAG_RR_ratios[i] >= AAG_RR_MIN 
                and GTs[i] == AAG):

                if ABs[i] < MIN_VAF:
                    VAF_FILT = True
                
                for j in range(len(samples)):
                    if i == j:
                        continue

                    #default RR/AAG min is 1.
                    ratio_min = 1
                    #if same animal
                    if same_animal(sample_map, samples, i, j):
                        if not max_depth > DEPTHS[j] > min_depth:
                            unique = False
                            break

                        #if different case, set ratio_min to control RR_AAG min
                        #thus, SCNT lines for the same animal are treated as the control for controls
                        if sample_map[samples[i]]['Case'] != sample_map[samples[j]]['Case']:
                            ratio_min = RR_AAG_MIN

                    #criteria for failing or presence in other samples 
                    if (ABs[j] > MAX_VAF or
                        RR_AAG_ratios[j] < ratio_min):
                        unique = False
                        break

                if unique:
                    #get trinucleotide context (VCF coords are 1-based)
                    if var.is_snp:
                        context = genome[str(var.CHROM)][var.POS-2:var.POS+1]
                        #ts or tv?
                        tstv = 'Tv'
                        if var.is_transition:
                            tstv = 'Tr'

                        var.INFO['CONTEXT'] = context.seq
                        var.INFO['TYPE'] = tstv

                    var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                    var.INFO['CASE'] = sample_map[samples[i]]['Case']
                    var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                    var.INFO['UNIQ'] = samples[i]
                    var.INFO['UAB'] = str(numpy.around(ABs[i], 3))
                    var.INFO['UDP'] = str(DEPTHS[i])
                    var.INFO['AAGR'] = str(AAG_RR_ratios[i])
                    var.INFO['UGQ'] = str(QUALS[i])
                    var.INFO['ANIMAL'] = sample_map[samples[i]]['Animal']

                    filters = []
                    f = var.FILTER
                    if f:
                        filters = f.split(";")

                    if VAF_FILT:
                        filters.append("LowVAF")

                    try:
                        var.INFO["MGP"]
                        filters.append("MGP")
                    except KeyError:
                        pass

                    if filters:
                        var.FILTER = filters

                    #write record
                    writer.write_record(var)

    writer.close()



def snp_fnr(args): 

    #only concerned with autosome VAF cutoff
    VAF=0.30


    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)

    reader.update("ANIMAL", "String", ".", "Animals with this GSS var")
    reader.update("PRESENT", "String", ".", "SCNT lines detecting this GSS var")

    #open writerr
    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    counts_out = open(args.c, 'w')

    #get list of samples
    samples = reader.samples

    #sample to index map
    stoi = {s: samples.index(s) for s in samples}

    #load sample map and get animal sample groups
    sample_map, animal_map = load_sample_map(args.m)


    AAG_RR_MIN = numpy.power(10.,10.)

    min_depth = 10
    max_depth = 250

    counter = Counter()
    HIcounter = Counter()

    #iterate over vars
    for var in reader:

        #false by default
        gss = False
        PASS = False
        animals = []
        present = []

        #max alt allele balance for control samples
        if not (var.is_snp or var.is_indel):
            sys.stderr.write("Skipping Variant: Not SNP/Indel")
            continue

        #get RR and AAG genotype likelihoods
        RR_PLs = unphred(var.gt_phred_ll_homref)
        AAG_PLs = unphred(var.gt_phred_ll_het)

        #get AAG/RR ratios
        AAG_RR_ratios = numpy.true_divide(AAG_PLs, RR_PLs)

        #get genotypes, depths, and alt allele depths
        GTs = var.gt_types
        DEPTHS = var.gt_depths
        ALT_DEPTHS = var.gt_alt_depths

        #get allele balances
        ABs = numpy.true_divide(ALT_DEPTHS, DEPTHS)

        if not var.FILTER: PASS = True

        for animal, group in animal_map.items():

            #get sample name of control
            control = group['Control'][0]
            SCNTs = group['SCNT']

            #if var not called in the control, continue:
            # if GTs[stoi[control]] == [0,3]:
            #     continue

            ALL = [control] + SCNTs

            #if at least one sample was called het, 
            #   var is present in mouse
            if 1 in [GTs[stoi[x]] for x in ALL]:
                counter[control] += 1
                if PASS: 
                    HIcounter[control] += 1
                gss = True
                control_depth = False
                if (DEPTHS[stoi[control]] >= min_depth and
                    DEPTHS[stoi[control]] <= max_depth):
                    control_depth = True

                animals.append(animal)

                for sample in SCNTs:
                    i = stoi[sample]

                    if (DEPTHS[i] >= min_depth and
                        DEPTHS[i] <= max_depth and
                        control_depth and
                        AAG_RR_ratios[i] >= AAG_RR_MIN and
                        ABs[i] >= VAF and
                        GTs[i] == 1):


                        counter[sample] += 1
                        present.append(sample)
                        if PASS:
                            HIcounter[sample] += 1

        if gss:
            var.INFO['ANIMAL'] = ",".join(animals)
            var.INFO['PRESENT'] = ",".join(present)
            writer.write_record(var)

    counts_out.write("#COUNTS\tSample\tCase\tCount\tHQCount\n")
    for sample in sorted(counter.keys()):
        outstr = "\t".join(["#COUNT", sample, sample_map[sample]['Case'], str(counter[sample]), str(HIcounter[sample])])
        counts_out.write(outstr+"\n")

    counts_out.write("#FNR\tAnimal\tFNR\tHQFNR\n")

    for animal, group in sorted(animal_map.items()):

        #get sample name of control
        control = group['Control'][0]
        SCNTs = group['SCNT']

        present = counter[control]

        rates = []
        hrates = []

        for sample in SCNTs:
            called = counter[sample]
            hcalled = HIcounter[sample]
            rate = 1.0-(called/float(present))
            hrate = 1.0-(hcalled/float(present))
            rates.append(rate)
            hrates.append(hrate)

        a_rate = numpy.mean(rates)
        a_hrate = numpy.mean(hrates)

        counts_out.write("\t".join(["#FNR", animal, str(a_rate), str(a_hrate)])+"\n")


    # writer.close()
    counts_out.close()


def sv(args):


    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)

    #get list of samples
    samples = reader.samples  
    
    #read sample map
    sample_map, animal_map = load_sample_map(args.m)

    #add the new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")
    reader.update("ANIMAL", "String", 1, "ID of origin animal")
    # reader.update("FILTER", "String", 1, "VAF Filter PASS/FAIL")

    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    min_su = 5
    for var in reader:
        SUs = var.format('SU')

        ABs = var.format('AB')

        for i in range(len(samples)):
            if SUs[i][0] >= min_su and ABs[i][0] >= 0.20:
                unique = True

                for j in range(len(samples)):
                    if j != i:

                        if ABs[j][0] > 0.0:
                            unique = False
                            break

                if unique:
                    #set new info fields
                    var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                    var.INFO['CASE'] = sample_map[samples[i]]['Case']
                    var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                    var.INFO['UNIQ'] = samples[i]
                    var.INFO['UAB'] = str(numpy.around(ABs[i][0], 3))
                    var.INFO['ANIMAL'] = sample_map[samples[i]]['Animal']
                    # var.INFO['FILTER'] = filt

                    writer.write_record(var)
    writer.close()


def mei(args):

    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)

    #get list of samples
    samples = reader.samples  
    
    #read sample map
    sample_map, animal_map = load_sample_map(args.m)


    #add the new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    # reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")
    reader.update("ANIMAL", "String", 1, "ID of origin animal")
    # reader.update("PL", "String", 1, "RR")

    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    min_su = 3
    for var in reader:
        LP = var.INFO['LP']
        RP = var.INFO['RP']
        if not var.FILTER and (LP > min_su and RP > min_su):
        # if (LP > min_su and RP > min_su):
            unique = True
            GTs = var.gt_types

            #strange behaviour... MELT PLs read as the input value *10. divide by 10 to correct.
            #should have MELT return positive integers rather than negative floats.
            RR_PLs = unphred(numpy.divide(var.gt_phred_ll_homref, 10.))
            AAG_PLs = unphred(numpy.divide(var.gt_phred_ll_het, 10.))

            #get AAG/RR and RR/AAG likelihood ratios
            AAG_RR_ratios = numpy.true_divide(AAG_PLs, RR_PLs)
            RR_AAG_ratios = numpy.true_divide(RR_PLs, AAG_PLs)

            for i in range(len(samples)):
                if not unique:
                    break

                #heterozygote
                if GTs[i] == 1:
                    for j in range(len(samples)):
                        if j != i:
                
                            if GTs[j] != 0 or RR_PLs[j] < 0.90:
                                unique = False
                                break

                    if unique:
                        # print RR_PLs
                        #set new info fields
                        var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                        var.INFO['CASE'] = sample_map[samples[i]]['Case']
                        var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                        var.INFO['UNIQ'] = samples[i]
                        var.INFO['ANIMAL'] = sample_map[samples[i]]['Animal']
                        writer.write_record(var)

    writer.close()


if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


