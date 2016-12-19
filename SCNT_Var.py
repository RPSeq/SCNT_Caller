#!/usr/bin/env python

import sys, numpy, pyfaidx
from argparse import RawTextHelpFormatter, ArgumentParser
from collections import defaultdict, Counter
import cyvcf2


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
                        help='Autosomal VAF Minimum for SNP call')

    parser_snp.add_argument('--svaf', 
                        help='Sex VAF Minimum for SNP call')

    parser_snp.add_argument('--low', action='store_true',
                        help='Annotate instead of discarding low VAF calls')

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
            assert line == ['Sample', 'Animal', 'Source', 'Experiment', 'Case'], \
                "Invalid sample map header, must match: \"Animal Sample Source Experiment Case\""
            header = line
            continue

        sample, animal, source, experiment, case = line

        assert case == 'SCNT' or case == 'Control', \
            "Invalid sample map. Case must be one of \"SCNT\" or \"Control\""

        sample_map[sample]['Animal'] = animal
        sample_map[sample]['Source'] = source
        sample_map[sample]['Experiment'] = experiment
        sample_map[sample]['Case'] = case

        animal_map[animal][case].append(sample)

    infile.close()
    return sample_map, animal_map



def same_animal(smap, samples, i , j):

    try:
        if smap[samples[j]]['Animal'] == smap[samples[i]]['Animal']:
            return True
    except KeyError:
        sys.stderr.write("Error: IDs in sample map do not match input .vcf samples")
        exit(0)

    return False

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
    SVAF=0.95
    if args.vaf:    VAF=float(args.vaf)
    if args.svaf:   SVAF=float(args.svaf)

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
    reader.update("FILTER", "String", 1, "VAF Filter PASS/FAIL")
    reader.update("CONTEXT", "String", 1, "Trinucleotide Context")

    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    AAG_RR_MIN = numpy.power(10.,10.)
    RR_AAG_MIN = numpy.power(10.,5.)

    min_depth = 10
    max_depth = 250

    #iterate over vars
    for var in reader:

        #don't care about untigs
        if 'GL' in var.CHROM or 'JH' in var.CHROM:
            continue

        #for now, site must be genotyped in ALL samples.
        #could lose some sensitivity here:
        #often have multiple mouse strains in single experiment.
        #if a region is deleted in one strain but not another,
        # there will never be genotypes in that region in the del strain.
        if var.call_rate < 1:
            continue

        # GT ids:
        #   0=HOM_REF
        #   1=HET
        #   2=HOM_ALT
        #   3=UNKNOWN
        #alternate allele genotype is het
        AAG = 1

        #IF FEMALE, X AAG should still be 0/1
        if var.CHROM == "X" or var.CHROM == "Y":
            AAG = 2
            VAF=SVAF

        #max alt allele balance for control samples
        if var.is_snp:
            MAX_VAF = 0.05
        elif var.is_indel:
            MAX_VAF = 0

        else:
            sys.stderr.write("Skipping Variant: Not SNP/Indel")
            continue

        unique = True

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

        #get allele balances
        ABs = numpy.true_divide(ALT_DEPTHS, DEPTHS)

        for i in range(len(samples)):
            vaf_filt = "PASS"
            if not unique:
                break

            #criteria for presence in given sample
            if (DEPTHS[i] >= min_depth and
                DEPTHS[i] <= max_depth and
                AAG_RR_ratios[i] >= AAG_RR_MIN and
                GTs[i] == AAG):

                if ABs[i] < VAF:
                    if args.low:
                        vaf_filt = "FAIL"
                    else:
                        continue
                
                for j in range(len(samples)):
                    if i == j:
                        continue

                    #check for case/control/other
                    # control is only control from same animal
                    # other is any other sample
                    #if not control, RR_AAG_MIN is 1.
                    ratio_min = 1
                    if same_animal(sample_map, samples, i, j):
                        ratio_min = RR_AAG_MIN


                    #criteria for presence in other samples 
                    #   (need to make this control vs ALL j_ratio)
                    if (DEPTHS[j] < min_depth or
                        DEPTHS[j] > max_depth or
                        ABs[j] > MAX_VAF or
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
                    var.INFO['FILTER'] = vaf_filt

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

    #iterate over vars
    for var in reader:

        #false by default
        gss = False

        animals = []
        present = []

        #max alt allele balance for control samples
        if not var.is_snp:
            sys.stderr.write("Skipping Variant: Not SNP")
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

        if gss:
            var.INFO['ANIMAL'] = ",".join(animals)
            var.INFO['PRESENT'] = ",".join(present)
            writer.write_record(var)

    for sample in sorted(counter.keys()):
        outstr = "\t".join([sample, sample_map[sample]['Case'], str(counter[sample])])
        counts_out.write(outstr+"\n")

    writer.close()
    counts_out.close()


def sv(args):


    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)

    #get list of samples
    samples = reader.samples  
    
    #read sample map
    sample_map = load_sample_map(args.m)

    #add the new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")
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

                        if ABs[j][0] > 0.05:
                            unique = False
                            break

                if unique:
                    #set new info fields
                    var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                    var.INFO['CASE'] = sample_map[samples[i]]['Case']
                    var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                    var.INFO['UNIQ'] = samples[i]
                    var.INFO['UAB'] = str(numpy.around(ABs[i][0], 3))
                    # var.INFO['FILTER'] = filt
                    writer.write_record(var)
    writer.close()


def mei(args):

    #gts012 sets the uncalled GTs to 3
    reader = cyvcf2.VCF(args.i, gts012=True)

    #get list of samples
    samples = reader.samples  
    
    #read sample map
    sample_map = load_sample_map(args.m)

    #add the new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")
    reader.update("FILTER", "String", 1, "VAF Filter PASS/FAIL")

    if not args.o:
        writer = cyvcf2.Writer("/dev/stdout", reader)
    else:
        writer = cyvcf2.Writer(args.o, reader)

    min_su = 9
    for var in reader:
        LP = var.INFO['LP']
        RP = var.INFO['RP']

        if var.FILTER == 'PASS' and (LP + RP) > min_su:
            GTs = var.gt_types
            RR_PLs = var.gt_phred_ll_homref
            AAG_PLs = var.gt_phred_ll_het

            for i in range(len(samples)):
                unique = True

                #heterozygote
                if GTs[i] == 1:
                    for j in range(num):
                        if j != i:
                
                            if GTs[j] != 0 or RR_PLs[j] < -0.05:
                                unique = False
                                break

                    if unique:
                        #set new info fields
                        var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                        var.INFO['CASE'] = sample_map[samples[i]]['Case']
                        var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                        var.INFO['UNIQ'] = samples[i]
                        var.INFO['UAB'] = str(numpy.around(ABs[i], 3))
                        # var.INFO['FILTER'] = filt
                        writer.write_record(var)

    writer.close()



if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


