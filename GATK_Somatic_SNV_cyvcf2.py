import sys, math, numpy
from argparse import RawTextHelpFormatter, ArgumentParser
from collections import defaultdict
import cyvcf2


__author__ = "Ryan Smith (ryanpsmith@wustl.edu)"
__version__ = "$Revision: 0.0.1 $"
__date__ = "$Date: 2016-2-8 13:45 $"



def get_args():
    parser = ArgumentParser(
        formatter_class=RawTextHelpFormatter, add_help=False)

    parser.add_argument('-i', metavar='VCF_IN', 
                        help='Input VCF (stdin)')

    parser.add_argument('-o', metavar='VCF_OUT', 
                        help='Output VCF (stdout)')

    parser.add_argument('-s', metavar='SAMPLE_MAP',
                        help='Sample Map File: \n (< Animal >  < Sample >  < Source >  < Experiment >  < Case (SCNT/Control) >)',
                        required=True)

    parser.add_argument('--indel', 
                        help='Indel mode', action='store_true')

    parser.add_argument('--fnr', 
                        help='FNR test mode', action='store_true')

    parser.add_argument('--vaf', 
                        help='Autosomal VAF Minimum for call')

    parser.add_argument('--svaf', 
                        help='Sex VAF Minimum for call')

    parser.add_argument('-h','--help', action='help') 


    # parse the arguments
    args = parser.parse_args()
    
    # bail if no BAM file
    if args.i is None:
        if sys.stdin.isatty():
            parser.print_help()
            exit(1)
    
    # send back the user input
    return args

def unphred(array):
    return numpy.power(10., numpy.divide(array, -10.))

def get_PL_ratio(PL, AAG_IDX, RR_IDX, case=True):

    #get AAG and RR PL values
    AAG_PL = unphred(int(PL[AAG_IDX]))
    RR_PL = unphred(int(PL[RR_IDX]))

    num = AAG_PL
    denom = RR_PL

    if not case:
        num = RR_PL
        denom = AAG_PL

    try:
        ratio = num / denom
    except ZeroDivisionError:
        ratio = float('inf')

    return ratio

def load_sample_map(infile):
    infile = open(infile, 'r')
    sample_map = defaultdict(dict)
    header = False

    for line in infile:

        line = line.strip().split("\t")

        if not header:
            assert line == ['Sample', 'Animal', 'Source', 'Experiment', 'Case'], \
                "Invalid sample map header, must match: \"Animal Sample Source Experiment Case\""
            header = line
            continue

        assert line[4] == 'SCNT' or line[4] == 'Control', \
            "Invalid sample map. Case must be one of \"SCNT\" or \"Control\""

        sample_map[line[0]]['Animal'] = line[1]
        sample_map[line[0]]['Source'] = line[2]
        sample_map[line[0]]['Experiment'] = line[3]
        sample_map[line[0]]['Case'] = line[4]

    infile.close()
    return sample_map




# ============================================
# driver
# ============================================

class Usage(Exception):
    def __init__(self, msg):
        self.msg = msg

def main():

    #parse arguments
    args = get_args()
    indel = args.indel
    fnr = args.fnr

    #should look specifically at X chrom VAFs in males
    VAF=0.30
    SVAF=0.95
    if args.vaf:    VAF=float(args.vaf)
    if args.svaf:   SVAF=float(args.svaf)

    try:
      sample_map = load_sample_map(args.s)
    except:
        print("Error: Misformatted Sample map (-s)")
        exit(1)

    reader = cyvcf2.VCF(args.i, gts012=True)
    #gts012 sets the uncalled GTs to 3

    #get list of samples
    samples = reader.samples

    #add the two new INFO tags
    reader.update("UNIQ", "String", 1, "Sample(s) with unique somatic variant")
    reader.update("UAB", "Float", 1, "Allele Balance in UNIQ sample")
    reader.update("TYPE", "String", 1, "Varant Type (SNPS: TS/TV) (INDELS: INS/DEL)")
    reader.update("TISSUE", "String", 1, "Source Tissue Type")
    reader.update("CASE", "String", 1, "Control or SCNT?")
    reader.update("EXPT", "String", 1, "Experiment")

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
        MAX_VAF = 0.05
        if indel:
            MAX_VAF = 0

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
            if not unique:
                break

            #GT call must match GT of interest
            if GTs[i] != AAG:
               continue

            #criteria for presence in given sample
            if DEPTHS[i] >= min_depth \
                and DEPTHS[i] <= max_depth \
                and ABs[i] >= VAF \
                and AAG_RR_ratios[i] >= AAG_RR_MIN:
                

                for j in range(len(samples)):
                    control = False
                    #other = False
                    if i == j:
                        continue

                    #check for case/control/other
                    # control is only control from same animal
                    # other is any other sample
                    if sample_map[samples[j]]['Animal'] == sample_map[samples[i]]['Animal']:
                        control = True

                    ratio_min = 1
                    if control:
                        ratio_min = RR_AAG_MIN

                    #criteria for presence in other samples (need to make this control vs ALL j_ratio)
                    if DEPTHS[j] < min_depth \
                        or DEPTHS[j] > max_depth \
                        or ABs[j] > MAX_VAF \
                        or RR_AAG_ratios[j] < ratio_min:
                        unique = False
                        break

                if unique:
                    #set new info fields and write record
                    var.INFO['TISSUE'] = sample_map[samples[i]]['Source']
                    var.INFO['CASE'] = sample_map[samples[i]]['Case']
                    var.INFO['EXPT'] = sample_map[samples[i]]['Experiment']
                    var.INFO["UNIQ"] = samples[i]
                    var.INFO["UAB"] = str(numpy.around(ABs[i], 3))
                    writer.write_record(var)

if __name__ == "__main__":
    try:
        sys.exit(main())
    except IOError, e:
        if e.errno != 32:  # ignore SIGPIPE
            raise


