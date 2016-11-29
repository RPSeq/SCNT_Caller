import pysam, sys


infile = sys.argv[1]
vcf_file = pysam.VariantFile(infile, 'r')



samples = vcf_file.header.samples
num = len(samples)
min_su = 5


vcf_out = pysam.VariantFile("-", 'w', header=vcf_file.header)

sample_outs = [pysam.VariantFile(sample, 'w', header=vcf_file.header) for sample in samples]

matrix = [[0 for i in samples] for j in samples]
uniques = [0 for i in samples]

# count = 0
for var in vcf_file:
    LP = int(var.info['LP'])
    RP = int(var.info['RP'])

    if var.filter.keys()[0] == 'PASS' and (LP + RP) > 9:
        
        for i in range(num):
            unique = True
            i_GT = str(var).split()[9+i].split(":")[0]

            i_GL = var.samples[samples[i]]['GL']
            i_GL = [float(x) for x in i_GL]


            if i_GT == "0/1":
                for j in range(num):

                    if j != i:
                        
                        j_GL = var.samples[samples[j]]['GL']
                        j_GL = [float(x) for x in j_GL]

                        j_GT = str(var).split()[9+j].split(":")[0]

                        if j_GT != "0/0" or j_GL[0] < -0.05:
                            unique = False
                            break

                if unique:
                    vcf_out.write(var)
                    sample_outs[i].write(var)





for out in sample_outs:
    out.close()

vcf_file.close()
vcf_out.close()








