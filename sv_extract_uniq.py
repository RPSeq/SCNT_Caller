import pysam, sys

infile = sys.argv[1]
vcf_file = pysam.VariantFile(infile, 'r')
samples = vcf_file.header.samples
num = len(samples)
min_su = 5
vcf_name = "SV_SU"+str(min_su)+"_uniqs.vcf"
vcf_out = pysam.VariantFile(vcf_name, 'w', header=vcf_file.header)
sample_outs = [pysam.VariantFile(sample, 'w', header=vcf_file.header) for sample in samples]

for var in vcf_file:
    for i in range(num):

        i_SU = int(var.samples[samples[i]]['SU'])

        try:
            i_AB = float(var.samples[samples[i]]['AB'])
        except:
            i_AB = 0

        if i_SU >= min_su and i_AB >= 0.15:
            unique = True

            for j in range(num):
                if j != i:

                    try:
                        j_AB = float(var.samples[samples[j]]['AB'])
                    except:
                        j_AB = 0

                    if j_AB > 0:
                        unique = False
                        break

            if unique:
                vcf_out.write(var)
                sample_outs[i].write(var)





for out in sample_outs:
    out.close()
vcf_file.close()
vcf_out.close()








