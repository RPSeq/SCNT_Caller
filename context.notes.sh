# Get trinucleotide contexts from SNP vcf using bedtools getfasta
bcftools view -f "PASS" delivs/calls/SNP.SCNT.vcf.gz -i "CASE=='SCNT'" \
    | bcftools query -f '%CHROM\t%POS\n' \
    | awk '{OFS="\t"} {print $1, $2-2, $2+1}' \
    | bedtools getfasta -fi ~/genomes/mouse/mm10/GRCm38_68.fa \
        -bed /dev/stdin -tab -fo /dev/stdout \
    | cut -f 2 \
    | paste - \
        <(bcftools view -f "PASS" delivs/calls/SNP.SCNT.vcf.gz -i "CASE=='SCNT'" \
            | bcftools query -f '%REF%ALT\t%CHROM\t%POS\t%EXPT\n') \
    | cat <(echo -e "CONTEXT\tTYPE\tCHROM\tPOS\tEXPT") - > MT.NRL-NT.df

# germline calls
# this is more stringent than the paper methods....
# using the MGP overlap hi-conf germline calls from FNR calculations
bcftools query -f '%CHROM\t%POS\n' rerun/gatk/gss/SCNT.ES.SNPS.BI.RECAL.GSS.vcf.gz \
    | awk '{OFS="\t"} {print $1, $2-2, $2+1}' \
    | bedtools getfasta -fi ~/genomes/mouse/mm10/GRCm38_68.fa \
        -bed /dev/stdin -tab -fo /dev/stdout \
    | cut -f 2 \
    | paste - \
        <(bcftools query -f '%REF%ALT\t%CHROM\t%POS\tGERM\n' \
            rerun/gatk/gss/SCNT.ES.SNPS.BI.RECAL.GSS.vcf.gz) \
    | cat <(echo -e "CONTEXT\tTYPE\tCHROM\tPOS\tEXPT") - > GERM.df



# try python from cmdline for strand normalizing
cat MT.NRL-NT.df | ./bin/snorm > MT.NRL-NT.norm.df
cat GERM.df | ./bin/snorm > GERM.norm.df