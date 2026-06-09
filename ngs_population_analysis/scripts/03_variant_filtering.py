import subprocess

hwe = 'vcftools --gzvcf chr16.vcf.gz --hwe 0.05 --min-alleles 2 --max-alleles 2 --recode --stdout | gzip -c > chr16_hwe.vcf.gz'
subprocess.run(hwe, shell = True)

het = 'vcftools --gzvcf chr16_hwe.vcf.gz --het'
subprocess.run(het, shell = True)

allele = 'vcftools --gzvcf chr16_hwe.vcf.gz --freq'
subprocess.run(allele, shell = True)

with open('out.het','r') as file:
    OHOM = []
    N_sites = []
    for line in file:
        readlist = line.split('\t')
        OHOM.append(readlist[1])
        N_sites.append(readlist[3])

OHOM.pop(0)
N_sites.pop(0)

het = []

for i in range(0,len(OHOM)):
    hom = int(OHOM[i])/int(N_sites[i])
    HET = 1-hom
    het.append(HET)

print('Heterozygosity:')
print(het)




