# SNP Analysis
## Filtering for Mendelian Inheritance
After calling 1,908,356 SNPs as per [06_SNP_genotyping.sl](github.com/nuzla/kakapo_pangenome/blob/main/01_genome_assemblies/-6_SNP_genotyping.sl), all SNPs with a minimum genotype quality of 10 were retained then tested for mendelian inheritance.
```
bcftools view -i 'GQ>=10' -O z -o bird_10minGQ.vcf.gz bird_clair3/merge_output.vcf.gz
```
We then extracted SNP calls for the sire and dam of each of our cell line individuals, then merged with our filtered SNPs for each individual before removing all sites with missing data.  
```
bcftools view -s Sire,Dam -O z -o bird_parents.vcf.gz kakapo125_pop_filter_snps_chrs.vcf.gz

bcftools merge -O z -o bird_family.vcf.gz bird_10minGQ.vcf.gz bird_parents.vcf.gz

bcftools view -i 'N_PASS(GT="mis")=0' -O z -o bird_family_noMiss.vcf.gz bird_family.vcf.gz
```
Finally, SNPs were tested for Mendelian Inheritance with BCFtools.  
```
bcftools +mendelian bird_family_noMiss.vcf.gz --mode + -t Sire,Dam,Child -O z -o bird_family_noMiss_mendelian.vcf.gz
```
This left 680,616 SNPs for bird A and 645,658 SNPs for bird B.  

To get a sense of SNP density on chromosome 7 (NC_044283.2), SNPs were merged for a final call set of XXX,XXX. This meant that XXX SNPs were missing for bird A, while XXX were missing for bird B. Files were merged, and SNPs for the sire and dam for the two cell line individuals was extracted from the kākāpō125+ SNP call set to ultimately plot SNP density was plotted in R with the `CMplot` package as per below.  
```
bcftools merge -O z -o 

```

```
tool=sniffles
printf "SNP\tChromosome\tPosition\tEnding\tSVtype\tGenotype\tTool\n" > SV_density.tsv
for tool in cuteSV sniffles
    do
    while read -r line
        do
        counter=$(echo $line | awk '{print $3}')
        end=$(echo $line | awk '{print $4}')
        snp=$(echo $line | awk '{print $1}')
        chr=$(echo $line | awk '{print $2}')
        type=$(echo $line | awk '{print $5}')
        geno=$(echo $line | awk '{print $6}')
        printf "${snp}:${pos}\t${chr}\t${counter}\t${end}\t${type}\t${geno}\t${tool}\n" >> SV_density.tsv
        while (($counter <= $end))
            do
            counter=$((counter + 1))
            printf "${snp}:${counter}\t${chr}\t${counter}\t${end}\t${type}\t${geno}\t${tool}\n" >> SV_density.tsv
        done
    done < "${tool}"_density.tsv
done
```