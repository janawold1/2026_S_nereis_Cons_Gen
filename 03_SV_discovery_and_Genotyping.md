# Structural Variant Discovery and Analysis
Structural variant (SV) analyses followed a similar pipeline where variant discovery, variant quality filtering, and initial genotyping was performed under recommended pipelines. High quality genotypes representing homozygous reference, heterozygous and homozygous alternate calls were then plotted with SAMplot and evaluated using PlotCritic for call refinement. Calls passing evaluation with PlotCritic were then merged using SURVIVOR and used to construct a genome graph with the VG tools suite. Detailed methods for each of these steps using Delly, Manta and Smoove are outlined below.  

## Variant Discovery
### Manta
[Manta](https://github.com/Illumina/manta) v1.6.0 was used to call SVs for Australian fairy tern and tara iti. Three samples had to be excluded for Manta to run, AU13, TI06, TI34 & TI35. The errors indicated issues with the proportion of reads and read depth statistics. Each chromosome was called independently with the `--callRegions` flag to save computational resource and increase efficiency. Five samples did not pass [Manta's read-pair orientation threshold](https://github.com/Illumina/manta/issues/168) of 90% and were excluded from SV discovery. This included two Australian samples and three tara iti samples (AU08, AU13, SND11, TI35).

Running Manta is relatively simple, with the initial configuration setup for each individual chromosome and executed each run. Because Manta outputs inversions as breakends ([see here](https://github.com/Illumina/manta/blob/master/docs/userGuide/README.md)), Inversion calls were converted from Breakends using the `convertInversions.py` script supplied by Manta. It is best to do this step on the individual output files before merging into a single file.  
```
while read -r line
    do
    CHROM=$(echo $line | awk '{print $1}')
    REGION=$(echo $line | awk '{print $1":1-"$3}')
    echo $REGION
    configManta.py --referenceFasta $REF --region $REGION --runDir manta/${CHROM}/ --bam AU01_nodup_autosomes.bam \
        --bam AU03_nodup_autosomes.bam --bam AU04_nodup_autosomes.bam --bam AU06_nodup_autosomes.bam --bam AU09_nodup_autosomes.bam \
        --bam AU14_nodup_autosomes.bam --bam AU17_nodup_autosomes.bam --bam AU20_nodup_autosomes.bam --bam AU21_nodup_autosomes.bam \
        --bam AU23_nodup_autosomes.bam --bam AU24_nodup_autosomes.bam --bam AU25_nodup_autosomes.bam --bam AU27_nodup_autosomes.bam \
        --bam AU28_nodup_autosomes.bam --bam AU29_nodup_autosomes.bam --bam AU30_nodup_autosomes.bam --bam AU33_nodup_autosomes.bam \
        --bam SND04_nodup_autosomes.bam --bam SND05_nodup_autosomes.bam --bam SND06_nodup_autosomes.bam --bam SND15_nodup_autosomes.bam \
        --bam SP02_nodup_autosomes.bam --bam SP03_nodup_autosomes.bam --bam SP07_nodup_autosomes.bam --bam TI21_nodup_autosomes.bam \
        --bam TI36_nodup_autosomes.bam --bam TI37_nodup_autosomes.bam --bam TI38_nodup_autosomes.bam --bam TI40_nodup_autosomes.bam \
        --bam TI41_nodup_autosomes.bam
    manta/${CHROM}/runWorkflow.py
    convertInversion.py ~/anaconda3/envs/samtools/bin/samtools $REF \
        manta/${CHROM}/results/variants/diploidSV.vcf.gz | bgzip -c > manta/${CHROM}/results/variants/diploidINV.vcf.gz
done < reference/SP01_autosomes2.bed
``` 
Once individual chromosome runs were complete, VCFs were concatenated into a single file for SV quality filtering.  
```
echo CM0204{37..58}.1_RagTag | tr ' ' '\n' | while read chrom; do if [[ ! -d manta/${chrom} ]]; then continue; fi; find manta/${chrom}/results/variants/ -name diploidINV.vcf.gz | sort; done > manta/fairy_vcf_list

echo scaffold_{1..28} | tr ' ' '\n' | while read chrom; do if [[ ! -d manta/${chrom} ]]; then continue; fi; find manta/${chrom}/results/variants/ -name diploidINV.vcf.gz | sort; done > manta/kaki_vcf_list

bcftools concat --naive --file-list manta/fairy_vcf_list \
    -O z -o manta/01_fairy_SVcalls.vcf.gz

bcftools concat --naive --file-list manta/kaki_vcf_list \
    -O z -o manta/01_kaki_SVcalls.vcf.gz
```

## Filtering
For filtering raw Manta calls, all reads had to pass all 'hard' filtering thresholds and have 'PRECISE' breakpoints.  
```
bcftools view -i 'FILTER=="PASS" & IMPRECISE==0' \
    -O z -o manta/03_fairy_filtered.vcf.gz manta/02_fairy_SVcalls_INVconverted.vcf
bcftools index manta/03_fairy_filtered.vcf.gz

bcftools view -i 'FILTER=="PASS" & IMPRECISE==0' \
    -O z -o manta/03_kaki_filtered.vcf.gz manta/02_kaki_SVcalls_INVconverted.vcf
bcftools index manta/03_kaki_filtered.vcf.gz
```

## Genotyping
Then performed genotyping as per:
```
while read -r line
    do
    printf "STARTED GENOTYPING $line AT "
    date
    graphtyper genotype_sv reference/SP01_5kb_ragtag.fa manta/03_fairy_filtered.vcf.gz --sams=fairy.list --threads=16 --region=$line
done < fairy_chroms.tsv
```
Genotypes for individual chromosomes were then merged and filtered to include SVs that passed in the `AGGREGATED` genotyping model for graphtyper and had a mean GQ>=25.
```
echo CM0204{37..58}.1_RagTag | tr ' ' '\n' | while read chrom; do if [[ ! -d graphtyper/sv_results/${chrom} ]]; then continue; fi; find graphtyper/sv_results/${chrom} -name "*.vcf.gz" | sort; done > graphtyper/fairy_vcf_list

bcftools concat --naive --file-list graphtyper/fairy_vcf_list \
    -O z -o graphtyper/01_fairy_genotypes.vcf.gz

bcftools view -i 'SVMODEL="AGGREGATED" & mean(GQ)>=25 & N_PASS(GT="mis")=0 & FILTER="PASS" & AC<68 & SVTYPE!="INS"' \
    -O z -o graphtyper/02_fairy_genotypes_filtered_noINS.vcf.gz \
    graphtyper/01_fairy_genotypes.vcf.gz 

bcftools view -i 'AGGREGATED" & mean(GQ)>=25 & N_PASS(GT="mis")=0 & FILTER="PASS" & AC<68 & SVTYPE="INS" & INFO/SEQ!="."' \
    -O z -o graphtyper/02_fairy_genotypes_filtered_INS.vcf.gz \
    graphtyper/01_fairy_genotypes.vcf.gz

bcftools index grapthyper/02_fairy_genotypes_INS.vcf.gz
bcftools index grapthyper/02_fairy_genotypes_noINS.vcf.gz

bcftools concat -a -O z -o graphtyper/02_fairy_genotypes_filtered.vcf.gz 02_fairy_genotypes_filtered_noINS.vcf.gz 02_fairy_genotypes_filtered_INS.vcf.gz

rm graphtyper/02_fairy_genotypes_*
```
This genotyping and filtering approach was also used for kakī.  

|   SV Type   | Raw FT Genotypes | Filtered FT Genotypes | KĪ 50x Genotypes | KĪ 50x GenoFiltered | KĪ 10x Genotypes | KĪ 10x GenoFiltered |
|:-----------:|:----------------:|:---------------------:|:----------------:|:-------------------:|:----------------:|:-------------------:|
|  Deletion   |       5,948      |         1,087         |       6,912      |         533         |      6,299       |         308         |
| Duplication |        450       |           41          |         674      |          15         |       561        |           7         |
|  Insertion  |       3,659      |         1,544         |       5,304      |        1,192        |      4,431       |         724         |
|  Inversion  |        486       |           0           |         61       |          0          |        49        |          0          |
|    Total    |      10,543      |         2,672         |      12,951      |        1,740        |     11,340       |        1,039        |

### Smoove
SV calling
```
#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name initial_smoove_calls
#SBATCH --cpus-per-task=32
#SBATCH --mem=32G
#SBATCH --time 04:00:00  # Walltime (HH:MM:SS)
#SBATCH --output /nesi/nobackup/uc03718/smoove/slurms/TI_smoove_call_%A.%a.out
#SBATCH --error /nesi/nobackup/uc03718/smoove/slurms/TI_smoove_call__%A.%a.err
#SBATCH --profile=ALL

REF=/nesi/nobackup/uc03718/new_TI/hifiasm_longstitch_racon_ragtag/Katie_racon_ragtag_autosomes.fa
OUT=/nesi/nobackup/uc03718/smoove/calls/

ml purge
ml smoove/0.2.8-Miniconda3

SAMP=$(sed -n "$((SLURM_ARRAY_TASK_ID +1))p" /nesi/nobackup/uc03718/smoove/TI.list)
NAME=$(basename ${SAMP} _autosomes_nodup.bam)

smoove call --outdir ${OUT} --name ${NAME} --fasta ${REF} -p 1 --genotype ${SAMP}
```
Merging calls

```

```
Genotyping
```
#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name GLOBAL_smoove_genos
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time 04:00:00  # Walltime (HH:MM:SS)
#SBATCH --output /nesi/nobackup/uc03718/smoove/slurms/GLOBAL_smoove_genos_%A.%a.out
#SBATCH --error /nesi/nobackup/uc03718/smoove/slurms/GLOBAL_smoove_genos_%A.%a.err
#SBATCH --profile=ALL

REF=/nesi/nobackup/uc03718/new_TI/hifiasm_longstitch_racon_ragtag/Katie_racon_ragtag_autosomes.fa
OUT=/nesi/nobackup/uc03718/smoove/

ml purge
ml smoove/0.2.8-Miniconda3

SAMP=$(sed -n "$((SLURM_ARRAY_TASK_ID +1))p" /nesi/nobackup/uc03718/smoove/sample.list)
NAME=$(basename ${SAMP} _autosomes_nodup.bam)

smoove genotype -d -x -p 1 --name ${NAME}_GLOBALsites --outdir ${OUT}genotypes/ --fasta $REF --vcf ${OUT}00_GLOBAL_merged.sites.vcf.gz ${SAMP}
```
Filtering

```
bcftools view -i 'SVTYPE!="BND" & INFO/IMPRECISE!=1 & N_PASS(GT="mis")=0' 01_GLOBAL_genotypes.smoove.square.vcf.gz | \
    bcftools view -i '(SVTYPE = "DEL" & FMT/DHFFC[0-74] < 0.7) | (SVTYPE = "DUP" & FMT/DHBFC[0-74] > 1.3) | (SVTYPE = "INV")' -O b -o 02_GLOBAL_filtered.bcf
```

## SV Evaluation
Genotypes were visually inspected for support by plotting with SAMplot and scoring with [PlotCritic](https://github.com/jbelyeu/PlotCritic). As with Smeds et al. only one curator was used to prevent the overfiltering of genotypes. We used BCFtools to first extract regions of interest for each individual.
```
bcftools view -i 'SVTYPE=!"INS" & GT=="RR"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thomRef\n]' 02_fairy_genotypes_filtered.vcf.gz > fairy_samplot_regions.tsv
bcftools view -i 'SVTYPE=!"INS" & GT=="het"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thet\n]' 02_fairy_genotypes_filtered.vcf.gz >> fairy_samplot_regions.tsv
bcftools view -i 'SVTYPE=!"INS" & GT=="AA"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thomAlt\n]' 02_fairy_genotypes_filtered.vcf.gz >> fairy_samplot_regions.tsv

bcftools view -i 'SVTYPE=!"INS" & GT=="RR"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thomRef\n]' 02_kaki_genotypes_filtered.vcf.gz > kaki_samplot_regions.tsv
bcftools view -i 'SVTYPE=!"INS" & GT=="het"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thet\n]' 02_kaki_genotypes_filtered.vcf.gz >> kaki_samplot_regions.tsv
bcftools view -i 'SVTYPE=!"INS" & GT=="AA"' -f [%CHROM\t%POS\t%END\t%SVTYPE\t%SAMPLE\thomAlt\n]' 02_kaki_genotypes_filtered.vcf.gz >> kaki_samplot_regions.tsv
```

These regions were then plotted into individual folders with SAMPLOT.  
```
while read -r line
    do
    CHROM=$(echo $line | awk '{print $1})
    START=$(echo $line | awk '{print $2})
    ENDIN=$(echo $line | awk '{print $3})
    SVTYP=$(echo $line | awk '{print $4})
    SAMPL=$(echo $line | awk '{print $5})
    GENOS=$(echo $line | awk '{print $6})
    printf "STARTED SAMPLOT FOR $SAMPL ${CHROM}:${START} AT "
    date
    samplot plot -n ${SAMPL} -b alignments/${SAMPL}_nodup_autosomes.bam \
        -o fairy_SAMplot/${SVTYP}_${CHROM}_${START}_${ENDIN}_${SAMPLE}_${GENOS}.png \
        -c $CHROM \
        -s $START \
        -e $ENDIN \
        -t $SVTYP
done < fairy_samplot_regions.tsv
```
This was then repeated for kakī genotypes.  

We then created PlotCritic projects and manually scored all SV genotypes.  
```
plotcritic -p fairy_genos -i fairy_SAMplot/ -q "What is the Genotype? -A "a":"HomRef" "b":"Het" "c":"HomAlt"
plotcritic -p kaki_genos -i kaki_SAMplot/ -q "What is the Genotype? -A "a":"HomRef" "b":"Het" "c":"HomAlt"
```
Scored genotypes were then extracted from the PlotCritic Summary files upon completion of scoring.  
```
printf "Chromosome\tPosition\tEnd\tSV_type\tSample\tGraphTyper_GT\thomRef_score\thet_score\thomAlt_score\n" > {kaki,fairy}_scores.tsv
tail -n+4 fairy_genos/fairy_genos_summary_report.tsv | tr '_' '\t' | awk '{print $2"_"$3"\t"$4"\t"$5"\t"$1"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> fairy_scores.tsv

printf "Chromosome\tPosition\tEnd\tSV_Type\tSample\tGraphTyper_GT\tPlotCritic_GT\n" > fairy_genos.tsv
tail -n+2 fairy_scores.tsv | awk '{if ($7==100) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t0";
    else if ($8==100) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t1";
    else if ($9==100) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t2"}' >> fairy_genos.tsv

tail -n+4 kaki_genos/kaki_genos_summary_report.tsv | tr '_' '\t' | awk '{print $2"_"$3"\t"$4"\t"$5"\t"$1"\t$6"\t"$7"\t"$8"\t"$9"\t"$10}' >> kaki_scores.tsv
awk '{if ($7==100) print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t0";
    else if ($8==100) print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t1";
    else if ($9==100) print $1"\t"$3"\t"$3"\t"$4"\t"$5"\t"$6"\t2"}' kaki_scores.tsv >> kaki_genos.tsv
```
The number counts of each genotype for each individual were then acquired as below.  
```
tail -n+2 fairy_genos.tsv | cut -f5,7 | sort | uniq -c
```
## Population Structure and Differentiation
### MDS
We prepared the called SVs for ngsDist by using BCFtools to output the posterior probabilities (`FORMAT/PL`) for each individual's genotype.  
```
bcftools query -f '%CHROM\t%POS\t[%PL\t]\n' graphtyper/02_fairy_genotypes_filtered.vcf.gz | sed 's/,/\t/g' | gzip -c > graphtyper/fairy_mds/fairy_SVgenos.tsv.gz

NSITES=$(zcat graphtyper/fairy_mds/fairy_SVgenos.tsv.gz | wc -l)
ngsDist --geno graphtyper/fairy_mds/fairy_SVgenos.tsv.gz --probs --log_scale \
    -n_ind 34 -n_sites $NSITES -labels graphtyper/fairy_mds/mds_list -o graphtyper/fairy_mds/GLOBAL_svDist

tail -n +3 graphtyper/fairy_mds/GLOBAL_svDist | Rscript --vanilla --slave getMDS.R \
    --no_header --data_symm -n 4 -m "mds" -o graphtyper/fairy_mds/GLOBAL_sv.mds
```
### F<sub>ST</sub>
```
vcftools --gzvcf 02_fairy_genotypes.vcf.gz \
    --weir-fst-pop AU_samples.txt \
    --weir-fst-pop TI_samples.txt \
    --out fst/AU_TI
```

## Intersections with BUSCO genes
After calling BUSCO sites for the tara iti and kakī genome assemblies, we extracted regions corresponding with complete, single-copy BUSCOs and those regions corresponding to our autosomal chromosomes.
```
cut -f2,3,4,5 run_aves_odb10/full_table.tsv | grep Complete | awk '{print $2"\t"$3"\t"$4"\t"$5}' > SP01_scBUSCOs.bed

while read -r line
    do
    CHROM=(echo $line | awk '{print $1}')
    printf "EXTRACTING CHROMOSOME $CHROM...\n"
    grep "$CHROM" SP01_scBUSCOs.bed >> SP01_autosome_scBUSCOs.bed
done < SP01_autosomes2.bed

sort -k1,1 -k2,2n SP01_autosome_scBUSCOs.bed > SP01_autosome_scBUSCOs.sorted.bed
```
After extracting BUSCOs corresponding to our autosomal scaffolds, we identified SVs overlapping these genes by extracting SV regions and finding intersecting sites with BCFtools and BEDtools.  

```
bcftools query -f '%CHROM\t%POS\t%SVLEN\n' 02_fairy_genotypes_filtered.vcf.gz | awk '{print $1"\t"$2"\t"$2+$3}' | sort -k1,1 -k2,2n > fairy_sv.bed

bedtools intersect -wa -a fairy_sv.bed  -b SP01_autsome_scBUSCOs.sorted.bed | awk '{print $1"\t"$2}' > fairy_BUSCO_svs.txt

bcftools view -T fairy_BUSCO_svs.txt -O z -o 03_fairy_BUSCO.vcf.gz 02_fairy_genotypes_filtered.vcf.gz
```
## Site Frequency Spectrum
To generate a site frequency spectrum for SVs, we generated population specific VCFs for each group and attempted to polarise the VCFs with genotypes from reads used to polarise the SNP-base SFS. However, SVs were not reliably genotyped in either tern outgroup or the pied avocet for kakī. This made comparisons of allele counts intersecting BUSCO genes with SNP-based estimates of load challenging. However, given that the BUSCO regions were identified as single-copy and complete in our reference individuals, we assumed any SV intersecting with these regions would be under strong purifying selection and therefore is the 'derived' allele. To ensure comparisons between fairy terns was comparable, we excluded any SVs private to one group as below. For kakī, we repeated the process using the poaka.  

```
bcftools view -S TI_sampleIDs.txt 02_fairy_genotypes_filtered.vcf.gz | bcftools view -i 'GT=="alt"' -O z -o 03_TI_genotypes_filtered.vcf.gz
bcftools view -S AU_sampleIDs.txt 02_fairy_genotypes_filtered.vcf.gz | bcftools view -i 'GT=="alt"' -O z -o 03_AU_genotypes_filtered.vcf.gz

bcftools index 03_TI_genotypes_filtered.vcf.gz
bcftools index 03_AU_genotypes_filtered.vcf.gz

bcftools index 02_pied_genotypes_filtered.vcf.gz
bcftools index 02_kaki_genotypes_filtered.vcf.gz # high-coverage kakī genotypes
bcftools index 02_KI_10x_of_10x_SVs_filtered.vcf.gz # low-coverage _kakī genotypes

bcftools isec 03_TI_genotypes_filtered.vcf.gz 03_AU_gentoypes_filtered.vcf.gz -p fairy_intersection/
bcftools isec 02_kaki_genotypes_filtered.vcf.gz 02_pied_genotypes_filtered.vcf.gz -p stilt_intersection/
bcftools isec 02_KI_10x_of_10x_SVs_filtered.vcf.gz 02_pied_genotypes_filtered.vcf.gz -p stilt10x_intersection/
```
Finally, we created our individuals and population frequency files.  
```
bcftools query -T fairy_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tTI\tIntersecting\n]' fairy_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^fairy_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tTI\tNonintersecting\n]' fairy_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv

bcftools query -T fairy_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tAU\tIntersecting\n]' fairy_intersection/0003.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^fairy_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tAU\tNonintersecting\n]' fairy_intersection/0003.vcf >> indiv_harmful_SVallele_frequency.tsv

bcftools query -T kaki_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tIntersecting\n]' stilt_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^kaki_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tNonintersecting\n]' stilt_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv

bcftools query -T kaki_10x_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tIntersecting\n]' stilt10x_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
bcftools query -T ^kaki_10x_BUSCO_svs.txt -f '[%CHROM\t%POS\t%SVTYPE\t%SVLEN\t%SAMPLE\t%GT\tKI\tNonintersecting\n]' stilt10x_intersection/0002.vcf >> indiv_harmful_SVallele_frequency.tsv
```
## SV Profiling with Repeat Annotation
1) Extract SV START/END positions into bed and extract sequence information into fasta file by SVTYPE.  
```
count=1
bcftools query -f '%CHROM\t%POS\t%END\t%SVTYPE\n' 03_fairy_filteredSVs.vcf.gz | while read line
        do
        echo $line | awk -v var="$count" '{print $1"\t"$2-30"\t"$2+30"\t"$3-30"\t"$3+30"\t"$4"\t"var}' >> fairy_SV_intervals.txt
        ((count++))
done
```
We then divided these intervals into 'start' and 'end', and included variant and count information to identify specific variants associated with each interval.  
```
awk '{print $1"\t"$2"\t"$3"\t"$6"_"$7"_start"}' fairy_SV_intervals.txt > fairy_intervals.bed
awk '{print $1"\t"$4"\t"$5"\t"$6"_"$7"_end"}' fairy_SV_intervals.txt >> fairy_intervals.bed

bedtools getfasta -fi Katie_5kb_ragtag.fa -bed fairy_intervals.bed -fo fairy_intervals.fasta
```
Now we adjust the names of the relevant fasta files to relate to their specific SV IDs.  
```
awk '{print $1":"$2"-"$3"\t"$1":"$2"-"$3":"$4}' fairy_intervals.bed > fairy_rename.txt

seqkit replace --pattern '(.+)' --replacement '{kv}' --kv-file fairy_rename.txt fairy_intervals.fasta --keep-key > fairy_intervals_renamed.fasta
```
2) Get repeat sequences for query by converting GFF3 to BED from complex regions predicted by repeatmodeler
```
awk '{print $1"\t"$4"\t"$5"\t"$9}' tara-iti-repeat-masked/Katie_5kg_ragtag-renamed.complex_mask.gff3 | sed 's/Target=//g' > Katie_complex.bed
```
3) Rename chromosome names and create a file containing repeat motif ID
```
grep ">" tara-iti-repeat-masked/Katie_5kg_ragtag-renamed.simple_mask.soft.complex_mask_hard.fasta | sed 's/>//g' > seqkit.txt
grep ">" reference/Katie_5kb_ragtag.fasta | sed 's/>//g' > fasta.txt
paste seqkit.txt fasta.txt > fairy_chrom_name_conversion.tsv

awk '{print $1"_"$4"-"$5"#"$9}' tara-iti-repeat-masked/Katie_5kb_ragtag-renamed.complex_mask.gff3 | sed 's/Target=//g' | sed 's%/%%g' > fairy_gff_ids.txt

while read -r line
    do
    SEQID=$(echo $line | awk '{print $1}')
    CHROM=$(echo $line | awk '{print $2}')
    echo "CONVERTING $SEQID TO $CHROM..."
    sed -i "s/$SEQID/$CHROM/g" fairy_gff_ids.txt
    sed -i "s/$SEQID/$CHROM/g" Katie_complex.bed
done < fairy_chrom_name_conversion.tsv
```
To make sure the correct repeat IDs were aligned to the correct locations we checked as per:
```
awk '{FS="#"} {print $1} fairy_chrom_name_conversion.tsv | sed 's/RagTag_/RagTag:/g' | awk '{if ($1!=$2) print $0}'
```
Although there were overlaps that caused chomosome names to have numbers overlapping, we found that the bp coordinates matched between the two.  
4) Get fasta sequence interval for these repeats
```
bedtools getfasta -fi Katie_5kb_ragtag.fa -bed Katie_complex.bed -fo Katie_complex.fasta
```
5) Convert names in fasta file to contain repeat motifs
```
seqkit replace --pattern '(.+)' \
    --replacement '{kv}' \
    --kv-file fairy_rename_ids.txt Katie_complex.fasta \
    --keep-key > Katie_complex_repeat_names.fasta
```
Repeat masker has a maximum character limit of 50 for sequence names. `_RagTag` was excluded from the fasta file to run RepeatMasker.  
```
sed -i 's/_RagTag//g' Katie_complex_repeat_names.fasta
```
6) Query repeat sequences using SV breakend sequences with RepeatMasker
```
RepeatMasker -pa 8 -lib Katie_complex_repeat_names.fasta -dir ./ fairy_intervals_renamed.fasta
```
7) Extract hit locations and their repeat types for plotting.  

```
printf "Chromosome\tStart\tEnd\tSV Type\tVariant ID\tVariant Position\tRepeat\n" > fairy_repeat_hits.tsv

tail -n+4 fairy_intervals_renamed.fasta.out | \
    sed 's/[[:blank:]]\+/\t/g' | \
    cut -f 6,12 | \
    sed 's/Low_complexity/LowComplexity/g' | \
    sed 's/Simple_repeat/SimpleRepeat/g' | \
    sed 's/DNA?/DNA/g' | \
    sed -e 's/:\|-\|_/\t/g' | \
    sort -k1,1 -k2,2n | \
    awk '{FS="\t"} { if (NF>7) print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7; else print $0}' | \
    uniq >> fairy_repeat_hits.tsv
```
Now we need to identify which variants did not align to repeat elements. For this, we leverage the full `fairy_SV_intervals.tsv` file and the assigned variant IDs from step 1. We first extract the variant ID and positions representing SVs with positive hits into a new file.  
```
tail -n+2 fairy_repeat_hits.tsv | awk '{print $4"_"$5"_"$6}' > fairy_positive_hits.tsv
```
Now we extract those SVs that did not have hits.  
```
grep -v -f fairy_positive_hits.tsv fairy_intervals.bed | sed 's/_RagTag//g' | sed 's/_/\t/g' | awk '{print $0"\tNoRepeat"}'
```
After genotyping (as outlined below), we identified the relationship between SVs that were genotyped and the repeat landscape in the fairy tern and kakī genomes by extracting the start and ending positions from the relevant VCFs and `*_repeat_hits.tsv` files.  
```
bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 03_AU_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$2-30"\t"$2+30"\t"$4}' | \
    sed 's/_RagTag//g' > AU_genos_start.tsv
bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 03_AU_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$3-30"\t"$3+30"\t"$4}' | \
    sed 's/_RagTag//g' > AU_genos_end.tsv

bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 03_TI_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$2-30"\t"$2+30"\t"$4}' | \
    sed 's/_RagTag//g' > TI_genos_start.tsv
bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 03_TI_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$3-30"\t"$3+30"\t"$4}' | \
    sed 's/_RagTag//g' > TI_genos_end.tsv

bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 02_KI_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$2-30"\t"$2+30"\t"$4}' > KI_genos_start.tsv
bcftools query -f '%CHROM\tPOS\t%END\t%SVTYPE\n' 02_KI_genotypes_filtered.vcf.gz | \
    awk '{print $1"\t"$3-30"\t"$3+30"\t"$4}' > KI_genos_end.tsv
```
Now extracting the relevant repeat motifs.
```
grep -f AU_genos_start.tsv fairy_repeat_hits.tsv | grep -v end > fairy_geno_repeats.tsv
grep -f AU_genos_end.tsv fairy_repeat_hits.tsv | grep end >> fairy_geno_repeats.tsv
grep -f TI_genos_start.tsv fairy_repeat_hits.tsv | grep -v end >> fairy_geno_repeats.tsv
grep -f TI_genos_end.tsv fairy_repeat_hits.tsv | grep end >> fairy_geno_repeats.tsv

grep -f KI_genos_start.tsv kaki_repeat_hits.tsv | grep -v end > kaki_geno_repeats.tsv
grep -f KI_genos_start.tsv kaki_repeat_hits.tsv | grep end >> kaki_geno_repeats.tsv
```


### Association with Repeats
Here we leverage an approach implemented by Katarina Stewart in a manuscript under review with Molecular Ecology to identify SVs overlapping with repeats annotated in the tara iti and kakī genome assemblies. Please see [03_SV_discovery_and_PopGen.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/03_SV_discovery_and_PopGen.md) for how intersections between SVs and these repeats were identified.  

First we load the data frames and identify which SVs had repeats and which did not. Those that were not were then included in the analyses below. 
```
fairy_repeats = pd.read_csv('repeatmasker/tara_iti_all_repeat_summary_autosomes_simple.tsv', sep='\t')
fairy_repeats.drop_duplicates(subset=['Chromosome', 'Repeat'], inplace=True)
fairy_repeat_svs = pd.read_csv('repeatmasker/fairy_repeat_hits.tsv', sep='\t')

fairy_repeat_svs.drop_duplicates(subset=['Chromosome', 'Start', 'End', 'SV Type', 'Variant ID', 'Repeat'], inplace=True)
fairy_repeat_svs = fairy_repeat_svs[fairy_repeat_svs['Repeat']!='NoRepeat']
fairy_SVrepeat_counts = fairy_repeat_svs.groupby(['Chromosome', 'SV Type', 'Repeat']).size().reset_index(name='Count')

kaki_repeats = pd.read_csv('repeatmasker/kaki_all_repeat_summary_autosomes_simple.tsv', sep='\t')

kaki_repeat_svs = pd.read_csv('repeatmasker/kaki_repeat_hits_simple.tsv', sep='\t')
kaki_repeat_svs.drop_duplicates(subset=['Chromosome', 'Start', 'End', 'SV Type', 'Variant ID', 'Repeat'], inplace=True)

# Specify the order of kakī chromosomes
chromosome_order = ['scaffold_1', 'scaffold_2', 'scaffold_3', 'scaffold_5', 'scaffold_6', 'scaffold_7', 'scaffold_8', 'scaffold_9',
                    'scaffold_10', 'scaffold_11', 'scaffold_12', 'scaffold_13', 'scaffold_14', 'scaffold_15', 'scaffold_16', 'scaffold_17',
                    'scaffold_18', 'scaffold_19', 'scaffold_20', 'scaffold_21', 'scaffold_22', 'scaffold_23', 'scaffold_24', 'scaffold_25',
                    'scaffold_26', 'scaffold_27', 'scaffold_28']  # Custom order
kaki_repeats['Chromosome'] = pd.Categorical(kaki_repeats['Chromosome'], categories=chromosome_order, ordered=True)
kaki_repeat_svs['Chromosome'] = pd.Categorical(kaki_repeat_svs['Chromosome'], categories=chromosome_order, ordered=True)
kaki_SVrepeat_counts = kaki_repeat_svs.groupby(['Chromosome', 'SV Type', 'Repeat']).size().reset_index(name='Count')

pivot_kaki_repeats = kaki_repeats.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
pivot_fairy_repeats = fairy_repeats.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
```

```
pivot_kaki_repeats.head()
```

Here we count all repeat classes for each group and plot the type of repeat discovered on each assembled chromosome in the tara iti and kakī genome assemblies.  

```
repeat_color_map = {
    'DNA': '#E69F00', 
    'LINE': '#56B4E9', 
    'LTR': '#009E73', 
    'Low Complexity': '#F0E442', 
    'Rolling Circle': '#0072B2', 
    'SINE': '#D55E00', 
    'Satellite': '#CC79A7', 
    'Simple Repeat': '#000000', 
    'Unknown': '#999999', 
    'tRNA': '#44AA99'
    }

figs, axes = plt.subplots(1, 2, figsize=(20,5), sharex=False, sharey=False)

pivot_fairy_repeats.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in pivot_kaki_repeats.columns], ax=axes[0], legend=False)
axes[0].set_title('A)', fontsize=20, loc='left')
axes[0].set_xlabel('Chromosome', fontsize=18)
axes[0].set_ylabel('Count', fontsize=18)
axes[0].set_xticks([])
axes[0].tick_params(axis='y', labelsize=16)

pivot_kaki_repeats.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in pivot_kaki_repeats.columns], ax=axes[1])
axes[1].set_title('B)', fontsize=20, loc='left')
axes[1].set_xlabel('Chromosome', fontsize=18)
axes[1].set_xticks([])
axes[1].tick_params(axis='y', labelsize=16)

plt.savefig('plots/Figure_rep_content_partA.png', dpi=300, bbox_inches='tight')
```
And here we examine the number of repeats intersecting SVs filtered for variant quality. We specifically examine deletions (C & F), duplications (D & G) and inversions (E & H).  
```
fairy_del_counts = fairy_SVrepeat_counts[fairy_SVrepeat_counts['SV Type']=='DEL']
fairy_del_pivot = fairy_del_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
fairy_dup_counts = fairy_SVrepeat_counts[fairy_SVrepeat_counts['SV Type']=='DUP']
fairy_dup_pivot = fairy_dup_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
fairy_ins_counts = fairy_SVrepeat_counts[fairy_SVrepeat_counts['SV Type']=='INS']
fairy_ins_pivot = fairy_ins_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
fairy_inv_counts = fairy_SVrepeat_counts[fairy_SVrepeat_counts['SV Type']=='INV']
fairy_inv_pivot = fairy_inv_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)

kaki_del_counts = kaki_SVrepeat_counts[kaki_SVrepeat_counts['SV Type']=='DEL']
kaki_del_pivot = kaki_del_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
kaki_dup_counts = kaki_SVrepeat_counts[kaki_SVrepeat_counts['SV Type']=='DUP']
kaki_dup_pivot = kaki_dup_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
kaki_ins_counts = kaki_SVrepeat_counts[kaki_SVrepeat_counts['SV Type']=='INS']
kaki_ins_pivot = kaki_ins_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)
kaki_inv_counts = kaki_SVrepeat_counts[kaki_SVrepeat_counts['SV Type']=='INV']
kaki_inv_pivot = kaki_inv_counts.pivot(index='Chromosome', columns='Repeat', values='Count').fillna(0)

figs, axes = plt.subplots(2, 4, figsize=(20,10), sharex=False, sharey=False)

fairy_del_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in fairy_del_pivot.columns], ax=axes[0,0], legend=False)
axes[0,0].set_title('C)', fontsize=20, loc='left')
axes[0,0].set_xticks([])
axes[0,0].set_xlabel('')
axes[0,0].set_ylabel('Count', fontsize=18)
axes[0,0].tick_params(axis='y', labelsize=16)

fairy_dup_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in fairy_dup_pivot.columns], ax=axes[0,1], legend=False)
axes[0,1].set_title('D)', fontsize=20, loc='left')
axes[0,1].set_xticks([])
axes[0,1].set_xlabel('')
axes[0,1].tick_params(axis='y', labelsize=16)

fairy_ins_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in fairy_ins_pivot.columns], ax=axes[0,2], legend=False)
axes[0,2].set_title('E)', fontsize=20, loc='left')
axes[0,2].set_xticks([])
axes[0,2].set_xlabel('')
axes[0,2].tick_params(axis='y', labelsize=16)

fairy_inv_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in fairy_inv_pivot.columns], ax=axes[0,3], legend=False)
axes[0,3].set_title('F)', fontsize=20, loc='left')
axes[0,3].set_xticks([])
axes[0,3].set_xlabel('')
axes[0,3].tick_params(axis='y', labelsize=16)

kaki_del_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in kaki_del_pivot.columns], ax=axes[1,0], legend=False)
axes[1,0].set_title('G)', fontsize=20, loc='left')
axes[1,0].set_xticks([])
axes[1,0].set_ylabel('Count', fontsize=18)
axes[1,0].set_xlabel('Chromosome', fontsize=18)
axes[1,0].tick_params(axis='y', labelsize=16)

kaki_dup_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in kaki_dup_pivot.columns], ax=axes[1,1], legend=False)
axes[1,1].set_title('H)', fontsize=20, loc='left')
axes[1,1].set_xticks([])
axes[1,1].set_xlabel('Chromosome', fontsize=18)
axes[1,1].tick_params(axis='y', labelsize=16)

kaki_ins_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in kaki_ins_pivot.columns], ax=axes[1,2], legend=False)
axes[1,2].set_title('I)', fontsize=20, loc='left')
axes[1,2].set_xticks([])
axes[1,2].set_xlabel('Chromosome', fontsize=18)
axes[1,2].tick_params(axis='y', labelsize=16)

kaki_inv_pivot.plot(kind='bar', stacked=True, color=[repeat_color_map[repeat] for repeat in kaki_inv_pivot.columns], ax=axes[1,3], legend=False)
axes[1,3].set_title('J)', fontsize=20, loc='left')
axes[1,3].set_xticks([])
axes[1,3].set_xlabel('Chromosome', fontsize=18)
axes[1,3].tick_params(axis='y', labelsize=16)

plt.savefig('plots/Figure_rep_content_partB.png', dpi=300, bbox_inches='tight')

```