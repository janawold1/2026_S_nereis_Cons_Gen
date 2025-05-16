# Genome Assembly
We used [HiFiasm](https://github.com/chhylp123/hifiasm) vX.X to generate *de novo* genome assemblies for our focal male and kākāpō hatched during the 2022 breeding season. Parental kmers were generated from Illumina resequence data using [yak]() vX.X as input into HiFiasm. Prior to generating kmer hash tables, parental reads were were trimmed and adapters removed using [Trimmomatic]() vX.X.  

| Individual |        Dam     |   Sire   |
|:----------:|:--------------:|:--------:|
|   Anahera  |        JEM     | Waihopai |
|  Madeline  | Margaret-Maree |  Sinbad  |
|   Ōtepoti  |        Waa     | Takitimu |
|  Te Ariki  |     Gertrude   |   Ariki  |

Reads were filtered further with fastp and Yak was then run as below.  
```
READS=/path/to/parental/reads/
OUT=/path/to/hifiasm/

fastp -5 -3 -n 0 -t 5 -T 5 -q 20 \
  -i ${READS}${INDIV}_R1_val_1.fq \
  -I ${READS}${INDIV}_R2_val_2.fq \
  -o ${OUT}${INDIV}_cleaned_R1.fq \
  -O ${OUT}${INDIV}_cleaned_R2.fq

printf "\nRUNNING YAK FOR $INDIV AT "
date
yak count -b37 -t32 -o ${OUT}${INDIV}_trimmed.yak \
        <(zcat ${READS}${INDIV}_R1_val_1.fq.gz) <(zcat ${READS}${INDIV}_R2_val_2.fq.gz) 
printf "\nFINISHED RUNNING YAK FOR $INDIV AT "
date
```

Four runs of HiFiasm were performed:
  1) With and parental kmers and default duplicate purging.  
  2) Without parental kmers and with duplicate purging by HiFiasm (`-l0`). This setting is recommended for inbred genomes.  
  3) With parental kmers and with `-l0`.
  4) Without Parental kmers and with `-l0`.  

Our HERRO corrected HiFi-like reads were used as input and reads >50 kb with a minimum Q-score of 10 were used as input for the 'ultralong' setting.  
```
hifiasm \
    -o ${out}hifiasm/${indiv}/${indiv}_trio_50kb.asm -t 32 \
    -1 mat.yak \
    -2 pat.yak \
    --ul ${reads}${indiv}/${indiv}_50kb_q10.fastq \
    ${reads}${indiv}/${indiv}_corrected_10kb.fasta

awk '/^S/{header=">"$2; for(i=4; i<=NF; i++) {header=header" "$i}; print header; printf "%s", $3 | "fold -w 80"; close("fold -w 80"); print ""}' ${out}hifiasm/${indiv}/${indiv}_trio_50kb.asm.bp.hap1.p_ctg.gfa > ${out}hifiasm/${indiv}/${indiv}_trio_50kb.asm.bp.hap1.p_ctg.fasta
```

## Initial Assembly QC
We used [Merqury](https://github.com/marbl/merqury) to assess haplotype assembly correctness, [BUSCO](https://busco.ezlab.org/busco_userguide.html) v5.4.7 to assess completeness and estimated summary stats using [QUAST](https://github.com/ablab/quast) v5.2.  

For assessments of haplotype correctness, we generated merqury hapmers for Parental short-reads and HERRO corrected ONT reads for our assembly.  
```
meryl count k=19 threads=32 memory=32 \
  ${DIR}trimmed_illumina_reads/mat_R*.fq \
  output ${DIR}merqury/mat.k19.meryl
meryl count k=19 threads=32 memory=32 \
  ${DIR}trimmed_illumina_reads/pat_R*.fq \
  output ${DIR}merqury/pat.k19.meryl
meryl count k=19 threads=32 memory=32 \
  ${DIR}dorado/offspring*_10kb.fasta \
  output ${DIR}merqury/offspring_ONT.k19.meryl

hapmers.sh \
  ${DIR}merqury/mat_k19.meryl \
  ${DIR}merqury/pat_k19.meryl \
  ${DIR}merqury/offspring_k19.meryl 

hapmers.sh \
  ${DIR}merqury/mat_k19.meryl \
  ${DIR}merqury/pat_k19.meryl \
  ${DIR}merqury/offspring_k19.meryl -no-filt

cd ${DIR}merqury/

merqury.sh \
  offspring_k19.meryl \
  mat_k19.hapmer.meryl \
  pat_k19.hapmer.meryl \
  ${DIR}hifiasm/offspring_hap1.p_ctg.fasta \
  ${DIR}hifiasm/offspring_hap2.p_ctg.fasta \
  offspring_k19_merqury
```

## Polishing and Scaffolding
We polished reads with [NextPolish2](https://github.com/Nextomics/NextPolish2).  
```
# make meryl database of assembly
meryl count k=15 output merylDB offspring_dual.fasta

# Get repetitive kmer list
meryl print greater-than distinct=0.9998 merylDB > repetitive_k15.txt

#use winnowmap to process reads (here using simplex reads for correction; will probably also test whether using the corrected reads improves things)
winnowmap -t 16 -W repetitive_k15.txt \
  -ax map-ont Dual_purged_R0048_plus_organelles.fasta \
  ../T2Tpolish/T2T-Polish/automated_polishing/R0048_dorado0.7_sup.q15_10kb.fastq | samtools sort -o Nextpolish2_R48test.sorted.bam
```
I then installed NextPolish2 in to a docker container with Apptainer (config file here), and generated 21- and 31-mer datasets for input into NextPolish2
```
yak count -o k21.yak -k 21 -b 37 <(zcat sr.R*.fastq.gz) <(zcat sr.R*.fastq.gz)

yak count -o k31.yak -k 31 -b 37 <(zcat sr.R*.fastq.gz) <(zcat sr.R*.fastq.gz)

nextPolish2 -r -t 4 Nextpolish2_R48test.sorted.bam \
  Dual_purged_R0048_plus_organelles.fasta k19.yak k21.yak k31.yak > Dual_purged_R0048_np2.fasta
```