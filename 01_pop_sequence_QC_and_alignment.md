# Illumina Read Processing and Alignment
## Initial QC and Trimming
Initial sequencing results were visualised with FastQCv0.12.1 and MultiQCv1.13. Reads were trimmed using [TrimGalore](https://www.bioinformatics.babraham.ac.uk/projects/trim_galore/) v0.6.10.  
```
ref=/dir/to/reference.fasta
dir=/dir/to/reads/
reports=/dir/to/fastq_out/
trim_out=/dir/to/trimmmed_reads/

for samp in ${dir}*R1.fastq.gz
    do
    base=$(basename ${samp} _R1.fastq.gz)
    echo "Running Trim_galore for ${base}..."
    trim_galore --paired \
        --nextseq 28 \
        --2colour 20 \
        --cores 8 \
        --fastqc_args "--nogroup --outdir ${reports} --threads 32" \
        --length 50 \
        --output_dir ${trim_out} \
        --clip_R1 20 \
        --clip_R2 20 \
        --three_prime_clip_R1 5 \
        --three_prime_clip_R2 5 \
        --retain_unpaired \
        --length_1 55 \
        --length_2 55 \
        ${dir}${base}_R1.fastq.gz \
        ${dir}${base}_R2.fastq.gz
done
```
## Alignment
For population analyses, Illumina short-reads were aligned, PCR-duplicates removed, and alignment statistics estimated in the same manner for both fairy tern populations and kakī. All manipulation of alignment files was performed with [SAMtools](https://www.htslib.org/) v1.16.  
```
REF=/dir/to/reference.fasta
INPUT=/dir/to/reads/
OUTPUT=/dir/to/alignments/

for LIB in lib1 lib2 LIC001 LIC002
    do
    for SAMP in ${INPUT}${LIB}/*_R1.fq.gz
        do
        BASE=$(basename ${SAMP} _${LIB}_R1.fq.gz)
        INFOLINE=$(zcat ${SAMP} | head -n 1)
        INSTRUMENT=`echo ${INFOLINE} | cut -d ':' -f1`
        INSTRUMENTRUN=`echo ${INFOLINE} | cut -d ':' -f2`
        FLOWCELL=`echo ${INFOLINE} | cut -d ':' -f3`
        LANE=`echo ${INFOLINE} | cut -d ':' -f4`
        INDEX=`echo ${INFOLINE} | cut -d ':' -f10`

        RGID="ID:${INSTRUMENT}_${INSTRUMENTRUN}_${FLOWCELL}_${LANE}_${INDEX}"
        RGPL="PL:Illumina"
        RGPU="PU:${FLOWCELL}.${LANE}"
        RGLB="LB:${BASE}_${LIB}"
        RGSM="SM:${BASE}"

        printf "BEGAN ALIGNING READS FOR ID:${BASE}, LIB:${LIB} AT "
        date
        bwa mem -M -R @RG'\t'$RGID'\t'$RGPL'\t'$RGPU'\t'$RGLB'\t'$RGSM -t 46 $REF $SAMP ${INPUT}${LIB}/${BASE}_${LIB}_R2.fq.gz | samtools sort -@ 46 -o ${OUTPUT}${BASE}.bam
        printf "FINISHED ALIGNING READS FOR ID:${BASE}, LIB:${LIB} AT "
        date
    done
done
```
Files were then merged for each of the datasets, the reference sample (SP01) is the only individual sequenced at two different facilites.  
```
samtools merge -@ 16 \
    -o ${DIR}SP01_merged.bam \
    ${DIR}SP01_lib1.bam ${DIR}SP01_lib2.bam \
    ${DIR}SP01_LIC001.bam ${DIR}SP01_LIC002.bam
```
 
 This `for` loop below was used to merge all other samples. The sequence batches were run simultaneously by changing the suffix of the file (i.e., LIC002, lib2) in the initiation of the `for` loop and for the creation of the `$BASE` variable.  
```
for BAM in ${DIR}*_lib2.bam
    do
    BASE=$(basename $bam _lib2.bam)
    echo "MERGING BAMS FOR $BASE..."
    samtools merge -@64 -o ${DIR}${BASE}_merged.bam ${DIR}${BASE}*.bam
done
```
All alignments were then sorted and PCR duplicates marked and removed using `SAMtools` v.1.16.  
```
for BAM in ${DIR}*_merged.bam
    do
    BASE=$(basename $BAM _merged.bam)
    echo "NSORTING AND FIXING MATE PAIRS FOR ${BASE}\n"
    samtools sort -@64 -n ${BAM} | samtools fixmate -@64 -m -c - ${DIR}${BASE}.fixmate.bam
    samtools sort -@64 -o ${DIR}bam/${BASE}.fixmate.sorted.bam ${DIR}${BASE}.fixmate.bam
    printf "MARKING AND REMOVING PCR DUPLICATES FOR ${BASE} AT "
    date
    samtools markdup -r -@64 --write-index ${DIR}${BASE}.fixmate.sorted.bam ${DIR}${BASE}_nodup.bam
    printf "\nFINISHED MARKING AND REMOVING DUPLICATES FOR ${BASE} AT "
    date
done
```
Once BAMs were processed, autosomal chromosomes were extracted for population analyses. Mean mapping quality and alignment depth was then estimated using SAMtools and [mosdepth](https://github.com/brentp/mosdepth).  
```
cut -f1,2 SP01_5kb_ragtag.fa.fai | grep "CM020" | grep -v CM020462.1_RagTag | grep -v CM020463.1_RagTag | awk '{print $1"\t1\t"$2} > SP01_autosomes.bed

for SAMP in ${DIR}*_nodup.bam
    do
    BASE=$(basename $SAMP .bam)
    printf "EXTRACTING AUTOSOMES FOR ${BASE} AT "
    date
    samtools view -@64 -L SP01_genome/SP01_autosomes.bed \
        --write-index -o ${DIR}nodup/${BASE}_autosomes.bam ${SAMP}
    printf "FINISHED CONVERTING AND INDEXING FILES FOR ${BASE}, NOW ESTIMATING STATS AT "
    date
    samtools stats ${DIR}nodup/${BASE}_autosomes.bam > ${DIR}bam_stat/${BASE}_autosomes.stats
    mosdepth --threads 24 --fast-mode --by 50 ${DIR}bam_stat/${BASE}_autosomes ${DIR}nodup/${BASE}_autosomes.bam
    echo "FINISHED PROCESSING ${BASE}!"
done
```
For ease of comparisons, outputs were also plotted with:
```
multiqc --title "Fairy tern align stats" -o FT-alignments ${DIR}bam_stats/*
```

A similar process was performed for kakī alignments, with some minor differences. Because the kakī genome was assembled and scaffolded without the use of a chromosomally assembled genome, we identified scaffolds > 2Mb in length that demonstrated sequence coverage consistent with larger scaffolds. This indicated that Scaffolds 1-28 likely represent relatively well assembled chromosomes. These, exlucding the Z sex chromosome (Scaffold 4) were used in all downstream analyses.  