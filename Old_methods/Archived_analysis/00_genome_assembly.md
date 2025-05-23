# Genome Graph Assembly
### Basecalling
XXX
```
model=dna_r10.4.1_e8.2_400bps_sup@5.0.0

dorado download --model $model 

for indiv in X Y Z
  do
  printf "STARTED BASECALLING FOR ${indiv} AT "
  date
  dorado basecaller $model ${pod5_dir}${indiv}/ > ${out}${indiv}/${indiv}_SUPv5.bam
  printf "FINISHED BASECALLING FOR ${indiv} AT "
  date
done
```
We then filtered for all simplex reads >=10kb in length for error correction using HERRO.  
```
for indiv in X Y Z
  do
  printf "EXTRACTING AND FILTERING READS GREATER THAN 10KB FOR ${indiv}..."
  samtools fastq -@32 ${out}${indiv}/${indiv}_SUPv5.bam > ${out}${indiv}/${indiv}_SUPv5.fastq
  seqkit seq -m 10000 ${out}${indiv}/${indiv}_SUPv5.fastq > ${out}${indiv}/${indiv}_10kb.fastq
  printf "STARTED READ CORRECTION FOR ${indiv} AT "
  date
  dorado correct -v ${out}${indiv}/${indiv}_10kb.fastq > ${out}${indiv}/${indiv}_10kb_corrected.fasta
  printf "FINISHED READ CORRECTION FOR ${indiv} AT "
  date
done
```
### Read filtering 
XXX

### kmer counts
To enable phasing of cell line assemblies, kmer counts for each individual's parents were estimated. Fortunately, Illumina resequence data was available for all for parents.  
#### Meryl
XXX
```
for indiv in Ariki Gertrude Waihopai JEM
    do
    printf "STARTED RUNNING MERYL FOR $indiv AT "
    date
    meryl count compress k=30 threads=32 memory=96 ${DIR}${indiv}_L00*.fq.gz output ${DIR}${indiv}_compress.k30.meryl
    printf "FINISHED RUNNING MERYL AT "
    date
```

## Create the hapmer DBs
```
$MERQURY/trio/hapmers.sh \
  maternal_compress.k30.meryl \
  paternal_compress.k30.meryl \
  child_compress.k30.meryl
```
#### yak
Yak, otherwise known as yet another kmer counter, 
```
for indiv in Ariki Gertrude Waihopai JEM
    do
    printf "STARTED RUNNING YAK FOR $indiv AT "
    date
    yak count -o ${dir}${indiv}_total.yak -b 37 -t 32 <(cat ${dir}${indiv}_R1.fq) <(cat ${dir}${indiv}_R2.fq)
    printf "FINISHED RUNNING YAK AT "
    date
done
```

## Genome Assembly
### Hifiasm
XXX

## Assembly 