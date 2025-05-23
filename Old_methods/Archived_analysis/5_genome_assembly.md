# Kākāpō genome assembly for graph construction
Two assembly tools were initially trialled for the construction of kākāpō founder reference genomes. As MaSuRCA requires the mean and standard deviation for the insert length of untrimmed reads, we first aligned the unprocessed reads for each sequencing lane independently to Jane's genome as below.  

First aligned raw Illumina reads to estimate the mean read length and standard deviation for input into MaSuRCA config file with the below shell script.  This was repeated for each of the individuals used to build the graph.  

```
#!/bin/bash -e
ref=/scale_wlg_nobackup/filesets/nobackup/uc03718/Janes_genome/GCF_004027225.2_bStrHab1.2.pri_genomic.fna.gz #Reference genome for alignment
datadir=/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/ #Directory with fastq data
samdir=/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_alignments/sam #Sam file output
bamdir=/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_alignments/bam/ #Bam file output
fq1=_R1.fq.gz #Read 1 suffix
fq2=_R2.fq.gz #Read 2 suffix
platform="Illumina"

#First index the reference genome
time bwa index $ref

#Now, retrieving read group and instrument information.
for samp in ${datadir}*_R1.fq.gz #Remember to be explicit with file location
do
    base=$(basename ${samp} _R1.fq.gz)
    infoline=$(zcat ${samp} | head -n 1)
    instrument=`echo ${infoline} | cut -d ':' -f1`
    instrumentrun=`echo $infoline | cut -d ':' -f2`
    flowcell=`echo $infoline | cut -d ':' -f3`
    lane=`echo $infoline | cut -d ':' -f4`
    index=`echo $infoline | cut -d ':' -f10`

    #Now to incorporate this information into the alignment
    rgid="ID:${instrument}_${instrumentrun}_${flowcell}_${lane}_${index}"
    rgpl="PL:${platform}"
    rgpu="PU:${flowcell}.${lane}"
    rglb="LB:${base}_library1"
    rgsm="SM:${base}"

    echo "Aligning reads for $base" #Be explicit with file location for read 2 and the sam file output
    time bwa mem -M -R @RG'\t'$rgid'\t'$rgpl'\t'$rgpu'\t'$rglb'\t'$rgsm -t 64 $ref $samp ${datadir}${base}${fq2} > ${samdir}${base}.sam

    echo "Converting sam file to bam file for $base"
    time samtools view -T $ref -b ${samdir}${base}.sam > ${bamdir}${base}.bam
done
```

Insert mean and standard deviation were then estimated with:
```
ml load picard
ml load R

picard SortSam \
    -I ${datadir}sam/Bird_L001.sam \
    -O ${datadir}sam/Bird_L001.sorted.sam \
    -SO coordinate
picard CollectInsertSizeMetrics \
    -I ${datadir}sam/Bird_L001.sorted.sam \
    -O ${datadir}Bird_L001_insert_size_metrics.txt \
    -H ${datadir}Bird_L001_insert_size_histogram.pdf 
```

Then formatted config file as per:
```
DATA
PE = la 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L001_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L001_R2.fastq.gz
PE = lb 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L002_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L002_R2.fastq.gz
PE = lc 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L003_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L003_R2.fastq.gz
PE = ld 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L004_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L004_R2.fastq.gz
PE = le 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L005_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L005_R2.fastq.gz
PE = lf 369 84 /scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L006_R1.fastq.gz \
/scale_wlg_nobackup/filesets/nobackup/uc03718/illumina_reads/Bill_L006_R2.fastq.gz
NANOPORE = /scale_wlg_nobackup/filesets/nobackup/uc03718/ONT_trimmed_reads/Bill_q10_4kbtrim.fastq
END
PARAMETERS
EXTEND_JUMP_READS=0
GRAPH_KMER_SIZE = auto
USE_LINKING_MATES = 0
USE_GRID=1
GRID_ENGINE=SLURM
GRID_QUEUE=all
GRID_BATCH_SIZE=500000000
LHE_COVERAGE=25
LIMIT_JUMP_COVERAGE = 300
CA_PARAMETERS =  cgwErrorRate=0.15
CLOSE_GAPS=1
NUM_THREADS = 16
JF_SIZE = 12000000000
SOAP_ASSEMBLY=0
FLYE_ASSEMBLY=0
END
```

## MaSuRCA Slurm script

```
#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name Bill_masurca
#SBATCH --cpus-per-task=24 #28
#SBATCH --mem 16G #32
#SBATCH --time 08:00:00 #HH:MM:SS
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# MaSuRCA assemblies for kakapo genome graphs
# Jana Wold, 2023-18-01

## Modules
ml purge
ml load MaSuRCA/4.0.9-gimkl-2020a

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/

echo Beginning MaSuRCA with 24CPU and 16Gb mem at
date

bash ${dir}masurca/assemble.sh > ${dir}scripts/02_Bill_assemble.out

echo Finishing MaSuRCA at
date
```