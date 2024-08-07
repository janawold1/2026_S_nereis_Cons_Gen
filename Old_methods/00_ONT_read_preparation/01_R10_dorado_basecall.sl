#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name Katie_HERRO
#SBATCH --partition=gpu,hgx
#SBATCH --gpus-per-node=A100:1
#SBATCH --cpus-per-task=32
#SBATCH --mem=72G
#SBATCH --time 48:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Basecalling ONT simplex and duplex reads for genome assembly.
###############################################################

ml purge
ml load Dorado/0.7.0
ml load SeqKit/2.4.0

# And now beginning to run the associated shell script.
pod5_dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/R10_raw_reads/
out=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/
lib=PromethION_lib1
#model_Murph=dna_r10.4.1_e8.2_400bps_sup@v4.3.0
model=dna_r10.4.1_e8.2_400bps_sup@v5.0.0

#dorado download --model ${model}

for indiv in Katie
        do
        echo "BEGINNING BASECALLING FOR ${indiv} AT "
        date
        dorado basecaller $model ${pod5_dir}${indiv}/${lib}/ --resume-from ${out}${indiv}/${indiv}_${lib}_SUPv5.bam > ${out}${indiv}/${indiv}_SUPv5.bam 
        echo "FINISHED BASECALLING FOR ${indiv} AT "
        date
        ml purge; ml load SAMtools/1.19-GCC-12.3.0
        samtools fastq -@16 ${out}${indiv}/${indiv}_SUPv5.bam > ${out}${indiv}/${indiv}_SUPv5.fastq
        seqkit seq -m 10000 ${out}${indiv}/${indiv}_SUPv5.fastq > ${out}${indiv}/${indiv}_10kb.fastq
        printf "STARTED ERROR CORRECTION FOR $indiv AT "
        date
        dorado correct -v ${out}${indiv}/${indiv}_10kb.fastq > ${out}${indiv}/${indiv}_10kb_corrected.fasta
        printf "FINISHED ERROR CORRECTION FOR $indiv AT "
        date
done