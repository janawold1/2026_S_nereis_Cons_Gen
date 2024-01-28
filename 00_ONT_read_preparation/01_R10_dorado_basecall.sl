#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name basecall_R10_flowcells
#SBATCH --partition=gpu,hgx
#SBATCH --gpus-per-node=A100:1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time 07:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Basecalling ONT simplex and duplex reads for genome assembly.
###############################################################

ml purge
ml load Dorado/0.4.2

# And now beginning to run the associated shell script.
pod5_dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/R10_raw_reads/
out=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_basecalled_reads/
lib=lib1
model=dna_r10.4.1_e8.2_400bps_sup@v4.2.0

dorado download --model ${model}

for indiv in A B C
        do
        echo "BEGINNING SIMPLEX BASECALLING FOR ${indiv} AT "
        date
        dorado basecaller \
                --min-qscore 10 --emit-moves ${model} \
                ${pod5_dir}${indiv}/${lib}/ > ${out}${indiv}/${indiv}_${lib}_R10_moves.sam --device 'cuda:all'
        echo "FINISHED SIMPLEX BASECALLING FOR ${indiv} AT "
        date
done