#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name Te_Ariki_FLYE
#SBATCH --cpus-per-task=32
#SBATCH --mem=100G
#SBATCH --time 8:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load Flye/2.9.1-gimkl-2022a-Python-3.10.5

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/

for indiv in Te_Ariki
	do
#	mkdir ${dir}${indiv}/dorado_q10_flye
#	mkdir ${dir}${indiv}/dorado_q20_flye
#	printf "\nRUNNING FLYE USING DORADO Q10 READS FOR ${indiv}...\n"
#	flye --nano-raw ${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq \
#		--out-dir ${dir}${indiv}/dorado_q10_flye/ --genome-size 1.2g \
#		--threads 24 --debug
	printf "\nRUNNING FLYE USING DORADO Q20 READS FOR ${indiv}...\n"
	flye --nano-raw ${dir}${indiv}/fasta/${indiv}_dorado_q20_1kb.fq \
		--out-dir ${dir}${indiv}/dorado_q20_flye/ --genome-size 1.2g \
		--threads 24 --debug --resume
done
