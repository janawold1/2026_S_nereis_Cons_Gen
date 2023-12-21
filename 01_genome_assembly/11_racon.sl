#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name Te_Ariki_racon2
#SBATCH --cpus-per-task=30
#SBATCH --mem=90G
#SBATCH --time 8:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load minimap2/2.24-GCC-11.3.0
ml load Racon/1.5.0-GCC-11.3.0

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/

for indiv in Te_Ariki
	do
#	mkdir ${dir}${indiv}/dorado_q10_flye_racon
#	mkdir ${dir}${indiv}/dorado_q10_longstitch_racon
	printf "\nRUNNING RACON USING Q10 ASSEMBLIES FOR ${indiv} AT "
	date
#	minimap2 -t 64 -ax map-ont ${dir}${indiv}/dorado_q10_flye_racon/${indiv}_flye_racon1.fasta \
#		${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq > ${dir}${indiv}/dorado_q10_flye/${indiv}_q10_flye_racon1.sam
#	minimap2 -t 64 -ax map-ont ${dir}${indiv}/dorado_q10_longstitch_racon/${indiv}_q10_longstitch_racon1.fasta \
#		${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq > ${dir}${indiv}/dorado_q10_longstitch_racon/${indiv}_q10_longstitch_racon1.sam
	racon -t 30 ${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq \
		${dir}${indiv}/dorado_q10_flye_racon/${indiv}_q10_flye_racon1.sam \
		${dir}${indiv}/dorado_q10_flye_racon/${indiv}_flye_racon1.fasta > ${dir}${indiv}/dorado_q10_flye_racon/${indiv}_flye_racon2.fasta
	racon -t 30 ${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq \
		${dir}${indiv}/dorado_q10_longstitch_racon/${indiv}_q10_longstitch_racon1.sam \
		${dir}${indiv}/dorado_q10_longstitch_racon/${indiv}_q10_longstitch_racon1.fasta > ${dir}${indiv}/dorado_q10_longstitch_racon/${indiv}_q10_longstitch_racon2.fasta
	printf "\nFINISHED RUNNING RACON FOR Q10 ASSEMBLIES FOR ${indiv} AT "
	date
done
