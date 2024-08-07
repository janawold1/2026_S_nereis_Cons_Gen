#!/bin/bash -e
#SBATCH -A uc03718
#SBATCH -J Anahera_dorado_TGSGapClose
#SBATCH --partition=bigmem
#SBATCH --time=16:00:00
#SBATCH --cpus-per-task=32
#SBATCH --mem=160G
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=task

# TGS-GapCloser.sl
# Nat Forsdick, 2022-08-18
# Running TGS-GapCloser for kakī 
# Takes 2 params: 1: full path to assembly, and 2: output directory.

##########
# PARAMS #
TGSGapCloser=/home/jwo83/TGS-GapCloser/tgsgapcloser
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/

ml purge
ml load Racon/1.5.0-GCC-11.3.0
ml load minimap2/2.24-GCC-11.3.0

RACON=/opt/nesi/mahuika/Racon/1.5.0-GCC-11.3.0/bin/racon

# First need to convert HiFi fastq to fasta
#if [ ! -e ${HIFI}.fasta ]; then
#zcat $HIFIIN | sed -n '1~4s/^@/>/p;2~4p' > ${HIFI}.fasta
#fi

for indiv in Anahera
	do
#	mkdir ${dir}${indiv}/dorado_{q10,q20}_gapcloser
	echo "RUNNING GAPCLOSER FOR Q10 FLYE ASSEMBLY FOR ${indiv} AT "
	date
	cd ${dir}${indiv}/dorado_q10_gapcloser/
	$TGSGapCloser \
		--scaff assembly.fasta \
		--reads ${indiv}_dorado_q10_1kb.fasta \
		--output ${indiv}_dorado_q10_gapcloser \
		--minmap_arg '-x ava-ont' \
		--racon $RACON \
		--tgstype ont \
		--thread 32 >pipe.log 2>pipe.err
	echo "RUNNING GAPCLOSER FOR Q20 FLYE ASSEMBLY FOR ${indiv} AT "
	date
	cd ${dir}${indiv}/dorado_q20_gapcloser/
	$TGSGapCloser \
		--scaff ${dir}${indiv}/dorado_q20_flye/assembly.fasta \
		--reads ${dir}${indiv}/fasta/${indiv}_dorado_q20_1kb.fasta \
		--output ${indiv}_dorado_q20_gapcloser \
		--minmap_arg '-x ava-ont' \
		--racon $RACON \
		--tgstype ont \
		--thread 32 >pipe.log 2>pipe.err
	echo "FINISHED RUNNING GAPCLOSER FOR ${indiv} FLYE ASSEMBLIES AT "
	date
done
