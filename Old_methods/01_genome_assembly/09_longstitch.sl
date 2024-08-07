#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name Te_Ariki_dorado_longstitch
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time 16:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load LongStitch/1.0.4-Miniconda3

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/

for indiv in Te_Ariki
	do
#	mkdir ${dir}${indiv}/dorado_q10_longstitch
#	mkdir ${dir}${indiv}/dorado_q20_longstitch
#	cp ${dir}${indiv}/dorado_q10_flye/assembly.fasta ${dir}${indiv}/dorado_q10_longstitch/${indiv}_dorado_q10_1kb_assembly.fa
#	cp ${dir}${indiv}/fasta/${indiv}_dorado_q10_1kb.fq ${dir}${indiv}/dorado_q10_longstitch/
	cp ${dir}${indiv}/dorado_q20_flye/assembly.fasta ${dir}${indiv}/dorado_q20_longstitch/${indiv}_dorado_q20_1kb_assembly.fa
	cp ${dir}${indiv}/fasta/${indiv}_dorado_q20_1kb.fq ${dir}${indiv}/dorado_q20_longstitch/
#	printf "\nRUNNING LONGSTITCH FOR DORADO Q10 ASSEMBLY AND READS FOR ${indiv} AT\n"
#	date
#	longstitch run \
#		-C ${dir}${indiv}/dorado_q10_longstitch \
#		draft=${indiv}_dorado_q10_1kb_assembly \
#		reads=${indiv}_dorado_q10_1kb \
#		t=32 G=1.2g --debug
	printf "\nRUNNING LONGSTITCH FOR DORADO Q20 ASSEMBLY AND READS FOR ${indiv} AT\n"
	date
	longstitch run \
		-C ${dir}${indiv}/dorado_q20_longstitch \
		draft=${indiv}_dorado_q20_1kb_assembly \
		reads=${indiv}_dorado_q20_1kb \
		t=32 G=1.2g --debug
	printf "\nFINISHED RUNNING LONGSTITICH FOR ${indiv} AT"
	date
done
