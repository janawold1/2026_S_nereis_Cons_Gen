#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name Te_Ariki_dorado_quast
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time 05:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load QUAST/5.2.0-gimkl-2022a

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/
ref=/scale_wlg_nobackup/filesets/nobackup/uc03718/Janes_genome/Jane_chr.fa.gz

for indiv in Te_Ariki
	do
	mkdir ${dir}${indiv}/dorado_q10_quast
	printf "\nRUNNING QUAST FOR DORADO ASSEMBLY AND LONGSTITCH FOR ${indiv} AT "
	date
	quast ${dir}${indiv}/dorado_q10_longstitch/${indiv}_dorado_q10_1kb_assembly.fa \
		${dir}${indiv}/dorado_q10_longstitch/${indiv}_dorado_q10_longstitch.fa \
		--output-dir ${dir}${indiv}/dorado_q10_quast -t 32 -r ${ref}
	printf "\nFINISHED RUNNING QUAST FOR ${indiv} AT "
	date
done
