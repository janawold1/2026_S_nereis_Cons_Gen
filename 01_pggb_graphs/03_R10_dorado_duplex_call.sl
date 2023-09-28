#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name R10_duplex
#SBATCH --partition=gpu,hgx
#SBATCH --gpus-per-node=A100:1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time 03:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Basecalling ONT simplex and duplex reads for genome assembly.
###############################################################

ml purge
ml load Dorado/0.2.1

# And now beginning to run the associated shell script.
pod5_dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/R10_raw_reads/
out=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_basecalled_reads/
split=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_split_pod5/
lib=lib1

for indiv in A B C
	do
	echo "BEGINNING STEREO BASECALLING OF MAIN DUPLEX READS FOR ${indiv} AT"
	date
	dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.1.0 \
		${pod5_dir}${indiv}/${lib}/ \
		--pairs ${out}${indiv}/${lib}_pair_ids_filtered.txt > ${out}${indiv}/${indiv}_${lib}_duplex_orig.sam
	echo "FINISHED STEREO BASECALLING OF MAIN DUPLEX READS FOR ${indiv} AT"
	date
	echo "BEGINNING STEREO BASECALLING OF SPLIT DUPLEX READS FOR ${indiv} AT"
	date
#	cat ${split}${indiv}/${lib}/*_split_duplex_pair_ids.txt > ${out}${indiv}/${lib}_split_duplex_pair_ids.txt
	dorado duplex dna_r10.4.1_e8.2_400bps_sup@v4.1.0 \
		${split}${indiv}/${lib}/ \
		--pairs ${out}${indiv}/${lib}_split_duplex_pair_ids.txt > ${out}${indiv}/${indiv}_${lib}_duplex_splitduplex.sam
	echo "FINISHED STEREO BASECALLING OF SPLIT DUPLEX READS FOR ${indiv} AT"
	date
done
