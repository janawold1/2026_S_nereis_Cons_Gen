#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name dorado_duplex_pairs
#SBATCH --partition=milan
#SBATCH --cpus-per-task=16
#SBATCH --mem=30G
#SBATCH --time 04:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Basecalling ONT simplex and duplex reads for genome assembly.
###############################################################

ml purge
ml load duplex-tools/0.2.20-gimkl-2022a-Python-3.10.5
ml load SAMtools/1.16.1-GCC-11.3.0

# And now beginning to run the associated shell script.
pod5_dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/R10_raw_reads/
out=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_basecalled_reads/
dup_pod=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_split_pod5/
lib=lib1

for indiv in A B C
	do
#	echo Converting ${indiv} SAM to bam at
#	date
#	samtools view \
#		-@ 16 -b -o ${out}${indiv}/${indiv}_${lib}_R10_moves.bam \
#		${out}${indiv}/${indiv}_${lib}_R10_moves.sam
#	echo Finished converting SAM to BAM for ${indiv} at
#	date
#	echo Indexing ${indiv} SAM at
#	date
#	samtools index -@ 16 ${out}${indiv}/${indiv}_${lib}_R10_moves.bam
#	echo Finished ${indiv} indexing SAM at 
#	date
	cd ${out}${indiv}
	echo Finding ${indiv} main duplex pairs for Dorado stereo basecalling at
	date
	duplex_tools pair --output_dir ${out}${indiv}/ \
		--prefix ${indiv}_${lib}_R10 --verbose --threads 32 \
		${out}${indiv}/${indiv}_${lib}_R10_moves.bam
	echo Finished finding ${indiv} main duplex pairs for stereo basecalling at
	date
	echo Begin finding split duplex reads for $indiv at
	date
	duplex_tools split_pairs --debug --threads 32 \
		${out}${indiv}/${indiv}_${lib}_R10_moves.bam \
		${pod5_dir}${indiv}/${lib} \
		${dup_pod}${indiv}/${lib}
	echo Finished finding split duplex pairs for $indiv at
	date
done
