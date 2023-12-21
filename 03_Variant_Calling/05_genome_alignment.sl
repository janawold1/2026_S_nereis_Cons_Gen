#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name Te_Ariki_Variant_alignments
#SBATCH --partition=milan
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time 16:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Basecalling ONT simplex and duplex reads for genome assembly.
###############################################################

ml purge
ml load minimap2/2.24-GCC-11.3.0
ml load SAMtools/1.16.1-GCC-11.3.0

# And now beginning to run the associated shell script.
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/variants/
ref=/scale_wlg_nobackup/filesets/nobackup/uc03718/Janes_genome/GCF_004027225.2_bStrHab1.2.pri_genomic.fna

for indiv in Te_Ariki
	do
	echo "BEGIN ALIGNING SNP READS FOR ${indiv} AT "
	date
	minimap2 -ax map-ont ${ref} ${dir}reads/${indiv}_dorado_q20.fq | samtools view -@32 -b > ${dir}bam/${indiv}_dorado_q20.bam
	samtools sort -@32 ${dir}bam/${indiv}_dorado_q20.bam > ${dir}bam/${indiv}_dorado_q20.sorted.bam
	echo "FINSIHED ALIGNING SNP READS FOR ${indiv} AT "
	date
	echo "BEGAN ALIGNING SV READS FOR ${indiv} AT "
	date
	minimap2 -ax map-ont ${ref} ${dir}reads/${indiv}_dorado_q15_1kb.fq | samtools view -@32 -b > ${dir}bam/${indiv}_dorado_q15_1kb.bam
	wait
	samtools sort -@32 ${dir}bam/${indiv}_dorado_q15_1kb.bam > ${dir}bam/${indiv}_dorado_q15_1kb.sorted.bam
	echo "ESTIMATING COVERAGE FOR ${indiv} AT "
	date
	samtools coverage ${dir}bam/${indiv}_dorado_q15_1kb.sorted.bam > ${dir}bam/${indiv}_dorado_q15_1kb_coverage.tsv 
	samtools coverage ${dir}bam/${indiv}_dorado_q20.sorted.bam > ${dir}bam/${indiv}_dorado_q20_coverage.tsv
	echo "FINISHED ALL ALIGNMENTS AND SORTING BAMS FOR $indiv AT "
	date
done
