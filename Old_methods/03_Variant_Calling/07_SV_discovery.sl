#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name 06_SV_discovery
#SBATCH --partition=milan
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G
#SBATCH --time 04:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# ONT SV discovery.
###############################################################

ml purge
ml load cuteSV/2.0.2-gimkl-2020a-Python-3.8.2 

# And now beginning to run the associated shell script.
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/variants/
ref=/scale_wlg_nobackup/filesets/nobackup/uc03718/Janes_genome/GCF_004027225.2_bStrHab1.2.pri_genomic.fna

for indiv in Te_Ariki
	do
	echo "BEGINNING RUNNING CUTESV USING DORADO READS FOR ${indiv} AT "
	date
	cuteSV -t 16 -S ${indiv}_dorado_q15 -r 1000 -l 50 \
		--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
		--max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
		${dir}bam/${indiv}_dorado_q15.sorted.bam \
		${ref} \
		${dir}SVs/${indiv}_dorado_q15_cuteSV.vcf \
		${dir}SVs/
	echo "FINISHED RUNNING CUTESV FOR DORADO READS AND STARTED GUPPY READS FOR ${indiv} AT "
	date
	cuteSV -t 16 -S ${indiv}_guppy_q15 -r 1000 -l 50 \
		--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
		--max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
		${dir}bam/${indiv}_guppy_q15.sorted.bam \
		${ref} \
		${dir}SVs/${indiv}_guppy_q15_cuteSV.vcf \
		${dir}SVs/
	echo "FINISHED RUNNING CUTESV FOR $indiv AT "
	date
done

echo "FINISHED RUNNING ALL CUTESV SAMPLES... STARTING TO RUN SNIFFLES..."
ml purge
ml load Sniffles/2.0.7-gimkl-2022a-Python-3.10.5

for indiv in Anahera Bill Blades Gulliver Huhu Mati-ma Te_Ariki
	do
	echo "BEGINNING RUNNING SNIFFLES USING DORADO READS FOR ${indiv} AT "
	date
	sniffles --input ${dir}bam/${indiv}_dorado_q15.sorted.bam \
		--snf ${dir}SVs/${indiv}_dorado_q15.snf \
		 --minsvlen 50 --sample-id ${indiv}_dorado_q15
	echo "FINISHED RUNNING SNIFFLES FOR DORADO READS AND STARTED GUPPY READS FOR ${indiv} AT "
	date
	sniffles --input ${dir}bam/${indiv}_guppy_q15.sorted.bam \
		--snf ${dir}SVs/${indiv}_guppy_q15.snf \
		--minsvlen 50 --sample-id ${indiv}_guppy_q15
	echo "FINISHED RUNNING SNIFFLES FOR $indiv AT "
	date
done

sniffles --input ${dir}SVs/*_dorado_q15.snf --vcf ${dir}SVs/sniffles_dorado_SVs.vcf
sniffles --input ${dir}SVs/*_guppy_q15.snf --vcf ${dir}SVs/sniffles_guppy_SVs.vcf

echo "FINISHED SV DISCOVERY AT "
date