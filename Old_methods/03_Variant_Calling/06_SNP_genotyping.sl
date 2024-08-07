#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name Te_Ariki_CLAIR3_SNP_genos
#SBATCH --partition=milan
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time 02:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# ONT SNP discovery.
###############################################################

ml purge
ml load Clair3/1.0.0-Miniconda3

# And now beginning to run the associated shell script.
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/variants/
ref=/scale_wlg_nobackup/filesets/nobackup/uc03718/Janes_genome/GCF_004027225.2_bStrHab1.2.pri_genomic.fna
deepSNPS=/scale_wlg_nobackup/filesets/nobackup/uc03718/variants/SNPs/kakapo125_pop_filter_snps.vcf.gz
model=/scale_wlg_nobackup/filesets/nobackup/uc03718/scripts/r1041_e82_400bps_sup_v400

for indiv in Te_Ariki
	do
	echo "BEGINNING SNP DISCOVERY OF DORADO READS FOR ${indiv} AT "
	date
	run_clair3.sh \
		--bam_fn=${dir}bam/${indiv}_dorado_q20.sorted.bam \
		--ref_fn=${ref} \
		--sample_name=${indiv} \
		--print_ref_calls \
		--vcf_fn=${deepSNPS} \
		--threads=32 \
		--platform="ont" \
		--model_path=${model} \
		--output=${dir}SNPs/${indiv}_dorado_q20_R10_model/
	echo "FINISHED SNP DISCOVERY FOR $indiv AT "
	date
done
