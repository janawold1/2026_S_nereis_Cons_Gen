# Structural Variant discovery with a linear Graph
## CuteSV and Sniffles Calling
```
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

for indiv in Anahera Bill Blades Gulliver Huhu Mati-ma Richard_Henry Te_Ariki
	do
	echo "BEGINNING RUNNING CUTESV USING DORADO READS FOR ${indiv} AT "
	date
	cuteSV -t 16 -S ${indiv} -r 1000 -l 50 \
		--max_cluster_bias_INS 100 --diff_ratio_merging_INS 0.3 \
		--max_cluster_bias_DEL 100 --diff_ratio_merging_DEL 0.3 \
        --report_readid -s 3 \
		${dir}bam/${indiv}*.sorted.bam \
		${ref} \
		${dir}SVs/${indiv}_cuteSV.vcf \
		${dir}SVs/
	echo "FINISHED RUNNING CUTESV FOR $indiv AT "
	date
done

echo "FINISHED RUNNING ALL CUTESV SAMPLES... STARTING TO RUN SNIFFLES..."
ml purge
ml load Sniffles/2.0.7-gimkl-2022a-Python-3.10.5

for indiv in Anahera Bill Blades Gulliver Huhu Mati-ma Richard_Henry Te_Ariki
	do
	echo "BEGINNING RUNNING SNIFFLES USING DORADO READS FOR ${indiv} AT "
	date
	sniffles --input ${dir}bam/${indiv}*.sorted.bam \
		--vcf ${dir}SVs/${indiv}_sniffles.vcf --threads 16 --reference $ref \
		--sample-id ${indiv} --output-rnames --combine-consensus --allow-overwrite --minsvlen 50
	echo "FINISHED RUNNING SNIFFLES FOR DORADO READS AND STARTED GUPPY READS FOR ${indiv} AT "
	date
done

echo "FINISHED SNIFFLES SV DISCOVERY AT "
date
```

## DysGu
Need to figure out exactly how I ran dysgu.... Was super easy from memory....
```
dysgu code
```
## Merging with Jasmine
The tool [Jasmine v1.1.5](github.com/mkirsche/Jasmine) was used to merge calls across all individuals and tools as per below. For Jasmine, the input `file_list` corresponds to each individual's VCFs called with either Dysgu, CuteSV and Sniffles. the input `bam_list` must also correspond to each VCF denoted in the `file_list`.
```

```
After initial preprocessing with Jasmine, calls were refined with [Iris v1.0.4](github.com/mkirsche/Iris)

First had to make Sniffles calls compatible with Iris prior to call refinement. The below script was borrowed from [Laurie Lecomte](github.com/LaurieLecomte/SV_long_reads).  
```
#Extract old names
bcftools query -l sniffles/sniffles.vcf > sniffles.old_names

#Create new names
sed -E 's/[0-9]+\_([A-Za-z0-9]+)/\1/' sniffles.old_names > sniffles.new_names

#Put both into a single file
paste -d "\t" sniffles.old_names sniffles.new_names > sniffles.rename

bcftools annotate -x
```

## Filtering