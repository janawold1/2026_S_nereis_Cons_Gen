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
The `REF` allele for Sniffles output contained only N's... To fix this, we normalised the VCF.
```
bcftools norm --threads 16 -f $ref  -O z -o ${dir}SVs/sniffles/sniffles_SVs_norm.vcf.gz -c -s -D ${dir}SVs/sniffles/sniffles_SVs.vcf.gz
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

## Graph Construction with VG
Prior to graph construction, each chromosome from Jane's genome was extracted into individual `fasta` files and the naming convention for these chromosome was updated to `Jane#chromosome_${i}`, where i ranged from 1..25, Z, & W. This was to facilitate later graph construction with `minigraph` and `pggb`, and ensure read mapping was for individual chromosomes. Filtering of SV calls from `sniffles` and `cuteSV` was performed to retain only `PRECISE` SVs that had a minimum read support of 5x. To address a [1bp offset]() that occurrs in Sniffles calling, the VCF for normalised prior to filtering
```
bcftools norm \
	-f ${dir}Janes_genome/Jane_chr.fa \
	--threads 16 -O z -o ${dir}variants/SVs/sniffles/sniffles_SVs_norm.vcf.gz \
	-c s -D ${dir}variants/SVs/sniffles/sniffles_SVs.vcf.gz
bcftools view \
	-i 'INFO/PRECISE=1 & INFO/SUPPORT>=5' \
	-r NC_044283.2 -O v -o ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf \
	sniffles_SVs_norm.vcf

bcftools view \
	-i 'INFO/PRECISE=1 & INFO/RE>=5' \
	-r NC_044283.2 -O v -o ${dir}variants/SVs/cuteSV/cuteSV_chr7_PRECISE_minSUPP5.vcf \
	cuteSV.vcf
``` 
For consistency, the header of the VCF files for SV calls on chromosome 7 were augmented to remove all other chromosomes not included in the VCF and update the chromosome name annotation to `Jane#chromosome_7` to enable construction of graphs from these VCFs. Then VCFs were compressed with `bgzip` and indexed.  
```
bgzip ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf
tabix -p vcf ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf.gz

bgzip ${dir}variants/SVs/cuteSV/cuteSV_chr7_PRECISE_minSUPP5.vcf
tabix -p vcf ${dir}variants/SVs/cuteSV/cuteSV_chr7_PRECISE_minSUPP5.vcf.gz
```
**NOTE: VG DOES NOT ACCEPT INVERSION CALLS WITHOUT FULL SEQUENCE AND SVLEN=0**  
This is because if an inversion call has a length >0, it is assumed to be a complex variant (e.g., an inversion + INDEL). To address this in our trial for chromosome 7, only 1 inversion passed QC thresholds in the sniffles dataset. As first sanity check for inversion calls, the REF and ALT fields in the VCF were manually updated with the sequences as obtained below.  
```
REF=$(samtools faidx Janes_genome/Jane_chr.fa NC_044283.2:17819394-17821049 | sed 's/>NC_044283.2:17819394-17821049//g' | tr '\n' ' ' | sed 's/ //g')
INV=$(samtools faidx Janes_genome/Jane_chr.fa NC_044283.2:17819394-17821049 | sed 's/>NC_044283.2:17819394-17821049//g' | tr '\n' ' ' | sed 's/ //g' | rev)
```
Finally, graphs were constructed and indexed for two mapping tools: 1) `vg map`, a classic mapping tool suitable for both long- and short-reads; and 2) `vg giraffe` a mapping tool optimised for short-read data.  
#### Constructing graph for `vg giraffe`
Constructing and indexing a graph for `vg giraffe` can be done all in one step.  
```
vg autoindex \
	--workflow giraffe \
	-r ${dir}graphs/Jane_by_chromosome/Jane_chromosome_7.fa \
	-v ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf.gz \
	-p ${dir}graphs/sniffles_total_giraffe -t 16
```
And finally, mapping can be done by passing the required indices.  
```
vg giraffe \
	-Z ${dir}graphs/sniffles_total_giraffe.giraffe.gbz \
	-m ${dir}graphs/sniffles_total_giraffe.min \
	-d ${dir}graphs/sniffles_total_giraffe.dist \
	-f ${dir}graphs/illumina_by_chr/${indiv}_chr7_R1.fq \
	-f ${dir}graphs/illumina_by_chr/${indiv}_chr7_R2.fq > ${dir}graphs/giraffe_maps/${indiv}_sniffles_total_giraffe.gam
```
#### Constructing graph for `vg map`
First we constructed and indexed the graph.  
```
vg construct \
	-S -r ${dir}graphs/Jane_by_chromosome/Jane_chromosome_7.fa \
	-v ${dir}variants/SVs/sniffles/sniffles_chr7_PRECISE_minSUPP5.vcf.gz > ${dir}graphs/sniffles_total.vg

vg index -t 16 -x ${dir}graphs/sniffles_chr7.xg \
	-g ${dir}graphs/sniffles_chr7.gcsa \
	${dir}graphs/sniggles_chr7.vg
```
And finally, mapping was performed.  
```
vg map -t 16 -x ${dir}graphs/sniffles_total.xg \
	-g ${dir}graphs/sniffles_total.gcsa \
	-f ${dir}graphs/illumina_by_chr/${indiv}_chr7_R1.fq \
	-f ${dir}graphs/illumina_by_chr/${indiv}_chr7_R2.fq > ${dir}graphs/vg_maps/${indiv}_sniffles_chr7.gam
```
## Genotyping Graphs
### VG giraffe alignments
The linear graphs constructed from SV calls from Sniffles and CuteSV were genotyped. First, the output `.gbz` from `vg autoindex` was converted to `.xg` format for variant calls.
```
vg convert -x --drop-haplotypes \
	${dir}graphs/sniffles_giraffe_chr7.giraffe.gbz > ${dir}graphs/giraffe_maps/sniffles_giraffe_chr7.xg
```
Then individual SV calls were performed.  
```
vg pack --threads 16 -x ${dir}graphs/giraffe_maps/sniffles_giraffe_chr7.xg \
	-g ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.gam \
	-Q 5 -o ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.pack

vg call --genotype-snarls ${dir}graphs/giraffe_maps/sniffles_giraffe_chr7.xg \
	-k ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.pack > ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.vcf
```
### VG map alignments
As the indexing for these graphs is already available, no conversion step was needed.  
```
vg pack --threads 16 -x ${dir}graphs/sniffles_chr7.xg \
	-g ${dir}graphs/vg_maps/${indiv}_sniffles_chr7.gam \
	-Q 5 -o ${dir}graphs/vg_maps/${indiv}_sniffles_chr7.pack 

vg call --genotype-snarls ${dir}graphs/sniffles_chr7.xg \
	-k ${dir}graphs/vg_maps/${indiv}_sniffles_chr7.pack > ${dir}graphs/vg_maps/${indiv}_sniffles_chr7.vcf

```
## Mapping quality comparisons
The below is modified from [this](https://gtpb.github.io/CPANG18/pages/toy_examples) tutorial for assessing mapping scores of reads aligned to own assemblies and to the graph.

First loaded appropriate modules:
```
ml load vg/1.46.0
ml load jq/1.5
```
MapQ scores for reads aligned to genome graph were extracted with:
```
vg view -aj ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.gam | jq -cr '[.name, .mapping_quality] | @tsv' > ${dir}graphs/giraffe_maps/${indiv}_sniffles_chr7.tsv
```
To assess the quality of reads aligned to Jane:
```
samtools view -b ${dir}graphs/illumina_by_chr/${indiv}_chr7_Jane.bam | awk '{print $1"\t"$5}' > ${dir}graphs/illumina_by_chr/${indiv}_chr7_Jane_scores.tsv
```
Finally, the output of mapping quality for reads aligned to Jane's genome and the genome graph were merged with
```
join <(sort ${indiv}_self_score.tsv ) <(sort ${indiv}_vg_score.tsv ) | awk '{print $0"\t"$3 - $2}' | sort -n -k 4 > ${indiv}_chr7_align_self_compared.tsv
join <(sort ${indiv}_Jane_scores.tsv ) <(sort ${indiv}_vg_score.tsv ) | awk '{print $0"\t"$3 - $2}' | sort -n -k 4 > ${indiv}_chr7_align_comparisons.tsv
```
And evaluated at a high level with
```
cat ${indiv}_self_compared.tsv | awk '{ if ($4 < 0) print $1 }' | wc -l # Number of reads that aligned better to self
cat ${indiv}_self_compared.tsv | awk '{ if ($4 == 0) print $1 }' | wc -l # Number of alignments that were the same quality
cat ${indiv}_self_compared.tsv | awk '{ if ($4 > 0) print $1 }' | wc -l # Number of reads that aligned better in the graph

awk '{ if ($4 < 0) print $1 }' ${indiv}_chr7_align_comparisons.tsv | wc -l # Number of reads that aligned better to Jane
awk '{ if ($4 == 0) print $1 }' ${indiv}_chr7_align_comparisons.tsv | wc -l #Number of reads that had same alignment quality
awk '{ if ($4 > 0) print $1 }' ${indiv}_chr7_align_comparisons.tsv | wc -l #Number of reads that aligned better in the graph
```
### Visualising mapping quality
First augmented the input file
```
printf "indiv\tread_id\tscore\tdata\n" > giraffe_mapping_scores.tsv

awk '{print $1"\t"$2"\t"$3"\tJane_aligned"} *_chr7_align_comparisons.tsv >> giraffe_mapping_scores.tsv

awk '{print $1"\t"$2"\t"$4"\tgraph_aligned"} *_chr7_align_comparisons.tsv >> giraffe_mapping_scores.tsv

```
Setting up the R environment
```
library(ggplot2)
library(ggridges)
library(ggpubr)
library(gridExtra)
library(dplyr)
library(tidyverse)
library(Manu)

pal <- get_pal("Kakapo")
get_pal("Kakapo")
print_pal(pal)

knitr::opts_chunk$set(dev = c("svg", "png"),
                      dpi = 300,
                      echo = FALSE,
                      cache = TRUE)

setwd("C:/Users/Jana/Desktop/kakapo/linear_graphs")

LRmapQ <- read.table("vg_map/mapping_scores.tsv", sep = "\t", header = TRUE)
SRmapQ <- read.table("giraffe_mapping_scores.tsv", sep = "\t", header = TRUE)
```
Plotting mapping quality
```
pdf("LRmapQ_violin_plot.pdf")
ggplot(LRmapQ, aes(x = data, y = score, fill = data)) +
    geom_violin() +
    scale_fill_manual(values = c("Jane_aligned" = "#DCC949", "graph_aligned" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~indiv)
dev.off()

pdf("LRmapQ_density_plot.pdf")
ggplot(LRmapQ, aes(x = score, fill = data)) +
    geom_density() +
    scale_fill_manual(values = c("Jane_aligned" = "#DCC949", "graph_aligned" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~indiv)
dev.off()

pdf("Waihopai_SRmapQ_violin_plot.pdf")
SRmapQ %>% filter(across(indiv, ~grepl("Waihopai" , .))) %>%
ggplot(aes(x = data, y = score, fill = data)) +
    geom_violin() +
    scale_fill_manual(values = c("Jane_aligned" = "#DCC949", "graph_aligned" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~indiv)
dev.off()

pdf("SRmapQ_violin_plot.pdf")
SRmapQ %>%
ggplot(aes(x = data, y = score, fill = data)) +
    geom_violin() +
    scale_fill_manual(values = c("Jane_aligned" = "#DCC949", "graph_aligned" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~indiv)
dev.off()

pdf("SRmapQ_density_plot.pdf")
ggplot(SRmapQ, aes(x = score, fill = data)) +
    geom_density() +
    scale_fill_manual(values = c("Jane_aligned" = "#DCC949", "graph_aligned" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~indiv)
dev.off()
```
Plotting absolute number of read counts mapping 'better', 'same', or 'worse' in the graph.  
```
counts <- data.frame(indiv = c("Ariki", "Gertrude", "Waihopai", "JEM", "Ariki", "Gertrude", "Waihopai", "JEM","Ariki", "Gertrude", "Waihopai", "JEM"),
        	        quality = c("better", "better", "better", "better", "same", "same", "same", "same", "worse", "worse", "worse", "worse"),
                    counts = c(76665, 90827, 210023, 110129, 4129505, 4209521, 5589174, 7224515, 7563593, 10604047, 15556545, 7240617))

pdf("read_counts.pdf")
ggplot(counts, aes(x = quality, y = counts, fill = quality)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("worse" = "#CD8862", "same" = "#CED38C", "better" = "#7D9D33")) +
    facet_wrap(~indiv) +
    theme_light()
dev.off()
```
