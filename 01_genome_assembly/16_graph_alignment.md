# Steps for graph manipulation and read alignment
Need to document how chromosome 7 scaffolds were extracted out of fasta files for graph construction.  
_
Plan:
    Step 1: Align reads to individual genome assembly assembled as Ref to generate BAM.
        Identify all reads that align to chromosome 7 and extract with `samtools fastq`
    Step 2: Align read subset to graph for chromosome 7
        Attempted to use GraphAligner. Is buggy and won't align to `*.vg` file despite allegedly being compatible. Something wrong with the converstion from PGGB `*.gfa` to VG?
        Indexing .vg file is proving problematic. Potentially too RAM intensive?
    Step 3: Visualise with SequenceTubeMaps
        Requires GAM to be aligned to .vg file used as input. Won't work with GAM aligned to .gfa, despite the file format being the only difference.
        Attempting to align with `vg map`. Having trouble getting through indexing step.  

# Preparing Chromosome 7 Graph for Visualisation
PGGB does not natively write the required *.vg file format. First converted the output GFA with VGvXX.  
```
vg convert -t 8 --gfa-in chr7.*.smooth.final.gfa --hash-out > chr7.hash
vg mod -X 256 -M 32 chr7.hash > chr7_mod.hash
vg prune chr7_mod.hash > chr7_mod_prune.hash
```
This graph had too many components for generating a .gcsa index and the edges were too long. The graph was modified and pruned to help generate the `.gcsa` index as recommended.  
```
vg index -b ./ -g chr7_mod_prune.hash.gcsa chr7_mod_prune.hash
vg index -b ./ -x chr7_mod_prune.hash.xg chr7_mod_prune.hash
```
Generate snarls for genotyping the graph.  
```
vg snarls chr7.hash > chr7.hash.snarls
```
Finally, reads were mapped to the graph. A band-width `-w` of 120,000 was trialled to accommodate long-read lengths, but didn't seem to work well. Stuck with the default binwidth.  
```
vg map --xg-name chr7_mod_prune.hash.xg -g chr7_mod_prune.hash.gcsa -f Bird.fq.gz chr7.hash > Bird.gam
```


# GENOTYPING
```
vg pack -x chr7_mod_prune.hash.xg -o Bird.pack -g Bird.gam -Q 30
vg call -t 16 -k bird.pack -r chr7.hash.snarls -s Bird -a -p Jane#chromosome_7 chr7.hash > Bird.vcf
```

SequenceTubeMap only allows you to upload files <5Mb in size. As an alternative for larger graphs, you are able to 'mount' larger files by editing the config.json file. To help narrow down the large graph to regions of interest, the largest variants were identified from the VCF output by PGGB as per:
```
grep -v "#" *final.Jane.vcf | awk '{print $3"\t"$2}' | tr ">" "\t" | awk '{print $2-$1"\t"$3}' | sort -n | tail
```

 # Read Mapping
 Reads for aligning to the graph for chromosome 7 were identified by first aligning to back to the scaffolded assemblies.  
 ```
minimap2 minimap_ref
minimap2 mashmap_ref
 ```
Sorted BAM files were then indexed and subsetted for chromosome 7.  
```
samtools index Bird_racon2_longstitch_mashmap_asRef.sorted.bam

samtools view -b Bird_racon2_longstitch_mashmap_asRef.sorted.bam NC_044283.2 > Bird_mashmap_asRef_chr7.bam
```
Then extracted reads were written to fastq file.
```
samtools fastq Bird_racon2_longstitch_mashmap_asRef_chr7.bam > Bird_racon2_longstitch_mashmap_asRef_chr7_reads.fq
```
This amounted to 107,857 reads for Anahera kākāriki assembled as ref with MashMap as implemented in D-GENIES. Reads were then aligned to the graph with GraphAligner as per below. First attempted to align to vg file conversion. But GraphAligner would not align to this file

# Comparing mapping quality
Followed [this](https://gtpb.github.io/CPANG18/pages/toy_examples) tutorial for assessing mapping scores of reads aligned to own assemblies and to the graph.

First loaded appropriate modules:
```
ml load vg/1.46.0
ml load jq/1.5
```

MapQ scores for reads aligned to genome graph were extracte into TSVs with:
```
vg view -aj bird.gam | jq -cr '[.name, .mapping_quality] | @tsv' > bird_vg_score.tsv
```
For reads aligned to Jane
```
samtools view -b bird_Jane.bam NC_044283.2 | awk '{print $1"\t"$5}' > bird_Jane_scores.tsv
```

Files for comparison were merged with
```
join <(sort bird_self_score.tsv ) <(sort bird_vg_score.tsv ) | awk '{print $0, $3 - $2}' | tr ' ' '\t' | sort -n -k 4 > bird_self_compared.tsv
join <(sort bird_Jane_scores.tsv ) <(sort bird_vg_score.tsv ) | awk '{print $0, $3 - $2}' | tr ' ' '\t' | sort -n -k 4 > bird_Jane_compared.tsv
```

And evaluated with:
```
cat bird_self_compared.tsv | awk '{ if ($4 < 0) print $1 }' | wc -l # Number of reads that aligned better to self
cat bird_self_compared.tsv | awk '{ if ($4 == 0) print $1 }' | wc -l # Number of alignments that were the same quality
cat bird_self_compared.tsv | awk '{ if ($4 > 0) print $1 }' | wc -l # Number of reads that aligned better in the graph

cat bird_Jane_compared.tsv | awk '{ if ($4 < 0) print $1 }' | wc -l # Number of reads that aligned better to Jane
cat bird_Jane_compared.tsv | awk '{ if ($4 == 0) print $1 }' | wc -l #Number of reads that had same alignment quality
cat bird_Jane_compared.tsv | awk '{ if ($4 > 0) print $1 }' | wc -l #Number of reads that aligned better in the graph
```

# R plots
There were more reads that aligned better to Jane than the graph, but the graph had a higher mapping Q score. To compare the mapping characteristics in each I used the mapQ scores and the number of reads with better (or worse) mapping scores in graph alignments.  
### Setting up the environment
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

setwd("C:/Users/Jana/Desktop/kakapo/graphs/kakapo_seg100kb_98perc_k79/vg_nesi")

LRmapQ <- read.table("vg_map/mapping_scores.tsv", sep = "\t", header = TRUE)
SRmapQ <- read.table("giraffe_mapping_scores.tsv", sep = "\t", header = TRUE)


counts <- data.frame(bird = c("Anahera", "Te_Ariki", "Anahera", "Te_Ariki", "Anahera", "Te_Ariki"),
        	        quality = c("better", "better", "same", "same", "worse", "worse"),
                    counts = c(912, 1247, 94619, 108313, 1914, 2341))

pdf("LRmapQ_violin_plot.pdf")
ggplot(LRmapQ, aes(x = data, y = score, fill = data)) +
    geom_violin() +
    scale_fill_manual(values = c("Jane_algn" = "#DCC949", "vg_algn" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~bird)
dev.off()

pdf("LRmapQ_density_plot.pdf")
ggplot(LRmapQ, aes(x = score, fill = data)) +
    geom_density() +
    scale_fill_manual(values = c("Jane_algn" = "#DCC949", "vg_algn" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~bird)
dev.off()

pdf("Waihopai_SRmapQ_violin_plot.pdf")
SRmapQ %>% filter(across(bird, ~grepl("Waihopai" , .))) %>%
ggplot(aes(x = data, y = score, fill = data)) +
    geom_violin() +
    scale_fill_manual(values = c("Jane_algn" = "#DCC949", "giraffe_algn" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~bird)
dev.off()

pdf("SRmapQ_density_plot.pdf")
ggplot(SRmapQ, aes(x = score, fill = data)) +
    geom_density() +
    scale_fill_manual(values = c("Jane_algn" = "#DCC949", "giraffe_algn" = "#7D9D33")) +
    theme_light() +
    facet_wrap(~bird)
dev.off()

pdf("read_counts.pdf")
ggplot(counts, aes(x = quality, y = counts, fill = quality)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    scale_fill_manual(values = c("worse" = "#CD8862", "same" = "#CED38C", "better" = "#7D9D33")) +
    facet_wrap(~bird) +
    theme_light()
dev.off()
```

# Comparing SNP precision
Tutorial used can be found [here](https://pggb.readthedocs.io/en/latest/rst/tutorials/small_variants_evaluation.html#variants-evaluation).  

Used rtf-tools, as below, to evaluate variant calling between the graph and kākāpō125+ SNP set for chromsome 7. First converted Jane's chromosome 7 to SDF format.   
```
mkdir graphs/kakapo_seg100kb_98perc_k79/variant_evaluation

rtg format -o graphs/kakapo_seg100kb_98perc_k79/variant_evaluation/Jane#chromosome_7.sdf Jane_chr7.fa
```

Comparisons against the kākāpō125+ SNP call set, and those successfully genotyped from ONT longread data were made as per:
```
bcftools view -G -O z -o SNPs/kakapo125_pop_filter_snps_chr7_nogenos.vcf.gz SNPs/kakapo125_pop_filter_snps_chr7.vcf.gz
tabix SNPs/kakapo125_pop_filter_snps_chr7_nogenos.vcf.gz
tabix

bcftools annotate --rename-chrs chr_rename \
    -O z -o graphs/kakapo_seg100kb_98perc_k79/variant_evaluation/pggb_chr7.vcf.gz \
    graphs/kakapo_seg100kb_98perc_k79/chr7.fa.gz.*.smooth.final.Jane.vcf.gz
tabix graphs/kakapo_seg100kb_98perc_k79/variant_evaluation/pggb_chr7.vcf.gz

rtg vcfeval -t graphs/kakapo_seg100kb_98perc_k79/variant_evaluation/Jane#chromosome_7.sdf \
        -b SNPs/kakapo125_pop_filter_snps_chr7.vcf.gz \
        -c graphs/kakapo_seg100kb_98perc_k79/chr7.fa.gz.*.smooth.final.Jane.vcf.gz \
        -T 8 \
        -o graphs/kakapo_seg100kb_98perc_k79/variant_evaluation/full_chr7_SNPset
```

# Short-read alignment
Creating index for giraffe.  
```
ml purge
ml load vg/1.46.0

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/graphs/racon_longstitch_mashmap_asRef/kakapo_seg100kb_98perc_k79/

vg autoindex -w giraffe -p ${dir}vg_giraffe/pggb_chr7_giraffe -g ${dir}*.gfa -T ./ -t 16

for fq in ${dir}parent_reads/*_R1.fq
        do
        indiv=$(basename $fq _R1.fq)
        printf "\nALIGNING CHROMOSOME 7 READS FOR ${indiv} AT "
        date
        vg giraffe -Z ${dir}vg_giraffe/pggb_chr7_giraffe.giraffe.gbz \
                -m ${dir}vg_giraffe/pggb_chr7_giraffe.min \
                -d ${dir}vg_giraffe/pggb_chr7_giraffe.dist \
                -f $fq -f ${dir}parent_reads/${indiv}_R2.fq > ${dir}${indiv}.gam
        printf "\nFINISHED ALIGNING CHROMOSOME 7 READS FOR ${indiv} AT "
        date
done
```

Read alignment for giraffe