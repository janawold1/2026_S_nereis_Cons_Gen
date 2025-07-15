# Population genomics with ANGSD
[ANGSDv0.935](https://bioconda.github.io/recipes/angsd/README.html?highlight=angsd) and ngsTools were used to estimate summary statistics for Australian fairy tern (*Sternual nereis nereis*) sampled from Western Australia and tara iti (*S. nereis davisae*) from Northland, NZ. The methods implemented below were modified from a helpful wiki provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md).  

Before progressing with analyses, an initial look at the distribution of quality scores and per-site depth on a global and individual-based basis. PCR duplicates were marked and sex chromosomes were excluded for SNP-based population analyses.  
```
ls ${DIR}{AU,SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}GLOBAL.list
ls ${DIR}{SND,SP,TI}*_markdup_autosomes.bam > ${ANGSD}TI.list
ls ${DIR}AU*_markdup_autosomes.bam > ${ANGSD}AU.list
ls ${DIR}H0*_markdup_autosomes.bam > ${ANGSD}KI.list

TREF=${DIR}SP01_genome/SP01_5kb_ragtag.fa
KREF=${DIR}kaki_genome/himNova-hic-scaff.fa

for POP in GLOBAL AU TI KI_10x
    do
    if (( $POP == KI_10x ))
    then
        angsd -P 32 -b ${ANGSD}${POP}.list -ref $KREF -out ${ANGSD}qc/${POP}.qc \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
    else
        angsd -P 32 -b ${ANGSD}${POP}.list -ref $REF -out ${ANGSD}qc/${POP}.qc \
            -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
            -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 800
    fi
done

angsd -P 16 -b ${ANGSD}KI.list -ref $KREF -out ${ANGSD}qc/KI.qc \
        -uniqueOnly 1 -remove_bads 1 -only_proper_pairs 1 -trim 0 -C 50 -baq 1 \
        -minMapQ 20 -doQsdist 1 -doDepth 1 -doCounts 1 -maxDepth 2000
```
Courtesy of scripts provided by [@mfumagalli](https://github.com/mfumagalli/ngsTools/blob/master/TUTORIAL.md), the distributions of these sites were visualised with `Rscript ~/ngsTools/Scripts/plotQC.R GLOBAL.qc` and meaningful filtering thresholds were identified. Below is a table outlining filtering thresholds for subsequent population analyses.

|        Population        |Minimum MapQ|Minimum Q|Minimum Depth|Maximum Depth|Number of Individuals|
| ------------------------ | ---------- | ------- | ----------- | ----------- | ------------------- |
|Australian Fairy Tern (WA)|     20     |    20   |     200     |     350     |         19          |
|         Tara iti         |     20     |    20   |     234     |     770     |         39          |
|     Global fairy tern    |     20     |    20   |     438     |    2,960    |         73          |
|         Kakī 10x         |     20     |    20   |     192     |     300     |         24          |
|         Kakī 50x         |     20     |    20   |     700     |    1,200    |         24          |

 Below is code run for analyses performed on alignments of fairy terns to the tara iti reference scaffolded using the common tern assembly (see [00_genome_assembly.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/00_genome_assembly.md)) and alignments of kakī to the high-quality reference assembly for this species (see [here](https://www.genomics-aotearoa.org.nz/our-work/completed-projects/high-quality-genomes) for details). All analyses for each group were performed in a similar manner for downstream comparisons. Quality thresholds were adjusted as per the table above.  

In addition, scaffolds that adhered to expected population sequence coverage (e.g., ~10x) and likely represented autosomal chromosomes were extracted from the reference fasta and individual BAM alignments before identifying putative coding regions as per below.  

## Excluding putative coding Regions
Neither the tara iti or kakī reference assemblies have transcriptomes. To assess the potential consequences of the inclusion of regions under selection, and to have a 'first look' at the genomic diversity around coding regions, we performed *ab initio* gene prediction using [AUGUSTUS](https://github.com/Gaius-Augustus/Augustus/tree/master) v3.5.0.  
```
augustus --sample=100 --alternatives-from-sampling=true --strand=both \
        --gff3=yes --progress=true --temperature=3 --species=chicken \
        --protein=on $TREF > ${DIR}augustus/TI_AUGUSTUS.gff

augustus --sample=100 --alternatives-from-sampling=true --strand=both \
        --gff3=yes --progress=true --temperature=3 --species=chicken \
        --protein=on $KREF > ${DIR}augustus/kaki_AUGUSTUS.gff
```
[BEDtools](https://bedtools.readthedocs.io/en/latest/content/tools/complement.html?highlight=complement) v2.31.1 was then used to sort, and merge these putative gene regions. An additional 1kb of sequence on either side of annotations were included to reduce linkage. Duplicate regions were then merged with `bedtools merge`.  
```
grep -v "#" reference/gene_predictions/SP01_AUGUSTUS.gff | \
    bedtools sort -i - > angsd/TI_augustus_autosomal_predictions.bed

bedtools merge -i angsd/TI_augustus_autosomal_predictions.bed > angsd/TI_augustus_autosomal_predictions_merged.bed

grep -v "#" kaki_genome/kaki_AUGUSTUS.gff | \
    grep -v scaffold_4 | \
    bedtools sort -i - > angsd/kaki_augustus_predictions.bed
```
 To account for linkage, the window size for these regions increased by 1kb on either side. The file was adusted in cases where the addition of a 1kb buffer extended beyond the start (n = 11) or end of the chromosome (n = 7).
```
awk '{print $1"\t"$2-1000"\t"$3+1000}' angsd/augustus_autosomal_predictions_merged.bed > angsd/augustus_autosomal_predictions_merged_add1kb.bed
awk '{print $1"\t"$2-1000"\t"$3+1000}' angsd/kaki_augustus_predictions_merged.bed > angsd/kaki_augustus_predictions_merged_add1kb.bed
```
To define putatively neutral sites for analyses, we extracted the complement regions as below.  
```
awk -v OFS='\t' '{print $1,$2}' ${TREF}$.fai > reference/SP01_autosome_lengths.txt

bedtools complement \
    -i angsd/augustus_autosomal_predictions_merged_add1kb.bed \
    -g reference/SP01_autosome_lengths.txt > angsd/SP01_neutral_sites.bed
```
For the fairy tern, this left 27,337 regions covering 562,663,009 bp for analyses, which roughly corresponds to 47% of the genome.  

This file was then indexed for ANGSD with `angsd sites`.  
```
angsd sites index TI_scaffolded_neutral_regions.bed
```