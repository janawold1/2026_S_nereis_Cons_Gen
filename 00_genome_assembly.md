# Tara iti Genome Assembly
Two tara iti clutchmates (1 male, 1 female) passed away as a result of a storm in the 2018/2019 breeding season. The female chick was sequenced with both short- and long-read technology. Short-reads, used for assembly polishing and assessing sequencing batch-effects, were generated using an Illumina NovaSeq 6000 at the [Institute of Clinical Molecular Biology](https://www.ikmb.uni-kiel.de/resources/sequencing/whole-genome-de-novo-sequencing) (Kiel, Germany) and Livestock Improvement Corp (Hamilton, NZ). Long-read sequencing was performed using a ONT PromethION (R10 flowcells, 5000Hz, LSK114 kit) at [Bragato Research Institute](https://bri.co.nz/) (Lincoln, NZ).  

The short-read data were processed as per [01_pop_sequence_QC_and_alignment.md](https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/01_pop_sequence_QC_and_alignment.md), while the ONT data was prepared as per below.

## Basecalling and Read Trimming
One library preparation was sufficient for all ONT sequencing. The PromethION flow cell was washed and loaded four times. All ONT `.pod5` files were basecalled using [Dorado](https://github.com/nanoporetech/dorado) v0.5.0 and the SUP model. First, simplex basecalling was performed, reads were then extracted to a fastq file with [SAMtools v1.19](), then trimmed with [SeqKit v2.4](https://github.com/shenwei356/seqkit) for a minimum length of 10kb and a minimum Q-score of 10 for read error correction (i.e,. [HERRO](https://github.com/lbcb-sci/herro)).  
```
dorado basecaller sup ${POD5}${indiv}/${lib}/ > ${BASECALL}Katie/Katie_SUPv5.bam

samtools fastq ${BASECALL}Katie/Katie_SUPv5.fastq

seqkit seq -m 10000 -Q 10 ${BASECALL}Katie/Katie_SUPv5.fastq > ${BASECALL}Katie/Katie_10kb_q10.fastq

dorado correct -v ${BASECALL}Katie/Katie_10kb_q10.fastq > ${BASECALL}Katie/Katie_corrected.fasta
```
Simplex reads were then filtered for a minimum length of 50kb for input into the ultralong read option of [HiFiasm v0.19.7-r598](https://github.com/chhylp123/hifiasm). HERRO corrected reads were also filtered for a minimum length of 10kb.
```
seqkit seq -m 50000 -Q 10 ${BASECALL}Katie/Katie_10kb_q10.fastq > ${BASECALL}Katie/Katie_50kb_q10.fastq
seqkit seq -m 10000 ${BASECALL}Katie/Katie_corrected.fasta > ${BASECALL}Katie/Katie_corrected_10kb.fasta
```
**Table 1. Summary of read statistics before and after filtering. Estimated sequence coverage is using an assumed genome size of 1.2Gb.**
|            Read Data Set           |    Basepairs    |    # Reads    | Read N50 (Mbp) | Q20(%) | Q30(%) | Estimated Sequence Coverage |
|:----------------------------------:|:---------------:|:-------------:|:--------------:|:------:|:------:|:---------------------------:|
|              Raw Reads             | 108,533,499,946 |  13,707,986   |      15,545    |  92.2  |  84.5  |             90.4            |
| Simplex Reads minLen 10kb, minQ 10 |  73,976,843,207 |   3,643,945   |      21,626    |  92.3  |  84.7  |             61.6            |
| Simplex Reads minLen 50kb, minQ 10 |   6,169,169,954 |      95,910   |      61,659    |  90.2  |  81.3  |              5.1            |
|            HERRO corrected         |  64,358,311,945 |   3,463,946   |      20,497    |  N/A   |  N/A   |             53.6            |
|      HERRO corrected minLen 10kb   |  62,742,877,123 |   3,192,127   |      20,714    |  N/A   |  N/A   |             52.3            |

## Initial Genome Assembly
Initial genome assemblies using reads trimmed to a minimum Q-score of 20 and a minimium length of either 1, 5 or 10kb were performed using [HiFiasm v0.19.7-r598](https://github.com/chhylp123/hifiasm). Given the likely high inbreeding of the tara iti population, we trialled assembly with and without duplicate purging `-l0`. Given the potentially higher depth provided by the simplex reads, we also trialled assembly using only the simplex reads, with and without duplicate purging.  
```
# Duplicate purging of assembly using HERRO corrected reads
hifiasm -o ${HIFIASM_OUT}Katie_full_HERRO/Katie_50kbUL_HERROfull \
        -t 32 \
        --ul ${BASECALL}Katie/Katie_50kb_q10.fastq \
        ${BASECALL}Katie/Katie_corrected_10kb.fasta

# No duplicate purging of assembly using HERRO corrected reads
hifiasm -o ${HIFIASM_OUT}Katie_full_HERRO_inbred/Katie_50kbUL_HERROfull \
        -t 32 \
        -l0 \
        --ul ${BASECALL}Katie/Katie_50kb_q10.fastq \
        ${BASECALL}Katie/Katie_corrected_10kb.fasta

# Duplicate puring of assembly using simplex reads
hifiasm -o ${HIFIASM_OUT}Katie_simplex/Katie_10kb_q10 \
        -t 32 \
        --ont ${BASECALL}Katie/Katie_10kb_q10.fastq 

# No duplicate purging of assembly using simplex reads
hifiasm -o ${HIFIASM_OUT}Katie_simplex_inbred/Katie_10kb_q10 \
        -t 32 \
        -l0 \
        --ont ${BASECALL}Katie/Katie_10kb_q10.fastq 
```
**Table 1. Initial assembly statistics:**
|                Assembly             |  Haplotype | # Scaffolds | Assembly Size | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) | % single copy BUSCO | % duplicated BUSCO | % missing BUSCO |
|:-----------------------------------:|:----------:|:-----------:|:-------------:|:---------:|:---:|:---------------------------:|:-------------------:|:------------------:|:---------------:|
|   HERRO without duplicate purging   |   primary  |     505     | 1,492,847,865 |   44.9    |  8  |             181             |        96.44%       |         1.18%      |      2.04%      |
|   HERRO without duplicate purging   |  alternate |    6,642    |  782,247,863  |  0.200    | 13  |             6.3             |        40.89%       |         7.10%      |     50.85%      |
|     HERRO with duplicate purging    |      1     |     377     | 1,345,034,214 |   44.9    |  9  |            116.5            |        93.63%       |         0.44%      |      5.47%      |
|     HERRO with duplicate purging    |      2     |     308     | 1,336,940,698 |   47.4    |  7  |            181.0            |        95.39%       |         0.65%      |      3.57%      |
|     HERRO with duplicate purging    |   primary  |     305     | 1,435,061,924 |   47.5    |  8  |            181.0            |        96.88%       |         0.77%      |      2.00%      |
|  Simplex without duplicate purging  |   primary  |     109     | 1,505,468,107 |  49.4     |  7  |            229.5            |        96.59%       |         1.06%      |      1.99%      |
|  Simplex without duplicate purging  |  alternate |    4,958    |  707,091,261  |  0.22     | 796 |             7.1             |        39.69%       |         4.63%      |     54.45%      |
|    Simplex with duplicate purging   |      1     |     106     | 1,403,185,722 |  49.4     |  7  |            229.6            |        94.59%       |         0.40%      |      4.65%      |
|    Simplex with duplicate purging   |      2     |     195     | 1,168,347,444 |  30.2     | 10  |            159.8            |        83.25%       |         0.48%      |     15.59%      |
|    Simplex with duplicate purging   |   primary  |     78      | 1,476,427,735 |  55.7     |  6  |            229.5            |        96.89%       |         0.74%      |      2.00%      |

From this table, we can see the most contiguous and complete genome assemblies are the primary assemblies using HERRO and simplex reads with duplicate purging. We will progress with these two assemblies for scaffolding and polishing.  

## Polishing and Scaffolding
A combination of [Racon v1.5.0](https://github.com/isovic/racon) and [longstitch v1.0.4](https://github.com/bcgsc/LongStitch) were used to polish and scaffold the genome. The below script aligns reads to the draft assembly using [Minimap2 v2.24](https://github.com/lh3/minimap2) for polishing with RACON and scaffolding with Longstitch. Polishing and scaffolding was repeated for two rounds.  
```
mkdir ${dir}${indiv}/${asm}_polish
printf "\nRUNNING RACON USING Q20 5kb ASSEMBLY FOR ${indiv} AT "
date
minimap2 -t 32 -ax map-ont ${dir}${indiv}/dorado_${asm}_simple/${indiv}_${asm}_flye_simple.fa ${reads}${indiv}/${indiv}_${asm}.fastq > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple.sam

printf "FINISHED RUNNING MINIMAP STARTING RACON FOR ${indiv} AT "
date

racon -t 32 ${reads}${indiv}/${indiv}_${asm}.fastq ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple.sam ${dir}${indiv}/dorado_${asm}_simple/${indiv}_${asm}_flye_simple.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple_racon1.fa

printf "FINISHED RACON, STARTING LONGSTITCH AT "
date
        
mkdir ${dir}${indiv}/${asm}_longstitch1
cp ${reads}${indiv}/${indiv}_${asm}.fastq ${dir}${indiv}/${asm}_longstitch1/
cp ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_flye_simple_racon1.fa ${dir}${indiv}/${asm}_longstitch1/

longstitch run -C ${dir}${indiv}/${asm}_longstitch1 draft=${indiv}_${asm}_flye_simple_racon1 reads=${indiv}_${asm} t=32 G=1.2g --debug

cat ${dir}${indiv}/${asm}_longstitch1/${indiv}_${asm}*.ntLink.scaffolds.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa

printf "FINISHED ROUND 1 LONGSTITCH, STARTING ROUND 2 MINIMAP AT "
date

minimap2 -t 64 -ax map-ont ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa ${reads}${indiv}/${indiv}_${asm}.fastq > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.sam

printf "FINISHED ROUND 2 MINIMAP, STARTING ROUND 2 RACON AT "
date
        
racon -t 64 ${reads}${indiv}/${indiv}_${asm}.fastq \
        ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.sam \
        ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1_racon2.fa

printf "FINISHED ROUND 2 RACON, STARTING ROUND 2 LONGSTITCH AT "
date

mkdir ${dir}${indiv}/${asm}_longstitch2
mv ${dir}${indiv}/${asm}_longstitch1/${indiv}_${asm}.fastq ${dir}${indiv}/${asm}_longstitch2/
cp ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch1_racon2.fa ${dir}${indiv}/${asm}_longstitch2/

longstitch run -C ${dir}${indiv}/${asm}_longstitch2 draft=${indiv}_${asm}_longstitch1_racon2 reads=${indiv}_${asm} t=32 G=1.2g --debug

cat ${dir}${indiv}/${asm}_longstitch2/${indiv}_${asm}*.ntLink.scaffolds.fa > ${dir}${indiv}/${asm}_polish/${indiv}_${asm}_longstitch2.fa

printf "\nFINISHED RUNNING 2 ROUNDS OF POLISHING AND SCAFFOLDING OF Q20 5kb ASSEMBLY FOR ${indiv} AT "
date
```

Finally, assembly quality was assessed using [BUSCO](https://busco.ezlab.org/) v5.4.7 to assess how complete it may be, and SeqKit to assess contiguity (Table 2).  

**Table 3. Scaffolded and polished assembly summary statistics:**

|                Assembly             |  Haplotype | # Scaffolds | Assembly Size | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) | % single copy BUSCO | % duplicated BUSCO | % missing BUSCO |
|:-----------------------------------:|:----------:|:-----------:|:-------------:|:---------:|:---:|:---------------------------:|:-------------------:|:------------------:|:---------------:|
|     HERRO with duplicate purging    |   primary  |     305     | 1,435,061,924 |   47.5    |  8  |            181.0            |        96.88%       |         0.77%      |      2.00%      |
|    Simplex with duplicate purging   |   primary  |     78      | 1,476,427,735 |  55.7     |  6  |            229.5            |        96.89%       |         0.74%      |      2.00%      |

[D-GENIES](https://dgenies.toulouse.inra.fr/) (figure below) to visualise structural differences with a HQ assembly for the common tern ([*Sterna hirundo*](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009819605.1/)).  

In the end, the genome assembly generated using a minimum read length of 5kb was used for population read alignment and analyses as it represented a good balance of coverage (Table 1) and contiguity (Table 2).

<figure>
        <div style="text-align: center;">
        <img src="https://github.com/janawold1/2024_MolEcol_ConsGen_Special_Issue/blob/main/Figures/Katie_q20_5kb_longstitch2_to_CommonTern.png"
             alt="Tara iti genome aligned against the common tern genome for comparing synteny and contiguity"
             width="600" height="600">
        </div>
        <figcaption>Tara iti draft assembly mapped to the VGP assembly for Common tern. In this dotplot, the tara iti assembly is represented along the y-axis and the common tern assembly is along the x-axis.</figcaption>
</figure>

Unsurprisingly, the *de novo* tara iti assembly was not chromosomally resolved. Because common terns and fairy terns are relatively related species and demonstrate high synteny, we used the common tern as a reference to scaffold the tara iti assembly with [RagTag](https://github.com/malonge/RagTag) v2.1.0 to maximise our ability to call structural variants.  
```
ragtag.py scaffold -o reference/SP01_ragtag/ reference/common_tern.fasta reference/SP01.fasta
```
**Table 3. Final assembly summary statistics:**
| Read Inputs | % Complete BUSCO | % Missing BUSCO | # Scaffolds | N50 (Mbp) | L50 | Largest Scaffold Size (Mbp) | # N's per 100 kbp |
|:-----------:|:----------------:|:---------------:|:-----------:|:---------:|:---:|:---------------------------:|:-----------------:|
|   Q20, 5kb  |       97.7       |       1.9       |     137     |    84.9   |  5  |            219.3            |       18.7        |

# Using SeqKit Properly
```
seqkit grep -p  "^[^]/|/(\D )(\w+) " -f autosome_names.txt ${REF} > ${OUT_REF}
```