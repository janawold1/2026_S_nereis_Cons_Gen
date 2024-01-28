### Creating Pangenome Graph for Kakapo Chr7 in Nesi server

Here we use Singularity and follow the procedure in https://github.com/pangenome/pggb#singularity to install pggb instance. 

Run the script https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_50k95.sh for segment size 50,000 and PCT 95. 
Please refer to the log file in this folder. 

This script executes pggb command inside the singularity container image. 

```
#export container to a variable for convenience
container=/nesi/project/ga03793/software/pggb/pggb_0.5.3.simg
data=/nesi/nobackup/ga03793/kakapo_pggb/8kakapo_chr7.fa

#Bind filesystem to container image 
export SINGULARITY_BIND="/nesi/project/ga03793/,/nesi/nobackup/ga03793/"

singularity exec ${container} pggb -i $data -s 50000 -p 98 -n 8 -k 79 -t 20 -S -m -o 8kakapoChr7_pggb0.5.3_50K95_D
```

In `pggb` `-i` is for specifying the sequence file. `-s` specifies the segment length for mapping and `-p` specifies percent identity for mapping/alignment. `-n` is for number of haplotypes. `-k` for minimum matching length. `-t` says number of threads to be used for the execution. `-S` will generate the stats. `-m` will generate MultiQC report of graphs' statistics and visualizations. `-o` specifies the output directory name. 

Chromosome 7 length details
| Name | Length |
|------|-------:|
|A_Jane | 70,131,234    |
|Bill_chr7_as_ref    |   72,941,507    | 
|Blades_chr7_as_ref   |   70,612,238   | 
|Gulliver_chr7_as_ref |   70,588,295    | 
|Huhu_chr7_as_ref     |   71,005,242    | 
|Margaret-Maree_chr7_as_ref  |    77,967,333   | 
|Mati-ma_chr7_as_ref   |  70,547,917    |  
|Sue_chr7_as_ref | 76,248,589    |  

| Alignment Size | % Identity | # Nodes | CPU Threads | Runtime | Mem Usage | No. of SNPs |
|---------------:|-----------:|--------:|--------:|--------:|----------:|------------:|
|      50000     |     98     |    1    |     20    |     01:45:16      |       6.44 GB      | 159,750|
|      50000     |     95     |    1    |     20   |      01:41:50    |       6.95 GB      | 160,620|

The output graphs are as below. 

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og lay draw](https://user-images.githubusercontent.com/8539123/225691565-af307bff-00dc-4ca1-af46-b7d697fb2bc2.png)


![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/225691627-49d67dd7-0133-4dea-a791-7d092d359608.png)


![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/225691680-cbb16c62-b60e-4fd0-af59-e49da4aaef4e.png)

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/225691717-70160b26-68d6-4d8d-81b9-34eaad0c57ea.png)

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/225691757-027eb627-a5d9-465d-8bca-03c396ad9ba9.png)

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/225691792-099d34ce-c79b-446e-879a-b459d00de6c2.png)

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/225691859-b1b68169-390a-4abe-b3f9-53c3404c8491.png)

![8kakapo_chr7 fa d9c9aba c2fac19 c09ef4a smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/225691912-02efcbd4-9650-4f85-bb9f-736b24671191.png)


Run the script https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_50k98.sh for segment size 50,000 and PCT 98. 
Please refer to the log file in this folder. 

The output graphs are as below. 

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og lay draw](https://user-images.githubusercontent.com/8539123/225694594-e9c0335e-a9b8-444a-a3ba-b124bd4aec3a.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/225694640-abde18b8-90c2-42e7-9d0b-5058ebc14e02.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/225694680-40ecdf39-daf9-4a52-9374-225bb14d941c.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/225694728-00ad363e-d8dc-478e-9cbd-721bff504ec1.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/225694779-a5976570-039a-42c9-93de-950dbd7cd002.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/225694842-2012c9ec-bd10-4746-a45c-5d419c6d5680.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/225694894-eafe8f2d-d29e-4eb8-a76c-03d932d70760.png)

![8kakapo_chr7 fa 5075e9f c2fac19 8b5f848 smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/225695167-ee54f498-ebfc-4633-9b44-355eae51fa36.png)

These scripts were executed again using Slurm Workload Manager and observed the resource utilization using the command `seff`

```
seff 33694477 # 50K98
Job ID: 33694477
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 19:32:47
CPU Efficiency: 55.71% of 1-11:05:20 core-walltime
Job Wall-clock time: 01:45:16
Memory Utilized: 6.44 GB
Memory Efficiency: 80.44% of 8.00 GB


seff 33694537 # 50K95
Job ID: 33694537
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 20
CPU Utilized: 19:39:35
CPU Efficiency: 57.92% of 1-09:56:40 core-walltime
Job Wall-clock time: 01:41:50
Memory Utilized: 6.95 GB
Memory Efficiency: 86.89% of 8.00 GB
```

### vg deconstruct and create vcf file using reference as A_Jane

Checking the indexs to get paths

```
cat /nesi/nobackup/ga03793/kakapo_pggb/8kakapo_chr7.fa.fai 
A_Jane  70131234        21      70131234        70131235
Bill_chr7_as_ref        72941507        70131274        60      61
Blades_chr7_as_ref      70612238        144288493       60      61
Gulliver_chr7_as_ref    70588295        216077624       60      61
Huhu_chr7_as_ref        71005242        287842409       60      61
Margaret-Maree_chr7_as_ref      77967333        360031100       60      61
Mati-ma_chr7_as_ref     70547917        439297910       60      61
Sue_chr7_as_ref 76248589        511021643       60      61
```

And also the files generated by `pggb`.

```
ls -ltrh 8kakapoChr7_pggb0.5.3_50K95_D
total 586M
-rw-rw-r-- 1 ismnu81p ismnu81p 1.5K Mar 16 18:02 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.03-16-2023_18:02:20.params.yml
-rw-rw-r-- 1 ismnu81p ismnu81p  19M Mar 16 18:24 8kakapo_chr7.fa.d9c9aba.alignments.wfmash.paf
-rw-rw-r-- 1 ismnu81p ismnu81p 155K Mar 16 22:24 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.fix.affixes.tsv.gz
-rw-rw-r-- 1 ismnu81p ismnu81p 259M Mar 16 22:31 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og
-rw-rw-r-- 1 ismnu81p ismnu81p 174M Mar 16 22:31 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.gfa
-rw-rw-r-- 1 ismnu81p ismnu81p  791 Mar 16 22:32 8kakapo_chr7.fa.d9c9aba.c2fac19.seqwish.og.stats.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p  785 Mar 16 22:32 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.stats.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p 4.0K Mar 16 22:32 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 6.3K Mar 16 22:35 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_pos_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.1K Mar 16 22:39 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_depth_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.6K Mar 16 22:42 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_inv_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 2.5K Mar 16 22:45 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_O_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 7.1K Mar 16 22:48 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.viz_uncalled_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p  89M Mar 16 23:10 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.lay.tsv
-rw-rw-r-- 1 ismnu81p ismnu81p  35M Mar 16 23:10 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.lay
-rw-rw-r-- 1 ismnu81p ismnu81p 298K Mar 16 23:10 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.lay.draw.png
-rw-rw-r-- 1 ismnu81p ismnu81p 714K Mar 16 23:10 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.og.lay.draw_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.1K Mar 16 23:10 multiqc_config.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p 2.1M Mar 16 23:10 multiqc_report.html
drwxrwxr-x 2 ismnu81p ismnu81p 4.0K Mar 16 23:10 multiqc_data
-rw-rw-r-- 1 ismnu81p ismnu81p 7.6M Mar 16 23:10 8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.03-16-2023_18:02:20.log
```

```
ls -ltrh 8kakapoChr7_pggb0.5.3_50K98_D
total 584M
-rw-rw-r-- 1 ismnu81p ismnu81p 1.5K Mar 15 21:35 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.03-15-2023_21:35:24.params.yml
-rw------- 1 ismnu81p ismnu81p  47K Mar 15 21:36 wfmash-6jaGd3
-rw-rw-r-- 1 ismnu81p ismnu81p 7.1K Mar 15 21:39 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.03-15-2023_21:35:24.log
-rw-rw-r-- 1 ismnu81p ismnu81p 1.5K Mar 16 09:57 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.03-16-2023_09:57:10.params.yml
-rw-rw-r-- 1 ismnu81p ismnu81p  19M Mar 16 10:16 8kakapo_chr7.fa.5075e9f.alignments.wfmash.paf
-rw-rw-r-- 1 ismnu81p ismnu81p 154K Mar 16 14:01 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.fix.affixes.tsv.gz
-rw-rw-r-- 1 ismnu81p ismnu81p 258M Mar 16 14:08 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og
-rw-rw-r-- 1 ismnu81p ismnu81p 174M Mar 16 14:08 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.gfa
-rw-rw-r-- 1 ismnu81p ismnu81p  787 Mar 16 14:08 8kakapo_chr7.fa.5075e9f.c2fac19.seqwish.og.stats.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p  784 Mar 16 14:09 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.stats.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p 4.2K Mar 16 14:09 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 6.1K Mar 16 14:12 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_pos_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.3K Mar 16 14:15 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_depth_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.8K Mar 16 14:19 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_inv_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 2.7K Mar 16 14:21 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_O_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 7.2K Mar 16 14:24 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.viz_uncalled_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p  88M Mar 16 14:46 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.lay.tsv
-rw-rw-r-- 1 ismnu81p ismnu81p  35M Mar 16 14:46 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.lay
-rw-rw-r-- 1 ismnu81p ismnu81p 307K Mar 16 14:46 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.lay.draw.png
-rw-rw-r-- 1 ismnu81p ismnu81p 754K Mar 16 14:47 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.og.lay.draw_multiqc.png
-rw-rw-r-- 1 ismnu81p ismnu81p 3.1K Mar 16 14:47 multiqc_config.yaml
-rw-rw-r-- 1 ismnu81p ismnu81p 2.1M Mar 16 14:47 multiqc_report.html
drwxrwxr-x 2 ismnu81p ismnu81p 4.0K Mar 16 14:47 multiqc_data
-rw-rw-r-- 1 ismnu81p ismnu81p 7.2M Mar 16 14:47 8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.03-16-2023_09:57:10.log
```

Executing the command for 50K98

```
vg deconstruct -p A_Jane 8kakapoChr7_pggb0.5.3_50K98_D/8kakapo_chr7.fa.5075e9f.c2fac19.8b5f848.smooth.final.gfa -a -e -K > 8kakapoChr7_A_Jane_pggb0.5.3_50K98_D.vcf

head -20 8kakapoChr7_A_Jane_pggb0.5.3_50K98_D.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=A_Jane,length=70131234>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Bill_chr7_as_ref        Blades_chr7_as_ref      Gulliver_chr7_as_ref   Huhu_chr7_as_ref Margaret-Maree_chr7_as_ref      Mati-ma_chr7_as_ref     Sue_chr7_as_ref
A_Jane  5       >1207>1210      G       C       60      .       AC=1;AF=0.142857;AN=7;AT=>1207>1209>1210,>1207>1208>1210;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  10      >1210>1213      C       A       60      .       AC=1;AF=0.142857;AN=7;AT=>1210>1212>1213,>1210>1211>1213;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  12      >1213>1216      T       C       60      .       AC=1;AF=0.142857;AN=7;AT=>1213>1215>1216,>1213>1214>1216;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  20      >1216>1219      G       C       60      .       AC=1;AF=0.142857;AN=7;AT=>1216>1218>1219,>1216>1217>1219;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  24      >1219>1222      G       A       60      .       AC=1;AF=0.142857;AN=7;AT=>1219>1220>1222,>1219>1221>1222;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  114     >1222>1225      AG      CA      60      .       AC=1;AF=0.142857;AN=7;AT=>1222>1224>1225,>1222>1223>1225;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  138     >1225>1227      C       CA      60      .       AC=2;AF=0.285714;AN=7;AT=>1225>1227,>1225>1226>1227;NS=7;LV=0   GT      0      00       1       0       1       0
A_Jane  160     >1227>1229      CA      C       60      .       AC=1;AF=0.142857;AN=7;AT=>1227>1228>1229,>1227>1229;NS=7;LV=0   GT      1      00       0       0       0       0
A_Jane  163     >1229>1232      A       C       60      .       AC=1;AF=0.142857;AN=7;AT=>1229>1230>1232,>1229>1231>1232;NS=7;LV=0      GT     10       0       0       0       0       0
```
Check stats using `bcftools`

```
bcftools stats 8kakapoChr7_A_Jane_pggb0.5.3_50K98_D.vcf
# This file was produced by bcftools stats (1.16+htslib-1.16) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  8kakapoChr7_A_Jane_pggb0.5.3_50K98_D.vcf
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       8kakapoChr7_A_Jane_pggb0.5.3_50K98_D.vcf
[W::vcf_parse_info] INFO 'CONFLICT' is not defined in the header, assuming Type=String
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      7
SN      0       number of records:      389018
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 159750
SN      0       number of MNPs: 12420
SN      0       number of indels:       215286
SN      0       number of others:       11612
SN      0       number of multiallelic sites:   23176
SN      0       number of multiallelic SNP sites:       680
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       102455  51418   1.99    102152  51103   2.00
# SiS, Singleton stats:
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent [9]repeat-inconsistent   [10]not applicable
SiS     0       1       99809   65565   34244   215716  0       0       215716
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent  [9]repeat-inconsistent  [10]not applicable
AF      0       0.000000        99809   65565   34244   215717  0       0       215717
AF      0       0.150000        1       1       0       0       0       0       0
AF      0       0.160000        1       1       0       0       0       0       0
AF      0       0.180000        6       5       1       4       0       0       4
AF      0       0.190000        18      12      6       18      0       0       18
AF      0       0.220000        13      9       4       11      0       0       11
AF      0       0.240000        70      49      21      40      0       0       40
AF      0       0.270000        2       1       1       1       0       0       1
AF      0       0.280000        19729   13379   6350    11608   0       0       11608
AF      0       0.290000        13      9       4       7       0       0       7
AF      0       0.330000        1590    1111    479     793     0       0       793
AF      0       0.360000        6       4       2       1       0       0       1
AF      0       0.370000        34      23      11      5       0       0       5
AF      0       0.390000        81      54      27      60      0       0       60
AF      0       0.410000        1       0       1       0       0       0       0
AF      0       0.420000        13441   9224    4217    2518    0       0       2518
AF      0       0.440000        12      9       3       1       0       0       1
AF      0       0.450000        5       4       1       2       0       0       2
AF      0       0.490000        1078    738     340     264     0       0       264
AF      0       0.530000        0       0       0       58      0       0       58
AF      0       0.540000        2       2       0       0       0       0       0
AF      0       0.550000        11      6       5       1       0       0       1
AF      0       0.560000        8515    5893    2622    1259    0       0       1259
AF      0       0.570000        0       0       0       12      0       0       12
AF      0       0.590000        70      43      27      25      0       0       25
AF      0       0.610000        22      16      6       1       0       0       1
AF      0       0.630000        1       1       0       1       0       0       1
AF      0       0.660000        621     440     181     112     0       0       112
AF      0       0.690000        10      8       2       1       0       0       1
AF      0       0.700000        5097    3460    1637    771     0       0       771
AF      0       0.740000        7       7       0       6       0       0       6
AF      0       0.770000        4       3       1       0       0       0       0
AF      0       0.790000        36      23      13      23      0       0       23
AF      0       0.810000        1       1       0       0       0       0       0
AF      0       0.820000        215     153     62      37      0       0       37
AF      0       0.840000        2500    1704    796     420     0       0       420
AF      0       0.860000        2       1       1       0       0       0       0
AF      0       0.880000        3       2       1       0       0       0       0
AF      0       0.890000        2       1       1       0       0       0       0
AF      0       0.900000        1       0       1       0       0       0       0
AF      0       0.990000        843     493     350     1361    0       0       1361
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       60.0    153255  102152  51103   215286
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -60     82      0       .
IDD     0       -59     1       0       .
IDD     0       -56     1       0       .
IDD     0       -54     2       0       .
IDD     0       -53     2       0       .
IDD     0       -52     2       0       .
IDD     0       -50     5       0       .
IDD     0       -49     1       0       .
IDD     0       -48     3       0       .
IDD     0       -46     1       0       .
IDD     0       -45     6       0       .
IDD     0       -44     3       0       .
IDD     0       -43     4       0       .
IDD     0       -42     3       0       .
IDD     0       -40     7       0       .
IDD     0       -39     2       0       .
IDD     0       -38     1       0       .
IDD     0       -37     1       0       .
IDD     0       -36     4       0       .
IDD     0       -35     1       0       .
IDD     0       -34     7       0       .
IDD     0       -33     5       0       .
IDD     0       -32     1       0       .
IDD     0       -31     5       0       .
IDD     0       -30     9       0       .
IDD     0       -29     4       0       .
IDD     0       -28     4       0       .
IDD     0       -27     5       0       .
IDD     0       -26     10      0       .
IDD     0       -25     9       0       .
IDD     0       -24     18      0       .
IDD     0       -23     7       0       .
IDD     0       -22     8       0       .
IDD     0       -21     13      0       .
IDD     0       -20     17      0       .
IDD     0       -19     12      0       .
IDD     0       -18     21      0       .
IDD     0       -17     10      0       .
IDD     0       -16     38      0       .
IDD     0       -15     42      0       .
IDD     0       -14     35      0       .
IDD     0       -13     33      0       .
IDD     0       -12     75      0       .
IDD     0       -11     74      0       .
IDD     0       -10     92      0       .
IDD     0       -9      125     0       .
IDD     0       -8      217     0       .
IDD     0       -7      279     0       .
IDD     0       -6      525     0       .
IDD     0       -5      948     0       .
IDD     0       -4      1965    0       .
IDD     0       -3      4420    0       .
IDD     0       -2      15763   0       .
IDD     0       -1      93465   0       .
IDD     0       1       85065   0       .
IDD     0       2       18393   0       .
IDD     0       3       5694    0       .
IDD     0       4       2466    0       .
IDD     0       5       1182    0       .
IDD     0       6       623     0       .
IDD     0       7       383     0       .
IDD     0       8       272     0       .
IDD     0       9       194     0       .
IDD     0       10      174     0       .
IDD     0       11      114     0       .
IDD     0       12      95      0       .
IDD     0       13      76      0       .
IDD     0       14      63      0       .
IDD     0       15      63      0       .
IDD     0       16      68      0       .
IDD     0       17      44      0       .
IDD     0       18      65      0       .
IDD     0       19      39      0       .
IDD     0       20      40      0       .
IDD     0       21      38      0       .
IDD     0       22      26      0       .
IDD     0       23      33      0       .
IDD     0       24      23      0       .
IDD     0       25      31      0       .
IDD     0       26      22      0       .
IDD     0       27      21      0       .
IDD     0       28      17      0       .
IDD     0       29      23      0       .
IDD     0       30      23      0       .
IDD     0       31      14      0       .
IDD     0       32      19      0       .
IDD     0       33      13      0       .
IDD     0       34      11      0       .
IDD     0       35      11      0       .
IDD     0       36      8       0       .
IDD     0       37      19      0       .
IDD     0       38      14      0       .
IDD     0       39      16      0       .
IDD     0       40      15      0       .
IDD     0       41      12      0       .
IDD     0       42      9       0       .
IDD     0       43      8       0       .
IDD     0       44      9       0       .
IDD     0       45      7       0       .
IDD     0       46      11      0       .
IDD     0       47      2       0       .
IDD     0       48      11      0       .
IDD     0       49      6       0       .
IDD     0       50      7       0       .
IDD     0       51      15      0       .
IDD     0       52      11      0       .
IDD     0       53      12      0       .
IDD     0       54      11      0       .
IDD     0       55      12      0       .
IDD     0       56      13      0       .
IDD     0       57      12      0       .
IDD     0       58      5       0       .
IDD     0       59      9       0       .
IDD     0       60      1053    0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     6695
ST      0       A>G     23757
ST      0       A>T     6149
ST      0       C>A     7175
ST      0       C>G     5525
ST      0       C>T     27370
ST      0       G>A     26953
ST      0       G>C     5537
ST      0       G>T     7326
ST      0       T>A     6166
ST      0       T>C     24375
ST      0       T>G     6845
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
```

Executing the command for 50K95

```
vg deconstruct -p A_Jane 8kakapoChr7_pggb0.5.3_50K95_D/8kakapo_chr7.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.gfa -a -e -K> 8kakapoChr7_A_Jane_pggb0.5.3_50K95_D.vcf

head -20 8kakapoChr7_A_Jane_pggb0.5.3_50K95_D.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=A_Jane,length=70131234>
#CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  Bill_chr7_as_ref        Blades_chr7_as_ref      Gulliver_chr7_as_ref   Huhu_chr7_as_ref Margaret-Maree_chr7_as_ref      Mati-ma_chr7_as_ref     Sue_chr7_as_ref
A_Jane  5       >2292>2295      G       C       60      .       AC=1;AF=0.142857;AN=7;AT=>2292>2294>2295,>2292>2293>2295;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  10      >2295>2298      C       A       60      .       AC=1;AF=0.142857;AN=7;AT=>2295>2297>2298,>2295>2296>2298;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  12      >2298>2301      T       C       60      .       AC=1;AF=0.142857;AN=7;AT=>2298>2300>2301,>2298>2299>2301;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  20      >2301>2304      G       C       60      .       AC=1;AF=0.142857;AN=7;AT=>2301>2303>2304,>2301>2302>2304;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  24      >2304>2307      G       A       60      .       AC=1;AF=0.142857;AN=7;AT=>2304>2306>2307,>2304>2305>2307;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  114     >2307>2310      AG      CA      60      .       AC=1;AF=0.142857;AN=7;AT=>2307>2309>2310,>2307>2308>2310;NS=7;LV=0      GT     10       0       0       0       0       0
A_Jane  138     >2310>2312      C       CA      60      .       AC=2;AF=0.285714;AN=7;AT=>2310>2312,>2310>2311>2312;NS=7;LV=0   GT      0      00       1       0       1       0
A_Jane  160     >2312>2314      CA      C       60      .       AC=1;AF=0.142857;AN=7;AT=>2312>2313>2314,>2312>2314;NS=7;LV=0   GT      1      00       0       0       0       0
A_Jane  163     >2314>2317      A       C       60      .       AC=1;AF=0.142857;AN=7;AT=>2314>2316>2317,>2314>2315>2317;NS=7;LV=0      GT     10       0       0       0       0       0
```

Check stats using `bcftools`

```
bcftools stats 8kakapoChr7_A_Jane_pggb0.5.3_50K95_D.vcf
# This file was produced by bcftools stats (1.16+htslib-1.16) and can be plotted using plot-vcfstats.
# The command line was: bcftools stats  8kakapoChr7_A_Jane_pggb0.5.3_50K95_D.vcf
#
# Definition of sets:
# ID    [2]id   [3]tab-separated file names
ID      0       8kakapoChr7_A_Jane_pggb0.5.3_50K95_D.vcf
[W::vcf_parse_info] INFO 'CONFLICT' is not defined in the header, assuming Type=String
# SN, Summary numbers:
#   number of records   .. number of data rows in the VCF
#   number of no-ALTs   .. reference-only sites, ALT is either "." or identical to REF
#   number of SNPs      .. number of rows with a SNP
#   number of MNPs      .. number of rows with a MNP, such as CC>TT
#   number of indels    .. number of rows with an indel
#   number of others    .. number of rows with other type, for example a symbolic allele or
#                          a complex substitution, such as ACT>TCGA
#   number of multiallelic sites     .. number of rows with multiple alternate alleles
#   number of multiallelic SNP sites .. number of rows with multiple alternate alleles, all SNPs
# 
#   Note that rows containing multiple types will be counted multiple times, in each
#   counter. For example, a row with a SNP and an indel increments both the SNP and
#   the indel counter.
# 
# SN    [2]id   [3]key  [4]value
SN      0       number of samples:      7
SN      0       number of records:      391294
SN      0       number of no-ALTs:      0
SN      0       number of SNPs: 160620
SN      0       number of MNPs: 12559
SN      0       number of indels:       216609
SN      0       number of others:       11626
SN      0       number of multiallelic sites:   23322
SN      0       number of multiallelic SNP sites:       676
# TSTV, transitions/transversions:
# TSTV  [2]id   [3]ts   [4]tv   [5]ts/tv        [6]ts (1st ALT) [7]tv (1st ALT) [8]ts/tv (1st ALT)
TSTV    0       102965  51720   1.99    102667  51403   2.00
# SiS, Singleton stats:
# SiS   [2]id   [3]allele count [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent [9]repeat-inconsistent   [10]not applicable
SiS     0       1       100577  66036   34541   217210  0       0       217210
# AF, Stats by non-reference allele frequency:
# AF    [2]id   [3]allele frequency     [4]number of SNPs       [5]number of transitions        [6]number of transversions      [7]number of indels     [8]repeat-consistent  [9]repeat-inconsistent  [10]not applicable
AF      0       0.000000        100577  66036   34541   217211  0       0       217211
AF      0       0.160000        0       0       0       3       0       0       3
AF      0       0.180000        10      8       2       7       0       0       7
AF      0       0.190000        18      12      6       24      0       0       24
AF      0       0.220000        15      12      3       2       0       0       2
AF      0       0.240000        58      37      21      46      0       0       46
AF      0       0.270000        2       1       1       4       0       0       4
AF      0       0.280000        19905   13483   6422    11667   0       0       11667
AF      0       0.290000        10      5       5       6       0       0       6
AF      0       0.330000        1449    1014    435     712     0       0       712
AF      0       0.360000        6       5       1       1       0       0       1
AF      0       0.370000        31      21      10      5       0       0       5
AF      0       0.390000        89      61      28      65      0       0       65
AF      0       0.410000        1       0       1       2       0       0       2
AF      0       0.420000        13565   9310    4255    2547    0       0       2547
AF      0       0.440000        6       5       1       1       0       0       1
AF      0       0.450000        3       2       1       1       0       0       1
AF      0       0.490000        995     687     308     233     0       0       233
AF      0       0.530000        0       0       0       30      0       0       30
AF      0       0.540000        1       1       0       0       0       0       0
AF      0       0.550000        5       3       2       1       0       0       1
AF      0       0.560000        8537    5910    2627    1289    0       0       1289
AF      0       0.590000        70      47      23      30      0       0       30
AF      0       0.610000        23      16      7       2       0       0       2
AF      0       0.630000        1       1       0       0       0       0       0
AF      0       0.660000        576     417     159     104     0       0       104
AF      0       0.690000        9       6       3       1       0       0       1
AF      0       0.700000        5142    3490    1652    782     0       0       782
AF      0       0.740000        6       6       0       9       0       0       9
AF      0       0.770000        4       4       0       0       0       0       0
AF      0       0.790000        22      12      10      10      0       0       10
AF      0       0.810000        1       1       0       0       0       0       0
AF      0       0.820000        219     152     67      37      0       0       37
AF      0       0.840000        2490    1708    782     422     0       0       422
AF      0       0.860000        2       1       1       0       0       0       0
AF      0       0.880000        2       1       1       0       0       0       0
AF      0       0.890000        3       2       1       0       0       0       0
AF      0       0.900000        2       1       1       0       0       0       0
AF      0       0.990000        830     487     343     1380    0       0       1380
# QUAL, Stats by quality
# QUAL  [2]id   [3]Quality      [4]number of SNPs       [5]number of transitions (1st ALT)      [6]number of transversions (1st ALT)    [7]number of indels
QUAL    0       60.0    154070  102667  51403   216609
# IDD, InDel distribution:
# IDD   [2]id   [3]length (deletions negative)  [4]number of sites      [5]number of genotypes  [6]mean VAF
IDD     0       -60     78      0       .
IDD     0       -59     1       0       .
IDD     0       -56     1       0       .
IDD     0       -54     3       0       .
IDD     0       -53     2       0       .
IDD     0       -52     2       0       .
IDD     0       -50     4       0       .
IDD     0       -49     1       0       .
IDD     0       -48     3       0       .
IDD     0       -47     1       0       .
IDD     0       -46     1       0       .
IDD     0       -45     3       0       .
IDD     0       -44     1       0       .
IDD     0       -43     4       0       .
IDD     0       -42     3       0       .
IDD     0       -40     11      0       .
IDD     0       -39     2       0       .
IDD     0       -38     2       0       .
IDD     0       -37     2       0       .
IDD     0       -36     4       0       .
IDD     0       -35     3       0       .
IDD     0       -34     6       0       .
IDD     0       -33     5       0       .
IDD     0       -32     1       0       .
IDD     0       -31     4       0       .
IDD     0       -30     10      0       .
IDD     0       -29     3       0       .
IDD     0       -28     6       0       .
IDD     0       -27     4       0       .
IDD     0       -26     11      0       .
IDD     0       -25     8       0       .
IDD     0       -24     20      0       .
IDD     0       -23     6       0       .
IDD     0       -22     11      0       .
IDD     0       -21     13      0       .
IDD     0       -20     13      0       .
IDD     0       -19     14      0       .
IDD     0       -18     20      0       .
IDD     0       -17     11      0       .
IDD     0       -16     35      0       .
IDD     0       -15     38      0       .
IDD     0       -14     32      0       .
IDD     0       -13     37      0       .
IDD     0       -12     72      0       .
IDD     0       -11     73      0       .
IDD     0       -10     99      0       .
IDD     0       -9      134     0       .
IDD     0       -8      216     0       .
IDD     0       -7      284     0       .
IDD     0       -6      520     0       .
IDD     0       -5      942     0       .
IDD     0       -4      1917    0       .
IDD     0       -3      4438    0       .
IDD     0       -2      15865   0       .
IDD     0       -1      93906   0       .
IDD     0       1       85733   0       .
IDD     0       2       18551   0       .
IDD     0       3       5768    0       .
IDD     0       4       2495    0       .
IDD     0       5       1188    0       .
IDD     0       6       641     0       .
IDD     0       7       379     0       .
IDD     0       8       277     0       .
IDD     0       9       197     0       .
IDD     0       10      171     0       .
IDD     0       11      118     0       .
IDD     0       12      96      0       .
IDD     0       13      75      0       .
IDD     0       14      63      0       .
IDD     0       15      68      0       .
IDD     0       16      68      0       .
IDD     0       17      41      0       .
IDD     0       18      65      0       .
IDD     0       19      41      0       .
IDD     0       20      41      0       .
IDD     0       21      37      0       .
IDD     0       22      29      0       .
IDD     0       23      34      0       .
IDD     0       24      23      0       .
IDD     0       25      34      0       .
IDD     0       26      21      0       .
IDD     0       27      20      0       .
IDD     0       28      18      0       .
IDD     0       29      24      0       .
IDD     0       30      22      0       .
IDD     0       31      14      0       .
IDD     0       32      22      0       .
IDD     0       33      10      0       .
IDD     0       34      13      0       .
IDD     0       35      10      0       .
IDD     0       36      11      0       .
IDD     0       37      20      0       .
IDD     0       38      13      0       .
IDD     0       39      18      0       .
IDD     0       40      13      0       .
IDD     0       41      13      0       .
IDD     0       42      8       0       .
IDD     0       43      8       0       .
IDD     0       44      10      0       .
IDD     0       45      8       0       .
IDD     0       46      12      0       .
IDD     0       47      1       0       .
IDD     0       48      11      0       .
IDD     0       49      9       0       .
IDD     0       50      8       0       .
IDD     0       51      15      0       .
IDD     0       52      12      0       .
IDD     0       53      11      0       .
IDD     0       54      14      0       .
IDD     0       55      15      0       .
IDD     0       56      12      0       .
IDD     0       57      11      0       .
IDD     0       58      7       0       .
IDD     0       59      9       0       .
IDD     0       60      1052    0       .
# ST, Substitution types:
# ST    [2]id   [3]type [4]count
ST      0       A>C     6731
ST      0       A>G     23935
ST      0       A>T     6180
ST      0       C>A     7167
ST      0       C>G     5569
ST      0       C>T     27472
ST      0       G>A     27104
ST      0       G>C     5559
ST      0       G>T     7411
ST      0       T>A     6234
ST      0       T>C     24454
ST      0       T>G     6869
# DP, Depth distribution
# DP    [2]id   [3]bin  [4]number of genotypes  [5]fraction of genotypes (%)    [6]number of sites      [7]fraction of sites (%)
```
