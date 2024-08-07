
### Creating Pangenome Graph for Kakapo Chr 7 + Chr Z 

Chromosome 7+Z length details

| Name | length | 
|:-----|-------:|
|A_Jane                       |   70,131,234   |
|Bill_chr7_as_ref             |   72,941,507   |
|Blades_chr7_as_ref           |   70,612,238   |
|Gulliver_chr7_as_ref         |   70,588,295   |
|Huhu_chr7_as_ref             |   71,005,242   |
|Margaret-Maree_chr7_as_ref   |   77,967,333   |
|Mati-ma_chr7_as_ref          |   70,547,917   |
|Sue_chr7_as_ref              |   76,248,589   |
|Jane:Chromosome              |   101,869,369  |
|Gulliver:Chromosome_Z        |   110,514,568  |
|Huhu:Chromosome_Z            |   91,693,272   |
|Margaret-Maree:Chromosome_Z  |   112,049,481  |
|Mati-ma:Chromosome_Z         |   85,189,130   |
|Sue:Chromosome_Z             |   86,075,096   |


| Alignment Size | % Identity | # Nodes | CPU Threads | Runtime | Mem Usage | No. of SNPs |
|---------------:|-----------:|--------:|--------:|--------:|----------:|------------:|
|      50,000     |     98     |    1    |     24    |      04:37:53      |       11.66 GB      | 388,944|
|      50,000     |     95     |    1    |     24   |      04:34:10    |       10.92 GB      | 445,519|

Job execution results
```
seff 33716924  #50K98
Job ID: 33716924
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 24
CPU Utilized: 2-08:19:36
CPU Efficiency: 50.67% of 4-15:09:12 core-walltime
Job Wall-clock time: 04:37:53
Memory Utilized: 11.66 GB
Memory Efficiency: 72.85% of 16.00 GB


seff 33716950 ##50K95
Job ID: 33716950
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 24
CPU Utilized: 2-11:56:39
CPU Efficiency: 54.66% of 4-13:40:00 core-walltime
Job Wall-clock time: 04:34:10
Memory Utilized: 10.92 GB
Memory Efficiency: 68.23% of 16.00 GB
```

The output graphs for `-s 50000` and `-p 98`

![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og lay draw](https://user-images.githubusercontent.com/8539123/229852863-d678aaad-4c68-420d-b70e-5121d0a157da.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/229852898-0b5a2d1d-f5ca-4b06-aa22-35ccbc242b73.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/229852926-794695ca-c16f-4fa0-a5e9-21f0164e388a.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/229852990-aeccd4b0-d0c3-4ce5-a6af-4bb7b99fbaa2.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/229853030-57766a59-218a-44ee-9633-b32929fc5074.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/229853057-fa8cb013-4fe6-4b88-81be-6116079eef17.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/229853084-1f1eda19-22b9-44ee-afe8-baed37da9922.png)
![8kakapo_chr7z fa 9789424 c2fac19 a393cd4 smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/229853128-ccde8c0a-7694-4bee-9a41-1ded236ae0bf.png)


The output graphs for `-s 50000` and `-p 95`
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og lay draw](https://user-images.githubusercontent.com/8539123/229853419-00d7fd46-7c65-4961-a0f3-70743e83ea8e.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/229853448-461cca94-2482-4c68-96ec-1a344267d8a5.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/229853485-b3d6e6de-f10f-4199-99ad-43716d62d58d.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/229853523-fc314830-8d9a-46e1-8e48-c33ab4168f10.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/229853551-3d4899ad-b301-4327-bfcc-5725070eadf9.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/229853590-fbda2304-a5dd-49b9-8dc8-4cd478428c95.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/229853614-a7eb9e83-4250-4319-9ff3-0117eb1dd69f.png)
![8kakapo_chr7z fa 7f5b031 c2fac19 9365171 smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/229853654-a569a3b2-ecdd-48c9-b258-ae174dce529e.png)


Create vcf file using `vg deconstruct` and the reference as Jane:Chromosome. 

```
export OMP_NUM_THREADS=1

vg deconstruct -p Jane:Chromosome  8kakapoChr7z_pggb0.5.3_50K98_D_2/8kakapo_chr7z.fa.9789424.c2fac19.a393cd4.smooth.final.gfa -a -e -K > 8kakapoChr7z_pggb0.5.3_50K98_D_2.vcf

vg deconstruct -p Jane:Chromosome  8kakapoChr7z_pggb0.5.3_50K95_D_2/8kakapo_chr7z.fa.7f5b031.c2fac19.9365171.smooth.final.gfa -a -e -K > 8kakapoChr7z_pggb0.5.3_50K95_D_2.vcf
```

```
head -20 8kakapoChr7z_pggb0.5.3_50K98_D_2.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=Jane:Chromosome,length=101869369>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A_Jane	Bill_chr7_as_ref	Blades_chr7_as_ref	Gulliver:Chromosome_Z	Gulliver_chr7_as_ref	Huhu:Chromosome_Z	Huhu_chr7_as_ref	Margaret-Maree:Chromosome_Z	Margaret-Maree_chr7_as_ref	Mati-ma:Chromosome_Z	Mati-ma_chr7_as_ref	Sue:Chromosome_Z	Sue_chr7_as_ref
Jane:Chromosome	7604	>902959>902961	C	CGGTTCCTGGAAGTCTGGCGGTTATCTCGTGATTGTGGTTATCTTATCTTACTCAAGCAGTGTCCGAGTTTCTGTC	60	.	AC=1;AF=1;AN=1;AT=<902961<902959,<902961<902960<902959;NS=1;LV=0	GT	.	.	.	1	.	.	.	.	.	.	.	.	.
Jane:Chromosome	7799	>902957>902959	G	GG	60	.	AC=1;AF=1;AN=1;AT=<902959<902957,<902959<902958<902957;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	9256	>902955>902957	GG	G	60	.	AC=1;AF=1;AN=1;AT=<902957<902956<902955,<902957<902955;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	9261	>902952>902955	A	G	60	.	AC=1;AF=1;AN=1;AT=<902955<902954<902952,<902955<902953<902952;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	11729	>902950>902952	GG	G	60	.	AC=1;AF=1;AN=1;AT=<902952<902951<902950,<902952<902950;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	11750	>902948>902950	GGG	G	60	.	AC=1;AF=1;AN=1;AT=<902950<902949<902948,<902950<902948;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	12689	>902944>902948	GCTTGGCCTATGCCTACAGAGCCCTGCTCAACACTATCCACCACCTTCAAAGGG	GGCTTGGCCTATGCCTACAGAGCCCTGCTCAACACTATCCACCACCTTCAAAGGGA	60	.	AC=1;AF=1;AN=1;AT=<902948<902947<902944,<902948<902946<902945<902944;NS=1;LV=0	GT	.	.	.	1	.	.	.	.	.	.	.	.	.
Jane:Chromosome	17785	>902942>902944	G	GG	60	.	AC=1;AF=1;AN=1;AT=<902944<902942,<902944<902943<902942;NS=1;LV=0	GT	.	.	.	1	.  .
Jane:Chromosome	18289	>902940>902942	CTCCCC	C	60	.	AC=1;AF=1;AN=1;AT=<902942<902941<902940,<902942<902940;NS=1;LV=0	GT	.	.	.	1	.  .
```

```
head -20 8kakapoChr7z_pggb0.5.3_50K95_D_2.vcf
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Total number of alternate alleles in called genotypes">
##INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1]">
##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">
##INFO=<ID=LV,Number=1,Type=Integer,Description="Level in the snarl tree (0=top level)">
##INFO=<ID=PS,Number=1,Type=String,Description="ID of variant corresponding to parent snarl">
##INFO=<ID=AT,Number=R,Type=String,Description="Allele Traversal as path in graph">
##contig=<ID=Jane:Chromosome,length=101869369>
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A_Jane	Bill_chr7_as_ref	Blades_chr7_as_ref	Gulliver:Chromosome_Z	Gulliver_chr7_as_ref	Huhu:Chromosome_Z	Huhu_chr7_as_ref	Margaret-Maree:Chromosome_Z	Margaret-Maree_chr7_as_ref	Mati-ma:Chromosome_Z	Mati-ma_chr7_as_ref	Sue:Chromosome_Z	Sue_chr7_as_ref
Jane:Chromosome	7604	>2321876>2321878	C	CGGTTCCTGGAAGTCTGGCGGTTATCTCGTGATTGTGGTTATCTTATCTTACTCAAGCAGTGTCCGAGTTTCTGTC	60	.	AC=1;AF=1;AN=1;AT=<2321878<2321876,<2321878<2321877<2321876;NS=1;LV=0	GT	.	.	.	1	.	.	.	.	.	.	.	.	.
Jane:Chromosome	7799	>2321874>2321876	G	GG	60	.	AC=1;AF=1;AN=1;AT=<2321876<2321874,<2321876<2321875<2321874;NS=1;LV=0	GT	.	.	.	1  .
Jane:Chromosome	9256	>2321872>2321874	GG	G	60	.	AC=1;AF=1;AN=1;AT=<2321874<2321873<2321872,<2321874<2321872;NS=1;LV=0	GT	.	.	.	1  .
Jane:Chromosome	9261	>2321869>2321872	A	G	60	.	AC=1;AF=1;AN=1;AT=<2321872<2321870<2321869,<2321872<2321871<2321869;NS=1;LV=0	GT	.	.	.  .
Jane:Chromosome	11729	>2321867>2321869	GG	G	60	.	AC=1;AF=1;AN=1;AT=<2321869<2321868<2321867,<2321869<2321867;NS=1;LV=0	GT	.	.	.	1  .
Jane:Chromosome	11750	>2321865>2321867	GGG	G	60	.	AC=1;AF=1;AN=1;AT=<2321867<2321866<2321865,<2321867<2321865;NS=1;LV=0	GT	.	.	.	1  .
Jane:Chromosome	12689	>2321862>2321865	GCTTGGCCTATGCCTACAGAGCCCTGCTCAACACTATCCACCACCTTCAAAGGG	GGCTTGGCCTATGCCTACAGAGCCCTGCTCAACACTATCCACCACCTTCAAAGGGA	60	.	AC=1;AF=1;AN=1;AT=<2321865<2321864<2321862,<2321865<2321863<2321862;NS=1;LV=0	GT	.	.	.	1	.	.	.	.	.	.	.	.	.
Jane:Chromosome	17785	>2321860>2321862	G	GG	60	.	AC=1;AF=1;AN=1;AT=<2321862<2321860,<2321862<2321861<2321860;NS=1;LV=0	GT	.	.	.	1  .
Jane:Chromosome	18289	>2321858>2321860	CTCCCC	C	60	.	AC=1;AF=1;AN=1;AT=<2321860<2321859<2321858,<2321860<2321858;NS=1;LV=0	GT	.	.	.	1  .
```

Counting SNPs using `bcftools stats`.

```
bcftools stats 8kakapoChr7z_pggb0.5.3_50K98_D_2.vcf | grep "SNPs:"
[W::vcf_parse] INFO 'CONFLICT' is not defined in the header, assuming Type=String
SN	0	number of SNPs:	388944

bcftools stats 8kakapoChr7z_pggb0.5.3_50K95_D_2.vcf | grep "SNPs:"
[W::vcf_parse] INFO 'CONFLICT' is not defined in the header, assuming Type=String
SN	0	number of SNPs:	445519
```
