### Creating Pangenome Graph for Kakapo Chr Z 

Chromosome Z length details 

| Name | Length |
|------|-------:|
|Jane:Chromosome | 101,869,369 |
|Gulliver:Chromosome_Z |  110,514,568 |
|Huhu:Chromosome_Z   |    91,693,272  |
|Margaret-Maree:Chromosome_Z  |   112,049,481  |
|Mati-ma:Chromosome_Z  |  85,189,130   |
|Sue:Chromosome_Z   |     86,075,096   |

| Alignment Size | % Identity | # Nodes | CPU Threads | Runtime | Mem Usage | No. of SNPs |
|---------------:|-----------:|--------:|--------:|--------:|----------:|------------:|
|      50,000     |     98     |    1    |     24    |      02:31:42      |       8.41 GB      | 388,444 |
|      50,000     |     95     |    1    |     24   |      02:51:00    |       8.04 GB      | 444,352 |

Please refer 2 scripts [https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_chZ_50k95.sh](https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_chZ_50k95.sh) and [https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_chZ_50k98.sh](https://github.com/nuzla/kakapo_pangenome/blob/main/scripts/pggb_0.5.3_8kakapo_chZ_50k98.sh)

Executing the SLURM job. 

```
sbatch ./pggb_0.5.3_8kakapo_chZ_50k98.sh
Submitted batch job 33689447
```

```
sbatch ./pggb_0.5.3_8kakapo_chZ_50k95.sh
Submitted batch job 33689553
```

Parameters used

```
#SBATCH --account       ga03793
#SBATCH --job-name      kakapo_chrZ_50K95
#SBATCH --cpus-per-task 24
#SBATCH --mem           16G
#SBATCH --time          24:00:00
```

We can see the job states using `squeue --me`

```
JOBID         USER     ACCOUNT   NAME        CPUS MIN_MEM PARTITI START_TIME     TIME_LEFT STATE    NODELIST(REASON)    
33686533      ismnu81p ga03793   spawner-jupy   4     32G bigmem  2023-03-27T1     6:00:56 RUNNING  wbl002              
33689447      ismnu81p ga03793   kakapo_chrZ_  20     16G large   2023-03-27T1    23:53:25 RUNNING  wbn074              
33689553      ismnu81p ga03793   kakapo_chrZ_  20     16G large   2023-03-27T1    23:57:24 RUNNING  wbn198 
```

Also more details using `seff` command. 

```
seff 33689447
Job ID: 33689447
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: RUNNING
Nodes: 1
Cores per node: 24
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 02:42:40 core-walltime
Job Wall-clock time: 00:08:08
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 16.00 GB (16.00 GB/node)
WARNING: Efficiency statistics may be misleading for RUNNING jobs.

seff 33689553
Job ID: 33689553
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: RUNNING
Nodes: 1
Cores per node: 24
CPU Utilized: 00:00:00
CPU Efficiency: 0.00% of 01:26:20 core-walltime
Job Wall-clock time: 00:04:19
Memory Utilized: 0.00 MB (estimated maximum)
Memory Efficiency: 0.00% of 16.00 GB (16.00 GB/node)
WARNING: Efficiency statistics may be misleading for RUNNING jobs.
```

`seff` output after completing the jobs

```
seff 33689447 #for 50K98
Job ID: 33689447
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 24
CPU Utilized: 1-10:35:44
CPU Efficiency: 57.01% of 2-12:40:48 core-walltime
Job Wall-clock time: 02:31:42
Memory Utilized: 8.41 GB
Memory Efficiency: 52.59% of 16.00 GB
```

```
seff 33689553 #for 50K95
Job ID: 33689553
Cluster: mahuika
User/Group: ismnu81p/ismnu81p
State: COMPLETED (exit code 0)
Nodes: 1
Cores per node: 24
CPU Utilized: 1-13:40:05
CPU Efficiency: 55.07% of 2-20:24:00 core-walltime
Job Wall-clock time: 02:51:00
Memory Utilized: 8.04 GB
Memory Efficiency: 50.28% of 16.00 GB
```

The output graphs for `-s 50000` and `-p 98`

![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og lay draw](https://user-images.githubusercontent.com/8539123/229793909-fbb7dfc5-48be-4b67-a192-ba913f9b27b9.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/229793942-a4a09695-22dd-4550-a54f-b6b5c4bd4501.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/229793970-22debc8e-0be8-496f-8aa9-2a3fdea144f4.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/229794002-51aae94a-43d5-4558-8948-b32447afaad3.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/229794038-7f4a9cab-ab7f-4259-a4b9-b3329d55ff26.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/229794058-9b137122-0939-453e-bd6b-7aab17e6fc3c.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/229794093-14f4657c-5ad0-4e39-9e59-4b406f9bda04.png)
![8kakapo_chrZ fa d4be715 c2fac19 f731fd2 smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/229794115-1c5e8bcb-31e4-4341-bc22-67e9512aa3fc.png)


The output graphs for `-s 50000` and `-p 95`

![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og lay draw](https://user-images.githubusercontent.com/8539123/229794211-30a45bb4-7834-49bb-8aa9-aff6e9ce0d98.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og lay draw_multiqc](https://user-images.githubusercontent.com/8539123/229794241-d618b452-1438-4dd8-a9cd-1b9fb9c65e8d.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_depth_multiqc](https://user-images.githubusercontent.com/8539123/229794279-99c43615-9ead-4fd2-9b97-510fb6fa5d99.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_inv_multiqc](https://user-images.githubusercontent.com/8539123/229794305-8f5ee0ce-43be-438b-ad6d-8f27390e2191.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_multiqc](https://user-images.githubusercontent.com/8539123/229794343-5f3fbc7e-e60e-4902-9d2b-081e9e585017.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_O_multiqc](https://user-images.githubusercontent.com/8539123/229794374-d527d544-609c-4219-baeb-e1d3aadde31f.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_pos_multiqc](https://user-images.githubusercontent.com/8539123/229794420-47973c1e-711b-4a36-8666-1212b91dc06f.png)
![8kakapo_chrZ fa 5a120b8 c2fac19 c47d9e7 smooth final og viz_uncalled_multiqc](https://user-images.githubusercontent.com/8539123/229794452-e4d6e34b-f269-44f8-888d-7e5f8376ef81.png)


Create vcf file using `vg deconstruct` and check stats. But the command got failed. This should be a bug of the vg tool. 

```
vg deconstruct -p Jane:Chromosome  8kakapoChrZ_pggb0.5.3_50K98_D/8kakapo_chrZ.fa.5075e9f.c2fac19.8b5f848.smooth.final.gfa -a -e -K > 8kakapoChrZ_pggb0.5.3_50K98_D.vcf
Bus Error
```

```
vg deconstruct -p Jane:Chromosome  8kakapoChrZ_pggb0.5.3_50K95_D/8kakapo_chrZ.fa.d9c9aba.c2fac19.c09ef4a.smooth.final.gfa -a -e -K > 8kakapoChrZ_pggb0.5.3_50K95_D.vcf
Bus Error
```

After few researches the command was successfull with option `export OMP_NUM_THREADS=1`.

```
head -20 8kakapoChrZ_pggb0.5.3_50K98_D.vcf
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Gulliver:Chromosome_Z	Huhu:Chromosome_Z	Margaret-Maree:Chromosome_Z	Mati-ma:Chromosome_Z	Sue:Chromosome_Z
Jane:Chromosome	7549	>3>5	G	GTTATCTCGCGATTGTGGTTATCTTATCT	60	.	AC=1;AF=1;AN=1;AT=>3>5,>3>4>5;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7550	>5>7	T	TAC	60	.	AC=1;AF=1;AN=1;AT=>5>7,>5>6>7;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7551	>7>9	T	TCA	60	.	AC=1;AF=1;AN=1;AT=>7>9,>7>8>9;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7552	>9>11	A	AGCAGTGTCCGAGTTTCTGTCGG	60	.	AC=1;AF=1;AN=1;AT=>9>11,>9>10>11;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7553	>11>13	T	TTC	60	.	AC=1;AF=1;AN=1;AT=>11>13,>11>12>13;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7554	>13>15	C	CTGGAAG	60	.	AC=1;AF=1;AN=1;AT=>13>15,>13>14>15;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7556	>15>17	C	CTGGCG	60	.	AC=1;AF=1;AN=1;AT=>15>17,>15>16>17;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7557	>17>19	G	GTTATCT	60	.	AC=1;AF=1;AN=1;AT=>17>19,>17>18>19;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7558	>19>21	C	CGT	60	.	AC=1;AF=1;AN=1;AT=>19>21,>19>20>21;NS=1;LV=0	GT	1	.	.	.	.
```

```
head -20 8kakapoChrZ_pggb0.5.3_50K95_D.vcf
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
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	Gulliver:Chromosome_Z	Huhu:Chromosome_Z	Margaret-Maree:Chromosome_Z	Mati-ma:Chromosome_Z	Sue:Chromosome_Z
Jane:Chromosome	7549	>3>5	G	GTTATCTCGCGATTGTGGTTATCTTATCT	60	.	AC=1;AF=1;AN=1;AT=>3>5,>3>4>5;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7550	>5>7	T	TAC	60	.	AC=1;AF=1;AN=1;AT=>5>7,>5>6>7;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7551	>7>9	T	TCA	60	.	AC=1;AF=1;AN=1;AT=>7>9,>7>8>9;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7552	>9>11	A	AGCAGTGTCCGAGTTTCTGTCGG	60	.	AC=1;AF=1;AN=1;AT=>9>11,>9>10>11;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7553	>11>13	T	TTC	60	.	AC=1;AF=1;AN=1;AT=>11>13,>11>12>13;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7554	>13>15	C	CTGGAAG	60	.	AC=1;AF=1;AN=1;AT=>13>15,>13>14>15;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7556	>15>17	C	CTGGCG	60	.	AC=1;AF=1;AN=1;AT=>15>17,>15>16>17;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7557	>17>19	G	GTTATCT	60	.	AC=1;AF=1;AN=1;AT=>17>19,>17>18>19;NS=1;LV=0	GT	1	.	.	.	.
Jane:Chromosome	7558	>19>21	C	CGT	60	.	AC=1;AF=1;AN=1;AT=>19>21,>19>20>21;NS=1;LV=0	GT	1	.	.	.	.
```

Counting SNPs using `bcftools stats`.

```
bcftools stats 8kakapoChrZ_pggb0.5.3_50K98_D.vcf | grep "SNPs:"
[W::vcf_parse] INFO 'CONFLICT' is not defined in the header, assuming Type=String
SN	0	number of SNPs:	388444

bcftools stats 8kakapoChrZ_pggb0.5.3_50K95_D.vcf | grep "SNPs:"
[W::vcf_parse] INFO 'CONFLICT' is not defined in the header, assuming Type=String
SN	0	number of SNPs:	444352
```
