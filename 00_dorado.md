# Simplex Basecalling with Dorado v0.7.2
First performed simplex basecalling for MinION and PromethION R10 flow cells. All of Anahera's reads were run on the 5Hz flowcell, while some of Te Ariki's MinION sequencing was run prior to the conversion from 4 to 5Hz. We performed all basecalling under the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model for all 5Hz `pod5` files and `XXX` model for Te Ariki's 4Hz reads.  

For either case, basecalling was performed as below.  
```
dorado download --model $model

dorado basecaller $model pod5_dir/ > ${out}/${indiv}_SUPv5.bam
```

The basecalled reads were then converted to fastq files with SAMtools v1.19 and filtered to a minmum length of 10kb and minimum Q-score of 10 with SeqKit v2.4 to prepare for read correction.  
```
samtools fastq -@32 ${out}/${indiv}_SUPv5.bam > ${out}/${indiv}_SUPv5.fastq
seqkit seq -m 10000 -Q 10 ${out}/${indiv}_SUPv5.fastq > ${out}/${indiv}_10kb_q10.fastq

dorado correct -v ${out}/${indiv}_10kb_q10.fastq > ${out}/${indiv}_corrected.fasta
```
## Read Preparation for Genome assembly
### HERRO Corrected Reads
As recommended by ONT representatives, HERRO corrected reads >30kb were trimmed to be between 10 and 30kb in length with the python script below. 
```
#!/usr/bin/env python3

import argparse
from Bio import SeqIO

def trim_to_length(record, max_length, min_length):
    trimmed_records = []
    while len(record.seq) > max_length:
        # Create a new record for the trimmed part
        trimmed_record = record[:]
        trimmed_record.seq = record.seq[max_length:]
        trimmed_record.id += "_part"
        trimmed_record.description += " (trimmed part)"
        trimmed_records.append(trimmed_record)

        # Trim the current record
        record.seq = record.seq[:max_length]

    if len(record.seq) >= min_length:
        trimmed_records.insert(0, record)
    return trimmed_records

def main():
    parser = argparse.ArgumentParser(description="Trim sequences in a FASTA file to a specified length range.")
    parser.add_argument("input_file", help="Input FASTA file")
    parser.add_argument("output_file", help="Output FASTA file for sequences within the length range")
    parser.add_argument("trimmed_file", help="Output FASTA file for trimmed parts of sequences")
    parser.add_argument("--max_length", type=int, default=30000, help="Maximum sequence length (default: 30000)")
    parser.add_argument("--min_length", type=int, default=10000, help="Minimum sequence length (default: 10000)")
    
    args = parser.parse_args()

    with open(args.output_file, "w") as output_handle, open(args.trimmed_file, "w") as trimmed_handle:
        for record in SeqIO.parse(args.input_file, "fasta"):
            if len(record.seq) < args.min_length:
                continue  # Skip sequences shorter than min_length
            
            trimmed_records = trim_to_length(record, args.max_length, args.min_length)
            
            for trimmed_record in trimmed_records:
                if len(trimmed_record.seq) < args.min_length:
                    SeqIO.write(trimmed_record, trimmed_handle, "fasta")
                else:
                    SeqIO.write(trimmed_record, output_handle, "fasta")

if __name__ == "__main__":
    main()
```
XXX
`./trim_fasta.py input.fasta output_trimmed.fasta trimmed_parts.fasta --max_length 30000 --min_length 10000`

### Ultralong ONT Reads 
Although ultralong reads are often defined as >100kb in length, we do not have suifficient depth. To alleviate this, we used reads >50kb in length and a Q-score >10 as input for HiFiasm. Below is a summary table of the basecalled and trimmed reads for each of our samples.  

|     Sample    |         Read Group       | Total Number of Reads | Total number of Base Pairs | Minimum Length | Maximum Length | Average Length | Read N50 | Q30%  |
|:-------------:|:------------------------:|:---------------------:|:--------------------------:|:--------------:|:--------------:|:--------------:|:--------:|:-----:|
|    Female 1   |     Raw Simplex Reads    |       X,XXX,XXX       |                            |       1        |                |                |          | XX.XX |
|     Male 1    |     Raw Simplex Reads    |       X,XXX,XXX       |                            |       1        |                |                |          | XX.XX |
| Richard Henry |     Raw Simplex Reads    |      13,775,740       |        45,945,262,782      |       1        |   1,012,251    |     3,335.2    |   5,635  | 82.71 |
|    Female 1   |    Simplex Reads >10kb   |       1,904,040       |        68,026,262,057      |     10,000     |     445,955    |    35,727.3    |  47,915  | 84.70 |
|     Male 1    |    Simplex Reads >10kb   |       X,XXX,XXX       |                            |     10,000     |                |                |          | XX.XX |
| Richard Henry |    Simplex Reads >10kb   |        653,301        |         8,890,372,442      |     10,000     |     198,129    |    13,608.4    |  13,182  | 84.48 |
|    Female 1   |    Simplex Reads >50kb   |        431,500        |        32,331,592,724      |     50,000     |     445,955    |    74,928.4    |  74,737  | 84.64 |
|     Male 1    |    Simplex Reads >50kb   |        383,007        |        28,567,762,973      |     50,000     |     391,680    |    74,588.1    |  74,642  | 85.61 |
| Richard Henry |    Simplex Reads >50kb   |          169          |          12,059,962        |     50,000     |     198,129    |    71,360.7    |  64,588  | 58.04 |
|    Female 1   |      Raw HERRO Reads     |       1,950,680       |        64,872,812,188      |       1        |     238,214    |    33,256.5    |  45,143  |   NA  |
|     Male 1    |      Raw HERRO Reads     |        XXX,XXX        |                            |       1        |                |                |          |   NA  |
| Richard Henry |      Raw HERRO Reads     |        XXX,XXX        |                            |       1        |                |                |          |   NA  |
|    Female 1   |     HERRO Reads >10kb    |       1,855,805       |        64,319,617,592      |     10,000     |     238,214    |    34,658.6    |  45,409  |   NA  |
|     Male 1    |     HERRO Reads >10kb    |       X,XXX,XXX       |        XX,XXX,XXX,XXX      |     10,000     |      XX,XXX    |    XX,XXX.X    |  XX,XXX  |   NA  |
| Richard Henry |     HERRO Reads >10kb    |       X,XXX,XXX       |                            |                |                |                |          |   NA  |
|    Female 1   | HERRO Reads >10kb, <30kb |       2,657,765       |        62,539,906,387      |     10,000     |      30,067    |    23,531.0    |  30,000  |   NA  |
|     Male 1    | HERRO Reads >10kb, <30kb |       2,403,372       |        55,888,871,709      |     10,000     |      30,000    |    23,254.4    |  30,000  |   NA  |
| Richard Henry | HERRO Reads >10kb, <30kb |       X,XXX,XXX       |                            |                |                |                |          |   NA  |