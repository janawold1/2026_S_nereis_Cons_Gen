# Simplex Basecalling with Dorado v0.7.2
First performed simplex basecalling for MinION and PromethION R10 flow cells. All of Anahera's reads were run on the 5Hz flowcell, while some of Te Ariki's MinION sequencing was run prior to the conversion from 4 to 5Hz. We performed all basecalling under the `dna_r10.4.1_e8.2_400bps_sup@v5.0.0` model for all 5Hz `pod5` files and `XXX` model for Te Ariki's 4Hz reads.  

For either case, basecalling was performed as below.  
```
dorado download --model $model

dorado basecaller $model pod5_dir/ > ${out}/${indiv}_SUPv5.bam
```

The basecalled reads were then converted to fastq files with SAMtools v1.19 and filtered to a minmum length of 10kb with SeqKit v2.4 to prepare for read correction.  
```
samtools fastq -@32 ${out}/${indiv}_SUPv5.bam > ${out}/${indiv}_SUPv5.fastq
seqkit seq -m 10000 ${out}/${indiv}_SUPv5.fastq > ${out}/${indiv}_10kb.fastq

dorado correct -v ${out}/${indiv}_10kb.fastq > ${out}/${indiv}_corrected.fasta
```
## Read Preparation for Genome assembly
### HERRO Corrected Reads
XXXX
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
XXXX

### Core-Simplex Reads
XXXX

Summary table of Reads