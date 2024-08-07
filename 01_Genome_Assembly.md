# Genome Assembly
## Kmer Analyses
XXXX
### YAK
XX

### Merqury
XX

## HiFiasm
Anahera maternal = JEM; paternal = Waihopai. Te Ariki maternal = Gertrude; paternal = Ariki
```
hifiasm \
    -o ${out}hifiasm/${indiv}/${indiv}_trio_50kb.asm -t 32 \
    -1 mat.yak \
    -2 pat.yak \
    --ul ${reads}${indiv}/${indiv}_50kb_q10.fastq \
    ${reads}${indiv}/${indiv}_corrected_10kb_to_30kb.fasta
```

## Verkko
XXXX
```
verkko -d ${out}verkko/${indiv} \
  --hifi ${reads}${indiv}/${indiv}_corrected_10kb_to_30kb.fasta \
  --nano  ${reads}${indiv}/${indiv}_10kb_q10.fastq \
  --hap-kmers pat_compress.k30.meryl \
              mat_compress.k30.meryl trio \
```

## Assembly QC
XXXX