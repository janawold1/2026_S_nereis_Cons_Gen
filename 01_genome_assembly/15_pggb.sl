#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name longstitch_asRef
#SBATCH --cpus-per-task=24
#SBATCH --mem=64G
#SBATCH --time 08:30:00
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load Singularity/3.11.3

# All chromosomes fasta need to be bgzipped and indexed with "samtools faidx"

# Export container to a variable for convenience
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/graphs/racon_longstitch_mashmap_asRef/
container=${dir}pggb_0.5.3.simg

# Bind filesystem to container image
export SINGULARITY_BIND="/scale_wlg_nobackup/filesets/nobackup/uc03718/graphs//racon_longstitch_mashmap_asRef/,/scale_wlg_nobackup/filesets/nobackup/uc03718/graphs/racon_longstitch_mashmap_asRef/"

cd ${dir}

pwd

for i in {1..24}
        do
        for segment in 30 50 100
                do
                for percent in 90 95 98
                        do
                        printf "RUNNING PGGB FOR A SEGMENT SIZE OF ${segment} AND A PRECENT IDENTITY OF ${percent} FOR CHROMOSOME ${i} AT "
                        date
                        singularity exec ${container} pggb -i ${dir}chromosome_${i}/chromosome_${i}.fa.gz \
                                -p ${percent} \
                                -s ${segment}000 \
                                -n 3 -k 79 -t 24 -T 24 -S -m -D ./ \
                                -V 'Jane:#' \
                                -o ${dir}chromosome_${i}/kakapo_seg${segment}kb_${percent}perc_k79
                        printf "FINISHED SEGEMENT ${segment} AND ${percent} PERCENT FOR CHROMOSOME ${i} AT "
                        date
                done
        done
done
