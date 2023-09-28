#!/bin/bash -e
#SBATCH --account ga03793
#SBATCH --job-name Te_Ariki_Graph_self_alignments
#SBATCH --partition=milan
#SBATCH --cpus-per-task=16
#SBATCH --mem=64G
#SBATCH --time 16:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Self alignment to genome assemblies.
###############################################################

ml purge
ml load minimap2/2.24-GCC-11.3.0
ml load SAMtools/1.16.1-GCC-11.3.0

# And now beginning to run the associated shell script.
dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/
graph=/scale_wlg_nobackup/filesets/nobackup/uc03718/graphs/

for indiv in bird
        do
        echo "BEGIN ALIGNING READS FOR MINIMAP ${indiv} AT "
        date
        minimap2 -ax map-ont ${graph}whole-genome_minimap_asRef/${indiv}_q10_longstitch_racon2_asRef.fasta.gz \
        ${dir}variants/reads/${indiv}_dorado_q15_1kb.fq | \
        samtools view -@32 -b > ${graph}05_longstitch_minimap_asRef/${indiv}_longstitch_minimap_asRef_algn.bam
     samtools sort -@32 ${graph}05_longstitch_minimap_asRef/${indiv}_longstitch_minimap_asRef_algn.bam > ${graph}05_longstitch_minimap_asRef/${indiv}_longstitch_minimap_asRef_algn.sorted.bam
        echo "FINSIHED ALIGNING READS FOR MINIMAP ${indiv} AT "
        date
        echo "BEGAN ALIGNING READS FOR MASHMAP ${indiv} AT "
        date
        minimap2 -ax map-ont ${graph}whole-genome_mashmap_asRef/${indiv}_racon2_longstitch_mashmap_asRef.fa \
        ${dir}variants/reads/${indiv}_dorado_q15_1kb.fq | \
        samtools view -@32 -b > ${graph}racon_longstitch_mashmap_asRef/${indiv}_longstitch_mashmap_asRef_algn.bam
        samtools sort -@32 ${graph}racon_longstitch_mashmap_asRef/${indiv}_longstitch_mashmap_asRef_algn.bam > ${graph}racon_longstitch_mashmap_asRef/${indiv}_longstitch_mashmap_asRef_algn.sorted.bam
        echo "ESTIMATING COVERAGE FOR ${indiv} AT "
        date
        samtools coverage ${graph}racon_longstitch_minimap_asRef/${indiv}_longstitch_minimap_asRef_algn.sorted.bam > ${graph}racon_longstitch_minimap_asRef/${indiv}_self_coverage.tsv
        samtools coverage ${graph}racon_longstitch_mashmap_asRef/${indiv}_longstitch_mashmap_asRef_algn.sorted.bam > ${graph}racon_longstitch_mashmap_asRef/${indiv}_self_coverage.tsv
        echo "FINISHED ALL ALIGNMENTS SORTING BAMS AND COVERAGE STATS FOR $indiv AT "
        date
done