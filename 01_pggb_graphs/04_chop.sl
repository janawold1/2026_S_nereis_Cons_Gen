#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --job-name dorado_chop
#SBATCH --partition=milan
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time 00:30:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

###############################################################
# Trimming ONT reads for genome assembly.
###############################################################

R10=/scale_wlg_nobackup/filesets/nobackup/uc03718/dorado/R10_basecalled_reads/
lib=lib1

for indiv in A B C
        do
        ml purge; ml load Porechop/0.2.4-gimkl-2020a-Python-3.8.2
        echo BEGAN R10 PORECHOP FOR ${indiv} AT
        date
        porechop -t 32 -i ${R10}${indiv}/${indiv}_${lib}_R10.fq --discard_middle \
                -o ${R10}${indiv}/${indiv}_${lib}_R10.chopped.fq
        echo FINISHED R10 PORECHOP FOR ${indiv} AT
        date
        ml purge; ml load chopper/0.5.0-GCC-11.3.0
        echo BEGAN R10 TRIMMING FOR ${indiv} AT
        date
        cat ${R10}${indiv}/${indiv}_${lib}_R10.chopped.fq | chopper -q 20 > ${R10}${indiv}/${indiv}_${lib}_q20_R10.fq &
        cat ${R10}${indiv}/${indiv}_${lib}_R10.chopped.fq | chopper -q 20 -l 1000 > ${R10}${indiv}/${indiv}_${lib}_q20_1kb_R10.fq &
        cat ${R10}${indiv}/${indiv}_${lib}_R10.chopped.fq | chopper -q 10 -l 1000 > ${R10}${indiv}/${indiv}_${lib}_q10_1kb_R10.fq &
        cat ${R10}${indiv}/${indiv}_${lib}_R10.chopped.fq | chopper -q 15 -l 1000 > ${R10}${indiv}/${indiv}_${lib}_q10_R10.fq
        echo FINISHED R10 TRIMMING FOR ${indiv} AT
        date
done