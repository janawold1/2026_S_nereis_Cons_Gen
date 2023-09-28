#!/bin/bash -e
#SBATCH --account uc03718
#SBATCH --partition=milan
#SBATCH --job-name Anahera_dorado_busco
#SBATCH --cpus-per-task=32
#SBATCH --mem=64G
#SBATCH --time 010:00:00  # Walltime (HH:MM:SS)
#SBATCH --output %x.%j.out
#SBATCH --error %x.%j.err
#SBATCH --profile=ALL

ml purge
ml load BUSCO/5.4.7-gimkl-2022a

dir=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/
aves=/scale_wlg_nobackup/filesets/nobackup/uc03718/initial_flye/busco_downloads/lineages/aves_odb10/

for indiv in Anahera
        do
        printf "\nRUNNING BUSCO FOR DORADO FLYE ASSEMBLY FOR ${indiv} AT "
        date
        cd ${dir}${indiv}
        busco --in dorado_q10_racon/${indiv}_q10_flye_racon2.fasta \    
                --out dorado_q10_flye_busco/ \
                --mode genome --lineage_dataset $aves --cpu 32
        printf "\nRUNNING BUSCO FOR DORADO LONGSTITCH ASSEMBLY FOR ${indiv} AT "
        date
        busco --in dorado_q10_racon/${indiv}_q10_longstitch_racon2.fasta \
                --out dorado_q10_longstitch_busco/ \
                --mode genome --lineage_dataset ${aves} --cpu 32 --restart
        printf "\nFINISHED RUNNING QUAST FOR ${indiv} AT "
        date
done
