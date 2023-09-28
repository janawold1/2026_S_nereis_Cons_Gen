# Read mapping with VG

Although indexing the `.gfa` directly is encouraged, it is extremely computationally intensive and takes a long time to prune sufficiently before it can be indexed. For instance, attempted to index chromosome 7 `.gfa` with `vg autoindex` ran for 48 hours and never pruned sufficiently to finish. To speed things along, the `.gfa` can be converted to packed to default HashGraph `.vg` format.  
```
vg convert -g chr7*.final.gfa > chr7.vg
```
Then, to reduce graph complexity and index, the graph was modified to reduce kmer offset to <=1024 and unlink nodes with an edge degree > 34.
```
vg mod -M 34 -X 1024 chr7.vg > chr7_M34_X1024.vg
vg index -p -b ./ -x chr7_M34_X1024.xg -g chr7_M34_X1024.gcsa chr7_M34_X1024.vg
```
Next, long-reads were mapped to each graph with `vg map`.
```
for chrom in {1..24}
do
    for graph in kakapo_seg{30,50,100}kb_perc{90,95,98}_k79
    do
        for indiv in bird1 bird2
        do
            cd ${dir}chromosome_${i}/${graph}/
            echo "STARTED MAPPING $indiv READS TO CHROMOSOME $i AND GRAPH ${graph}..."
            vg map --base-name chr${i}_M34_X1024 \
                --threads 32 --alignment-model long \
                --log-time \
                --fastq ${dir}reads/${indiv}_self_align_chr${i}.fq > ${indiv}_chr${i}.gam
            echo "FINISHED MAPPING ${indiv} READS TO CHROMOSOME $i AND GRAPH ${graph}..."
        done
    done
done
```

Average read mapping scores were assessed with:
```
vg view -h -a bird.gam | jq '.mapping_quality' | awk '{sum += $1}; END {print sum/NR}'
```