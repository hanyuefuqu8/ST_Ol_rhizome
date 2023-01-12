#input csv and annotated cluster
#if you don't have annotated cluster, use rds to generate one
#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Seurat_cluster.R -i data.RDS
#/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Giotto.R -i /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/soybean_seed/soybean_root/FP200000239BL_A2/roottip/bin40/BIN40.csv  -c 3 -g 100 -d 10 -l cluster.csv

/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript /hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/Script/Giotto.R -i /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/soybean_seed/soybean_root/FP200000239BL_A2/roottip/bin40/BIN40.csv -c 3 -g 100 -d 10 -l /hwfssz1/ST_EARTH/P20Z10200N0035/USER/zhongliyuan/zhongliyuan/soybean_seed/soybean_root/FP200000239BL_A2/roottip/bin40/BIN40.csv -k 12 
