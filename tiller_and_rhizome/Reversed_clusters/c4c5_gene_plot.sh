awk 'NR==1{print "geneid\tgenename"}NR!=1{print $1"\t"$3}' experiment_marker.txt>geneid
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.2/bin/Rscript \
diff_gene_plot.R \
-i ../data.RDS \
-g geneid \
-f ../diff_gene_SCT/info \
-m SCT



