awk -F'\t' '{if($7=="4"||$7=="5"){print $8"\t"$7"\t"$18}}' ../markergene/markers_0.15.info.txt|awk 'BEGIN{print "gene\tcluster\tgenename"}{if($3==""){$3=$1} ;split($3,name,",");$3=name[1];print $1"\t"$2"\t"$3}'>experiment_marker.txt
/hwfssz1/ST_EARTH/Reference/ST_AGRIC/USER/zhongliyuan/miniconda3/envs/R4.1/bin/Rscript \
known_marker.R \
 -i ../data.RDS\
 -v csv.file\
 -g experiment_marker.txt
 

