#!/bin/bash
wd=/var/www/nodejs/iris3-backend/plumber/1612742383433
exp_file=Zeisel_expression.csv
label_file=Zeisel_index_label.csv
gene_module_file=
jobid=1612742383433
motif_min_length=12
motif_max_length=12
perl /var/www/html/iris3/program/prepare_email1.pl $jobid
Rscript /var/www/html/iris3/program/genefilter.R $jobid $wd$exp_file , $label_file , No 0.5 20 5000 2 No
echo gene_module_detection > running_status.txt
/var/www/html/iris3/program/qubic2/qubic -i $wd$jobid\_filtered_expression.txt -q 0.06 -c 1.0 -k 20-o 500 -f 0.7
for file in *blocks
do
grep Conds $file |cut -d ':' -f2 >'$(basename $jobid_blocks.conds.txt)'
done
for file in *blocks
do
grep Genes $file |cut -d ':' -f2 >'$(basename $jobid_blocks.gene.txt)'
done
Rscript /var/www/html/iris3/program/ari_score.R $label_file $jobid , 2 
echo gene_module_assignment > running_status.txt
Rscript /var/www/html/iris3/program/cts_gene_list.R $wd $jobid 1000   
echo motif_finding_and_comparison > running_status.txt
/var/www/html/iris3/program/get_motif.sh $wd $motif_min_length $motif_max_length 1
Rscript /var/www/html/iris3/program/convert_meme.R $wd $motif_min_length
/var/www/html/iris3/program/get_motif.sh $wd $motif_min_length $motif_max_length 0
wait
cd $wd
Rscript /var/www/html/iris3/program/prepare_bbc.R $jobid $motif_min_length
mkdir tomtom
mkdir logo_tmp
mkdir logo
mkdir regulon_id
/var/www/html/iris3/program/get_logo.sh $wd
/var/www/html/iris3/program/get_tomtom.sh $wd
echo active_regulon_determination > running_status.txt
Rscript /var/www/html/iris3/program/merge_tomtom.R $wd $jobid $motif_min_length
echo regulon_inference > running_status.txt
Rscript /var/www/html/iris3/program/sort_regulon.R $wd $jobid
/var/www/html/iris3/program/get_atac_overlap.sh $wd
Rscript /var/www/html/iris3/program/prepare_heatmap.R $wd $jobid 2 
Rscript /var/www/html/iris3/program/get_alternative_regulon.R $jobid
Rscript /var/www/html/iris3/program/generate_rss_scatter.R $jobid
Rscript /var/www/html/iris3/program/process_tomtom_result.R $jobid
mkdir json
/var/www/html/iris3/program/build_clustergrammar.sh $wd $jobid 2 
zip -R $wd$jobid '*.regulon_gene_id.txt' '*.regulon_gene_symbol.txt' '*.regulon_rank.txt' '*_silh.txt' '*umap_embeddings.txt' '*.regulon_activity_score.txt' '*_cell_label.txt' '*.blocks' '*_blocks.conds.txt' '*_blocks.gene.txt' '*_filtered_expression.txt' '*_gene_id_name.txt' '*_marker_genes.txt' 'cell_cluster_unique_diffrenetially_expressed_genes.txt' '*_combine_regulon.txt'
perl /var/www/html/iris3/program/prepare_email.pl $jobid
echo 'finish'> done
chmod -R 777 .
rm $wd$jobid\_filtered_expression.txt
rm $wd$jobid\_filtered_expression.txt.chars
