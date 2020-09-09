BASE_PATH="/HOME_ann/BII/biipsashare/"
INPUT_FOLDER="/HOME_ann/BII/biipsashare/winston/input/"
OUTPUT_FOLDER="/HOME_ann/BII/biipsashare/winston/output/"
SCRIPTS="/HOME_ann/BII/biipsashare/winston/GISAID-Analysis-Update/"
PYTHON="~/miniconda3/envs/phylo/bin/python"

eval ${PYTHON} ${SCRIPTS}extractor.py ${BASE_PATH}raphael/treeupdate$(date --date="1 days ago" '+%Y%m%d') ${INPUT_FOLDER}

sleep 10

eval ${PYTHON} ${SCRIPTS}graph_plot.py ${INPUT_FOLDER}clade_progression.tsv ${INPUT_FOLDER}geoclade.tsv -co ${OUTPUT_FOLDER}clade.png -go ${OUTPUT_FOLDER}geoclade.png

eval ${PYTHON} ${SCRIPTS}nwk_to_fig.py -i ${INPUT_FOLDER}fnb_all_rapidnj.nwk -o ${OUTPUT_FOLDER}

eval ${PYTHON} ${SCRIPTS}phylo.py ${INPUT_FOLDER}fnb_all_rapidnj.nwk -o ${OUTPUT_FOLDER}phylo.svg >> phylo_gen.log

eval zip ${OUTPUT_FOLDER}tree_image_gen.zip ${OUTPUT_FOLDER}clade.png ${OUTPUT_FOLDER}geoclade.png ${OUTPUT_FOLDER}$(date '+%Y-%m_d')_Full ${OUTPUT_FOLDER}phylo.svg ${OUTPUT_FOLDER}phylo_gen.log