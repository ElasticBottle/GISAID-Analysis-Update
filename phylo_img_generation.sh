#! /bin/bash

# Folder paths
BASE_PATH="/HOME_ann/BII/biipsashare/"
INPUT_FOLDER="/HOME_ann/BII/biipsashare/winston/input/"
OUTPUT_FOLDER="/HOME_ann/BII/biipsashare/winston/output/"
SCRIPTS="/HOME_ann/BII/biipsashare/winston/GISAID-Analysis-Update/"
PYTHON="~/miniconda3/envs/phylo/bin/python"

#Script Paths
EXTRACTOR=${SCRIPTS}extractor.py
GRAPH_PLOTS=${SCRIPTS}graph_plot.py
NWK_TO_FIG=${SCRIPTS}nwk_to_fig.py
PHYLO=${SCRIPTS}phylo.py

# Input Paths
CLADE_PROGRESSION_INPUT=${INPUT_FOLDER}clade_progression.tsv
GEO_CLADES_INPUT=${INPUT_FOLDER}geoclade.tsv
NWK_FILE=${INPUT_FOLDER}fnb_all_rapidnj.nwk

# Output Paths
PHYLO_LOGS=${OUTPUT_FOLDER}phylo_gen.log
CLADE_OUTPUT=${OUTPUT_FOLDER}clade.png
GEO_CLADES_OUTPUT=${OUTPUT_FOLDER}geoclade.png
PHYLO_TREE_SVG_OUTPUT=${OUTPUT_FOLDER}phylo.svg
ZIP_FILE_OUTPUT${OUTPUT_FOLDER}tree_image_gen.zip

# Extracting the files
eval ${PYTHON} ${EXTRACTOR} ${BASE_PATH}raphael/treeupdate$(date --date="1 days ago" '+%Y%m%d').tar.xz ${INPUT_FOLDER}

sleep 10

# Plotting the graphs
eval ${PYTHON} ${GRAPH_PLOTS} ${CLADE_PROGRESSION_INPUT} ${GEO_CLADES_INPUT} -co ${CLADE_OUTPUT} -go ${GEO_CLADES_OUTPUT}

# Converting and cleaning the nwk file
eval ${PYTHON} ${NWK_TO_FIG} -i ${NWK_FILE} -o ${OUTPUT_FOLDER}

# Remove old log file, generating new Phylo tree and logs
eval rm ${PHYLO_LOGS}

eval ${PYTHON} ${PHYLO} ${NWK_FILE} -o ${PHYLO_TREE_SVG_OUTPUT} >> ${PHYLO_LOGS}

# Remove old zip files and generates a new one
eval  rm ${ZIP_FILE_OUTPUT}

eval zip ${ZIP_FILE_OUTPUT} ${CLADE_OUTPUT} ${GEO_CLADES_OUTPUT} ${OUTPUT_FOLDER}$(date '+%Y-%m-%d')_Full ${PHYLO_TREE_SVG_OUTPUT}

# Sends mail with Logs as content and zip file as attachment
eval mail -a ${ZIP_FILE_OUTPUT} -s "Tree auto generated" winstonyeo99@yahoo.com < ${PHYLO_LOGS}