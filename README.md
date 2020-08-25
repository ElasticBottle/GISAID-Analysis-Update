# Introduction

Global initiative on sharing all influenza data (GISAID) is a global science initiative and primary source that provides open-access to genomic data of influenza viruses and the novel coronavirus responsible for COVID-19.

This repository implements the workflow that translates newly scraped data from the GISAID database and summarizes them.

## Details

The summary includes:

- Phylogenetic tree categorized by the common clade type occurring in hcov-19 and color coded based on region
- Graph of new cases broken down by clade type
- Graph of the clade progression

## Using the Tool

### Setting up

1. Ensure that you have python 3.8 installed:
    - `pip install pipenv`
2. After which, clone the tool.
3. Then run  `pipenv install` in the `dir` of the tool.

### Generating Phylogenetic Tree

In the directory of the tool:

`python phylo.py "INPUT_FILE.nwk"`

Alternatively, run `python phylo.py -h` in the tool for help

### Generating Clade Progression and Geo Clade Representation

Run 'python graph_plot.py' in the dir of the tool for help.

## Todo(elasticBottle): fill this in
