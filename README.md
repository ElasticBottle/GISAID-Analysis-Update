# Introduction

Global initiative on sharing all influenza data (GISAID) is a global science initiative and primary source that provides open-access to genomic data of influenza viruses and the novel coronavirus responsible for COVID-19.

This repository implements the workflow that translates newly scraped data from the GISAID database and summarizes them.

![Phylogenetic tree with clade coloring](images/phylo.png)

## Details

The summary includes:

- Phylogenetic tree categorized by the common clade type occurring in hcov-19 and color coded based on region
- Graph of new cases broken down by clade type
- Graph of the clade progression
- Convertor from nwk file to custom figtree file
- Utility for extracting .tar files

## Setting up

1. Ensure that you have python 3.6 installed:
    - `pip install pipenv`
2. After which, clone the tool to your desired location
   - `git clone https://github.com/ElasticBottle/GISAID-Analysis-Update.git`
3. Then run  `pipenv install` in the `directory` of the tool.

## Using the tool

### Generating Phylogenetic Tree

In the directory of the tool:

Run in shell one of the following commands for `FastTree`, `RapinNj`, and `iqTreeFile` respectively.

```python
python phylo_updated.py path_to\fasttree.tree --output path_to_output\fast.svg --root-on "S_" --clade-color-density-coverage .\clade_del.json --display --rapid-fasttree-color-marker
```

```python
python phylo_updated.py path_to\fnb_all_rapidnj.nwk --output path_to_output\rapid.svg --root-on "S_" --clade-color-density-coverage .\clade_del.json --display --rapid-fasttree-color-marker --quoted
```

```python
python phylo_updated.py path_to\iqtree.treefile --output path_to_output\iq.svg --root-on "S_" --clade-color-density-coverage .\clade_del.json --display --iqtree-color-marker
```

For more options, run `python phylo.py -h` in shell.

Alternatively, you might want to use the web based [phyloTreeMaker](https://mendel3.bii.a-star.edu.sg/METHODS/corona/current/phyloTreeMaker/build/) built off this script.

### Generating Clade Progression and Geo Clade Graphs

In the directory of the tool:

Run in shell

```python
python graph_plot.py path_to\clade_progression.tsv path_to\geoclade.tsv -co path_to\clade.png -go path_to\geo.png
```

Run `python graph_plot.py -h` for more detials on the various options.

## Extras

### NWK to FigTree file convertor

[FigTree](http://tree.bio.ed.ac.uk/software/figtree/) saves it's file in a custom NEXUS file format

```python
python nwk_to_fig.py -i INPUT_FILE_LOC -o OUTPUT_FILE_LOC
```

`python nwk_to_fig.py -h` for the full details.

Due to confidentiality, I am unable to release the test files used to test the convertor.

### utils: .tar file extractor

`python extractor.py input_tar_file.tar output_folder` in the dir of the tool.

Will create a new folder if the folder path does not exist.
