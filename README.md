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

## Using the Tool

### Setting up

1. Ensure that you have python 3.8 installed:
    - `pip install pipenv`
2. After which, clone the tool to your desired location
   - `git clone https://github.com/ElasticBottle/GISAID-Analysis-Update.git`
3. Then run  `pipenv install` in the `dir` of the tool.

### Generating Phylogenetic Tree

In the directory of the tool:

Run in shell

```python
python phylo.py INPUT_FILE.nwk
```

Optionally output path can be specified with `-o OUTPUT_FOLDER` flag

For more options, run `python phylo.py -h` in shell

### Generating Clade Progression and Geo Clade Graphs

In the directory of the tool:

```python
python graph_plot.py CLADE_PROGRESSION_INPUT_FILE GEOCLADE_INPUT_FILE
```

Optionally `-co CLADE_PROG_OUT_FILE.png` and `-go  GEOCLADE_OUT_FILE.png` can be used to specify the output file path for the respective graphs.

Run `python graph_plot.py -h` for all the options.

### NWK to FigTree file convertor

[FigTree](http://tree.bio.ed.ac.uk/software/figtree/) saves it's file in a custom NEXUS file format

This tool eliminates the need for opening and saving the file to generate the proper formatting for the Phylogenetic tree to display properly.

This tool assumes that your tree **has quoted node names**. Should be a quick fix to cater to non-quoted names if need be.

```python 
python nwk_to_fig.py -i INPUT_FILE_LOC -o OUTPUT_FILE_LOC
```

both input and output file parameters are optional.

`python nwk_to_fig.py -h` for the full details.

Due to confidentiality, I am unable to release the test files used to test the convertor.

### utils: .tar file extractor

`python extractor.py input_tar_file.tar output_folder` in the dir of the tool.
