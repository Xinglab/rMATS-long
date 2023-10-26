# rMATS-long shiny

## About

A shiny server for visualizing rMATS-long results

## Table of Contents

* [Dependencies](#dependencies)
* [Usage](#usage)
  + [Data](#data)
  + [Shiny App](#shiny-app)
  + [Example](#example)

## Dependencies

For the shiny app, the dependencies for rMATS-long are needed as well as the r-shiny package. The rMATS-long dependencies can be installed to a conda environment by running [../install](../install). Then the r-shiny package can be installed
```
cd rMATS-long
./install
conda activate ./conda_env
conda install -c conda-forge -c bioconda r-shiny
```

The shiny server can be run when the conda environment is activated

## Usage

### Data

The app will look in `data/` for datasets. Each subdirectory is treated as a separate dataset and each dataset should have these files:

* `group_1.txt`: file used as `--group-1` with `rmats_long.py`
* `group_2.txt`: file used as `--group-2` with `rmats_long.py`
* `abundance.esp`: file used as `--abundance` with `rmats_long.py`
* `updated.gtf`: file used as `--updated-gtf` with `rmats_long.py`
* `reference.gtf`: file used as `--gencode-gtf` with `rmats_long.py`
* `differential_genes.tsv`: output from `rmats_long.py`
* `differential_transcripts.tsv`: output from `rmats_long.py`

### Shiny App

Run the app with: `Rscript ./run.R`

Then the app can be used with a web browser on the same machine that is running the app with the URL `http://localhost:8888/`. The shiny documentation has more details: https://shiny.posit.co/r/reference/shiny/latest/runapp

With the app open in a browser the plots for a gene can be generated:
* Select a dataset from the dropdown
* Fill out the input boxes (either "gene name" or "gene ID" is required)
* Click the "create plots" button
* Wait for the plot to be displayed

After plotting the gene, the "Isoform differences" section at the bottom of the page can be run by setting a transcript ID and clicking the "classify isoform differences" button

### Example

After running the example from the main rMATS-long [README.md](../README.md)
```
mkdir -p data/example
cp ../example/group_1.txt ./data/example/
cp ../example/group_2.txt ./data/example/
cp ../example/samples_N2_R0_abundance.esp ./data/example/abundance.esp
cp ../example/samples_N2_R0_updated.gtf ./data/example/updated.gtf
cp ../example/gencode.v43.annotation_filtered.gtf ./data/example/reference.gtf
cp ../example_out/differential_genes.tsv ./data/example/
cp ../example_out/differential_transcripts.tsv ./data/example/
```

Follow the instructions at [Shiny App](#shiny-app). Gene name "DDR1" can be used
