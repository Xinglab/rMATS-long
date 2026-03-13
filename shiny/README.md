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

For the shiny app, the dependencies for rMATS-long are needed as well as the `r-shiny` and `r-bslib` packages. The rMATS-long dependencies can be installed to a conda environment by running [../install](../install). Then the additional packages can be installed:
```
cd rMATS-long
./install
conda activate ./conda_env
conda install -c conda-forge -c bioconda --file ./shiny/conda_requirements.txt
```

The shiny server can be run when the conda environment is activated

## Usage

### Data

The app will prompt for an rMATS-long output directory. From that directory, it will attempt to find files and directories named like the output from the main rMATS-long README.md examples or like the snakemake workflow output. The required files can also be specified manually

### Shiny App

Run the app with: `Rscript ./run.R`

Then the app can be used with a web browser on the same machine that is running the app with the URL `http://localhost:8888/`. The shiny documentation has more details: https://shiny.posit.co/r/reference/shiny/latest/runapp

With the app open in a browser, plots of the results can be created and viewed

### Example

After running the "Isoform Analysis Example" from the main rMATS-long [README.md](../README.md) the output under `../example_out/` can be used

* Initially, the `Data` tab is selected
* Provide the full path to `/path/to/example_out/rmats_long/` in the `main-directory` field
* Click the `set-main-directory` button
  + Most of the remaining fields in the `Data` tab should fill in automatically
* Click on the `Summary` tab
* There are sub-tabs showing
  + A summary pie plot
  + The summary.txt values
  + A table of ASM statistical test results
  + A table of isoform statistical test results
* `DDR1` is one significant gene in the output
  + The significant `asm_id` is `2_0`
  + The `gene_id` is `ENSG00000204580.14`
* Click on the `ASM` tab
* Provide `2_0` in the `asm-id` field
* Click `find-gene-id-for-asm` to fill in the `gene-id` field
* There are sub-tabs showing results
* First, scroll down on the `Isoform abundance` sub-tab and click `create plot`
  + Creating that plot generates some files needed for the other tabs (like the `colors.tsv` field at the top of the `ASM` tab)
* Click on the `Isoform structure` sub-tab to see additional plots created along with the abundance plot
* Click the `Region within gene` sub-tab and scroll down and click `create plot`
* Click the `Splicing differences` sub-tab and click `create table`
* Click the `Simple splice graph` sub-tab and scroll down and click `create plot`
* Click the `Splice graph` sub-tab and scroll down and click `create plot`
  + This plot could take a long time to create for an ASM with a large number of possible splicing paths
