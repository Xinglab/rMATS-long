# rMATS-long

## About

rMATS-long is a collection of tools for analyzing long-read data

## Table of Contents

* [Dependencies](#dependencies)
* [Usage](#usage)
  + [Detect Differential Isoforms](#detect-differential-isoforms)
  + [Visualize Isoforms and Abundance](#visualize-isoforms-and-abundance)
  + [Classify Isoform Differences](#classify-isoform-differences)
  + [Example](#example)

## Dependencies

Dependencies can be installed to a conda environment by running [./install](./install). Then the scripts can be run when the conda environment is activated: `conda activate ./conda_env`

* [Python](https://www.python.org/) (v3.10.9)
  + [NetworkX](https://networkx.org/) (v2.8.8)
  + [NumPy](https://numpy.org/) (v1.24.2)
  + [pandas](https://pandas.pydata.org/) (v1.5.3)
* [R](https://www.r-project.org/) (v4.2.2)
  + [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) (v1.32.5)
  + [DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) (v1.26.0)
  + [ggplot2](https://ggplot2.tidyverse.org/) (v3.4.0)

Those versions were used during testing. Other versions may also be compatible

## Usage

First run [ESPRESSO](https://github.com/Xinglab/espresso) to detect and quantify isoforms using long-read data. The ESPRESSO output can be used to:
* [Detect differential isoform usage](#detect-differential-isoforms): [scripts/detect_differential_isoforms.py](scripts/detect_differential_isoforms.py)
* [Visualize isoforms and abundance](#visualize-isoforms-and-abundance): [scripts/visualize_isoforms.py](scripts/visualize_isoforms.py)
* [Classify isoform differences](#classify-isoform-differences): [scripts/classify_isoform_differences.py](scripts/classify_isoform_differences.py)

### Detect Differential Isoforms

[scripts/detect_differential_isoforms.py](scripts/detect_differential_isoforms.py) detects differential isoform expression using DRIMSeq. The samples in the ESPRESSO abundance file need to be separated into two groups. The two groups are written to the `--group-1` and `--group-2` input files as comma separated lists. The isoform counts from the abundance file are then used in the DRIMSeq pipeline

The main output file written to `--out-dir` is `differential_transcripts.tsv` which gives an adjusted p-value for transcripts that show differential usage between the two groups. The `--out-dir` contains:
* `differential_transcripts.tsv`
* `differential_genes.tsv`
* `transcripts_per_gene.png`
* `precision_by_gene_expression.png`
* `transcript_pvalues.png`
* `gene_pvalues.png`

```
python detect_differential_isoforms.py -h

usage: detect_differential_isoforms.py [-h] --abundance ABUNDANCE --out-dir
                                       OUT_DIR --group-1 GROUP_1 --group-2
                                       GROUP_2 [--num-threads NUM_THREADS]

Detect differential isoform expression using DRIMSeq

options:
  -h, --help            show this help message and exit
  --abundance ABUNDANCE
                        The path to abundance.esp file from ESPRESSO
  --out-dir OUT_DIR     The path to use as output directory
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list. The sample names
                        should match with the ESPRESSO abundance column names.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2.
  --num-threads NUM_THREADS
                        The number of threads to use (default: 1)
```

### Visualize Isoforms and Abundance

TODO: [scripts/visualize_isoforms.py](scripts/visualize_isoforms.py)

### Classify Isoform Differences

[scripts/classify_isoform_differences.py](scripts/classify_isoform_differences.py) compares the structures of isoforms within a gene by calling [scripts/FindAltTSEvents.py](scripts/FindAltTSEvents.py) multiple times using a "main" isoform and all other isoforms in a gene

[scripts/FindAltTSEvents.py](scripts/FindAltTSEvents.py) compares the structures of any two transcript isoforms. Local differences in transcript structure are classified into 7 basic alternative splicing categories:

![basic alternative splicing patterns](docs/basic_alt_splicing_patterns.jpg)

* Exon skipping (SE)
* Alternative 5'-splice site (A5SS)
* Alternative 3'-splice site (A3SS)
* Mutually exclusive exons (MXE)
* Intron retention (RI)
* Alternative first exon (AFE)
* Alternative last exon (ALE)

Any local differences in transcript structure that could not be classified as one of the 7 basic alternative splicing categories are classified as "complex" (COMPLEX). **Note:** It is possible to have combinations of alternative splicing events for any given pair of transcript isoforms.

The output is a tab-delimited file consisting of four fields:
* **Field 1**: ID for transcript isoform 1
* **Field 2**: ID for transcript isoform 2
* **Field 3**: Discovered alternative splicing events
* **Field 4**: Genomic coordinates for alternative splicing events

**Note:** Designation of transcript isoforms 1 and 2 is completely arbitrary. Moreover, if the two transcript isoforms contained in the input GTF file exhibit a combination of multiple alternative splicing events, each event will be reported as its own line in the output file.

```
python classify_isoform_differences.py -h

usage: classify_isoform_differences.py [-h] --main-transcript-id
                                       MAIN_TRANSCRIPT_ID --gtf GTF --out-tsv
                                       OUT_TSV

Compare the structures of isoforms within a gene

options:
  -h, --help            show this help message and exit
  --main-transcript-id MAIN_TRANSCRIPT_ID
                        The transcript_id of the main isoform in the .gtf file
  --gtf GTF             The path to a .gtf describing the isoforms
  --out-tsv OUT_TSV     The path of the output file
```

```
python FindAltTSEvents.py -h

usage: FindAltTSEvents.py [-h] -i /path/to/input/GTF -o /path/to/output/file

This is a script to enumerate all transcript structure differences between a
pair of transcript isoforms

options:
  -h, --help            show this help message and exit
  -i /path/to/input/GTF
                        path to GTF file describing structures of two
                        transcript isoforms
  -o /path/to/output/file
                        path to output file
```

### Example

<!-- TODO data for 1 or 2 genes ? -->

[test_data/](test_data/) includes:
* [TODO.gtf](test_data/TODO.gtf)
* [TODO_abundance.esp](test_data/TODO_abundance.esp)
* [group_1.txt](test_data/group_1.txt)
* [group_2.txt](test_data/group_2.txt)

First Detect differential isoform usage. The sample names from the abundance file are split into the group files. `group_1.txt`:
```
sample_1,sample_2
```
and `group_2.txt`
```
sample_3,sample_4
```

```
python ./scripts/detect_differential_isoforms.py --abundance ./test_data/TODO_abundance.esp --out-dir ./differential_usage_out --group-1 ./test_data/group_1.txt --group-2 ./test_data/group_2.txt --num-threads 1
```

From `./differential_usage_out/differential_transcripts.tsv` the top row is:
```
gene_id	feature_id	lr	df	pvalue	adj_pvalue
TODO_gene	TODO_transcript	0.00	1	1.00	1.00
```
The `gene_id` and `feature_id` from that row can be used with [scripts/visualize_isoforms.py](scripts/visualize_isoforms.py).
<!-- TODO more details -->

The differences among transcripts within that gene can be determined with:
```
python ./scripts/classify_isoform_differences.py --main-transcript-id TODO --gtf ./test_data/TODO.gtf --out-tsv isoform_differences.tsv
```
