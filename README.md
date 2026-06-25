# rMATS-long

[![Latest Release](https://img.shields.io/github/release/Xinglab/rMATS-long.svg?label=Latest%20Release)](https://github.com/Xinglab/rMATS-long/releases/latest)
[![Total GitHub Downloads](https://img.shields.io/github/downloads/Xinglab/rMATS-long/total.svg?label=Total%20GitHub%20Downloads)](https://github.com/Xinglab/rMATS-long/releases)
[![Total Bioconda Installs](https://img.shields.io/conda/dn/bioconda/rmats-long.svg?label=Total%20Bioconda%20Installs)](https://anaconda.org/bioconda/rmats-long)
[![Total Docker Pulls](https://img.shields.io/docker/pulls/xinglab/rmats-long.svg?label=Total%20Docker%20Pulls)](https://hub.docker.com/r/xinglab/rmats-long)
[![PyPI Installs](https://img.shields.io/pypi/dm/rmats-long.svg?label=PyPI%20Installs)](https://pypi.org/project/rmats-long/)

## About

rMATS-long is an integrated computational workflow for long-read RNA-seq data. Using transcript definitions and read alignments, rMATS-long enables differential isoform analysis between sample groups, as well as classification and visualization of isoform structure and abundance.

## Table of Contents

* [Dependencies](#dependencies)
  + [Tests](#tests)
* [Usage](#usage)
  + [Differential Isoform Analysis](#differential-isoform-analysis)
    - [Starting From Quantified Isoforms](#starting-from-quantified-isoforms)
  + [Alternative Splicing Module Analysis](#alternative-splicing-module-analysis)
  + [Snakemake](#snakemake)
  + [Comparing Multiple Groups](#comparing-multiple-groups)
  + [Shiny](#shiny)
  + [Examples](#examples)
    - [Creating an input GTF](#creating-an-input-gtf)
    - [Isoform Analysis Example](#isoform-analysis-example)
    - [From Abundance Example](#from-abundance-example)
    - [ASM Analysis Example](#asm-analysis-example)
    - [Basic Events Example](#basic-events-example)
  + [Individual Scripts](#individual-scripts)
    - [rmats_long.py](#rmats_longpy)
    - [detect_differential_isoforms.py](#detect_differential_isoformspy)
    - [count_significant_isoforms.py](#count_significant_isoformspy)
    - [visualize_isoforms.py](#visualize_isoformspy)
    - [classify_isoform_differences.py](#classify_isoform_differencespy)
    - [FindAltTSEvents.py](#findalttseventspy)
    - [organize_gene_info_by_chr.py](#organize_gene_info_by_chrpy)
    - [simplify_alignment_info.py](#simplify_alignment_infopy)
    - [organize_alignment_info_by_gene_and_chr.py](#organize_alignment_info_by_gene_and_chrpy)
    - [detect_splicing_events.py](#detect_splicing_eventspy)
    - [count_reads_for_asms.py](#count_reads_for_asmspy)
    - [plot_splice_graph.py](#plot_splice_graphpy)
    - [plot_simple_splice_graph.py](#plot_simple_splice_graphpy)
    - [create_gtf_from_asm_definitions.py](#create_gtf_from_asm_definitionspy)
    - [run_stat_model.py](#run_stat_modelpy)
    - [check_asm_coupling.py](#check_asm_couplingpy)

## Dependencies

The dependencies and the rmats-long Python package can be installed to a conda environment by running [./install](./install).

Another option is to install the rmats-long bioconda package:
```
conda install -c conda-forge -c bioconda rmats-long
```

In either case, when the conda environment is activated the scripts can be run with `rmats-long` and the name of the script. For example: `rmats-long rmats_long.py --help`

The PyPI package can be installed with `pip`, but non-Python dependencies will need to be installed separately

* [Python](https://www.python.org/) (v3.11.15)
  + [NetworkX](https://networkx.org/) (v2.8.8)
  + [NumPy](https://numpy.org/) (v1.24.4)
  + [pandas](https://pandas.pydata.org/) (v2.0.3)
  + [matplotlib](https://matplotlib.org/) (v3.7.3)
  + [pydot](https://github.com/pydot/pydot) (v3.0.4)
  + [rpy2](https://github.com/rpy2/rpy2) (v3.6.7)
  + [threadpoolctl](https://github.com/joblib/threadpoolctl) (v3.6.0)
  + [Cython](https://cython.org/) (v3.2.6)
* [R](https://www.r-project.org/) (v4.5.3)
  + [BiocParallel](https://bioconductor.org/packages/release/bioc/html/BiocParallel.html) (v1.44.0)
  + [DRIMSeq](https://bioconductor.org/packages/release/bioc/html/DRIMSeq.html) (v1.38.0)
  + [ggplot2](https://ggplot2.tidyverse.org/) (v3.5.2)
  + [cowplot](https://wilkelab.org/cowplot/index.html) (v1.2.0)
  + [ggrepel](https://ggrepel.slowkow.com/) (v0.9.8)
  + [ggVennDiagram](https://github.com/gaospecial/ggVennDiagram) (v1.5.7)
  + [this.path](https://github.com/ArcadeAntics/this.path) (v2.7.1)
  + [mclogit](https://github.com/melff/mclogit/) (v0.9.15)
  + [emmeans](https://github.com/rvlenth/emmeans) (v2.0.3)
  + [lme4](https://github.com/lme4/lme4) (v1.1-38)
* [samtools](https://www.htslib.org/) (v1.22.1)
* [GCC](https://gcc.gnu.org/) (v15.2.0)
* Optional dependencies for creating a GTF with [AMALGAM](https://github.com/RNA-ROB/amalgam):
  + [gffcompare](https://github.com/gpertea/gffcompare) (v0.12.6)
  + [pysam](https://github.com/pysam-developers/pysam) (v0.23.3)
  + [pytabix](https://github.com/slowkow/pytabix) (v0.1)
  + [stringtie](https://github.com/gpertea/stringtie) (v3.0.3)

Those versions were used during testing. Other versions may also be compatible.

### Tests

With the dependencies installed (and the conda environment activated if needed), the automated tests can be run with: `./run_tests`

An individual test can be run by supplying the path to the `test.py` file: `./run_tests ./tests/se_gene/test.py`

## Usage

rMATS-long can analyze full-length isoforms or alternative splicing modules (ASMs). Full-length isoforms are quantified from long-read data using the code in [src/rmats_long/](src/rmats_long/) ([Differential Isoform Analysis](#differential-isoform-analysis)). Another option is to detect ASMs from the full-length isoforms and quantify those ASM isoforms ([Alternative Splicing Module Analysis](#alternative-splicing-module-analysis)). Significant changes in isoform usage are detected and visualizations are created. An annotation file with full-length isoforms is required as input. The isoform definitions can be from any source, for example [AMALGAM](https://github.com/RNA-ROB/amalgam). The analysis can be restricted to the basic splicing event types used in [rMATS turbo](https://github.com/Xinglab/rmats-turbo) ([Basic Events Example](#basic-events-example)).

### Differential Isoform Analysis

Running [rmats_long.py](#rmats_longpy) with full-length isoforms requires:
* Two groups of samples to compare
* A `.gtf` with isoform definitions
* Either
  + Alignments for each sample
  + Isoform counts for each sample ([Starting From Quantified Isoforms](#starting-from-quantified-isoforms))

The main outputs are:
* Summary information: `summary.txt`, `summary_plot.png`
* List of differential isoforms: `differential_isoforms_filtered.tsv`
* For each gene with a significant isoform:
  + Visualization of the isoform abundance by sample: `results_by_gene/{gene}/{id}_abundance.png`
  + Visualization of the isoform structures: `results_by_gene/{gene}/{id}_structure.png`
  + List of splicing changes between the most significant isoforms: `results_by_gene/{gene}/{id}_isoform_differences_{isoform}_to_{isoform}.tsv`

The main steps are:
* [organize_gene_info_by_chr.py](#organize_gene_info_by_chrpy)
* [simplify_alignment_info.py](#simplify_alignment_infopy)
* [organize_alignment_info_by_gene_and_chr.py](#organize_alignment_info_by_gene_and_chrpy)
* [detect_splicing_events.py](#detect_splicing_eventspy) (`--output-full-gene-asm`)
* [count_reads_for_asms.py](#count_reads_for_asmspy)
* [rmats_long.py](#rmats_longpy) (`--no-splice-graph-plot`)
  + [detect_differential_isoforms.py](#detect_differential_isoformspy)
  + [visualize_isoforms.py](#visualize_isoformspy)
  + [classify_isoform_differences.py](#classify_isoform_differencespy)

See [Isoform Analysis Example](#isoform-analysis-example) for details.

#### Starting From Quantified Isoforms

With isoform counts for each sample already available, [rmats_long.py](#rmats_longpy) can be run directly using the `--abundance` argument.

The main outputs are:
* Summary information: `summary.txt`, `summary_plot.png`
* List of differential isoforms: `differential_transcripts_filtered.tsv`
* For each gene with a significant isoform:
  + Visualization of the isoform abundance by sample: `results_by_gene/{gene}/{gene}_abundance.png`
  + Visualization of the isoform structures: `results_by_gene/{gene}/{gene}_structure.png`
  + List of splicing changes between the most significant isoforms: `results_by_gene/{gene}/{gene}_isoform_differences_{isoform}_to_{isoform}.tsv`

See [From Abundance Example](#from-abundance-example) for details.

### Alternative Splicing Module Analysis

Running [rmats_long.py](#rmats_longpy) with ASMs requires:
* Two groups of samples to compare
* A `.gtf` with isoform definitions
* Alignments for each sample

The main outputs are:
* Summary information: `summary.txt`, `summary_plot.png`
* ASM isoform definitions: `asm.gtf`
* List of differential isoforms: `differential_isoforms_filtered.tsv`
* For each ASM with a significant isoform:
  + Visualization of the isoform abundance by sample: `results_by_gene/{gene}/{asm}_abundance.png`
  + Visualization of the isoform structures: `results_by_gene/{gene}/{asm}_structure.png`
  + Visualization of where the ASM region is within the gene: `results_by_gene/{gene}/{asm}_in_gene.png`
  + Visualization of the ASM splice graph: `results_by_gene/{gene}/{asm}_graph.png`, `results_by_gene/{gene}/{asm}_simple_graph.png`
  + List of splicing changes between the most significant isoforms: `results_by_gene/{gene}/{asm}_isoform_differences_{isoform}_to_{isoform}.tsv`
* If multiple significant ASMs were detected with the same splicing changes between the most significant isoforms then there will be a file listing the "duplicate" ASMs for that gene: `results_by_gene/{gene}/duplicate_asms.tsv`
* Statistics for ASM distal coupling: `distal_coupling_stats.tsv`

The main steps are:
* [organize_gene_info_by_chr.py](#organize_gene_info_by_chrpy)
* [simplify_alignment_info.py](#simplify_alignment_infopy)
* [organize_alignment_info_by_gene_and_chr.py](#organize_alignment_info_by_gene_and_chrpy)
* [detect_splicing_events.py](#detect_splicing_eventspy)
* [create_gtf_from_asm_definitions.py](#create_gtf_from_asm_definitionspy)
* [count_reads_for_asms.py](#count_reads_for_asmspy)
* [rmats_long.py](#rmats_longpy)
  + [detect_differential_isoforms.py](#detect_differential_isoformspy)
  + [visualize_isoforms.py](#visualize_isoformspy)
  + [classify_isoform_differences.py](#classify_isoform_differencespy)
  + [plot_splice_graph.py](#plot_splice_graphpy)
  + [plot_simple_splice_graph.py](#plot_simple_splice_graphpy)

See [ASM Analysis Example](#asm-analysis-example) for details.

### Snakemake

`snakemake` can be used to run all the steps in the workflow. After setting the configuration in [snakemake_config.yaml](snakemake_config.yaml) the workflow can be run with:
```
./run_snakemake
```

The main required configuration values are:
* `run_name`: used to name output files
* `gtf_name`: the `.gtf` file to use
* `group_1_samples` and `group_2_samples:`: each entry gives a sample name and a list of `sam`, `bam`, or `simplified` alignment files.
* `quantify_full_length_transcripts`: whether to analyze full-length transcripts
* `quantify_asms`: whether to analyze ASM isoforms
* `quantify_basic_events`: whether to do an analysis restricted to basic splicing events (SE, A5SS, A3SS, MXE, RI)
* `require_asms_to_be_strict`: if set, then `--output-strict-only` will be used during the `quantify_asms` analysis

A Venn diagram of significant genes will be output if run with both `quantify_full_length_transcripts` and `quantify_asms`.

If a file with `gtf_name` is found in `references/` then it will be used. The gtf can also be configured with a URL as explained in the config file.

There are also values to allocate resources like `{job}_mem_gb`, `{job}_threads`, and `{job}_time_hr`.

The files in [snakemake_profile/](snakemake_profile/) are used to allow submitting jobs to an HPC system. The main file is [snakemake_profile/config.v8+.yaml](snakemake_profile/config.v8+.yaml) which sets the commands snakemake will use to interact with the compute cluster.

The snakemake commands are run using [./conda_wrapper](./conda_wrapper) to use the dependencies installed to the conda environment. The [install](install) script writes the `conda_wrapper:` path to [snakemake_config.yaml](snakemake_config.yaml).

The default config is set to run the example files. First unpack the files as in [Examples](#examples). Then set `gtf_name: 'gencode.v43.annotation_filtered.gtf'` in [snakemake_config.yaml](snakemake_config.yaml) and copy the `.gtf` to `references/`:
```
mkdir references
cp example/gencode.v43.annotation_filtered.gtf references/
```

Replace `/path/to/example` with the actual path to [./example](./example) on your system. Finally run `./run_snakemake` to produce output at `./results/example/`

### Comparing Multiple Groups

If the input samples are from more than two groups then the workflow can be run multiple times to compare different pairs of groups. As an example, if the samples are grouped into three different conditions:

* Detect splicing events using samples from all three groups
  + List all samples under `group_1_samples:`
  + Set `group_2_samples: {}` or remove `group_2_samples` from [snakemake_config.yaml](snakemake_config.yaml)
  + Set `run_name:` with a descriptive name such as `'all_group_run'`
  + Run `./run_snakemake`
* Perform a two group comparison on groups 1 and 2 using the previously detected events
  + List the samples for group 1 under `group_1_samples:` and the samples for group 2 under `group_2_samples:`
    - The files like `results/all_group_run/simplified/{sample}_{i}_simplified.tsv` can be configured with `simplified:` instead of the original `sam:` or `bam:`
  + Set `group_1_name:` and `group_2_name:`
  + Set `annotation_dir: '/path/to/results/all_group_run/annotation'`
  + Set `predefined_asm_events: '/path/to/results/all_group_run/asm_events'`
  + Set `predefined_basic_events: '/path/to/results/all_group_run/basic_events'`
  + Set `predefined_full_length_events: '/path/to/results/all_group_run/full_length_events'`
  + Set `run_name:` with a descriptive name such as `'group_1_vs_group_2'`
  + Run `./run_snakemake`
* Perform a two group comparison on groups 1 and 3 using the previously detected events
  + Like above, but update `group_1_samples:`, `group_2_samples:`, `group_1_name:`, `group_2_name:` and `run_name:` to reflect the current comparison
* Perform a two group comparison on groups 2 and 3 using the previously detected events
  + Like above, but update `group_1_samples:`, `group_2_samples:`, `group_1_name:`, `group_2_name:` and `run_name:` to reflect the current comparison

### Shiny

* [shiny/](shiny/) provides a web interface to view rMATS-long output files

### Examples

Example data is provided in [example/data.tar.gz](example/data.tar.gz). Unpack that file with:
```
cd example/
tar -xvf ./data.tar.gz
```

The unpacked files are:
* `example/gencode.v43.annotation_filtered.gtf`
* `example/GRCh38.primary_assembly.genome_filtered.fa`
* `example/group_1.txt`
* `example/group_2.txt`
* `example/gs689_1_filtered.sam`
* `example/gs689_2_filtered.sam`
* `example/gs689_3_filtered.sam`
* `example/pc3e_1_filtered.sam`
* `example/pc3e_2_filtered.sam`
* `example/pc3e_3_filtered.sam`
* `example/gs689_1_corrected.sam`
* `example/gs689_2_corrected.sam`
* `example/gs689_3_corrected.sam`
* `example/pc3e_1_corrected.sam`
* `example/pc3e_2_corrected.sam`
* `example/pc3e_3_corrected.sam`
* `example/samples_N2_R0_abundance.esp`
* `example/samples_N2_R0_updated.gtf`

The example data is based on cell line data from [https://doi.org/10.1126/sciadv.abq5072](https://doi.org/10.1126/sciadv.abq5072). The 1D cDNA sequencing for GS689 and PC3E was processed to get .sam files. The reference data (gencode .gtf and GRCh38 .fa) and the .sam files were filtered to a few different regions to get a small dataset. ESPRESSO was run to get the corrected alignments, isoform abundance, and updated .gtf

#### Creating an input GTF

A GTF with high-confidence transcripts can be created using [AMALGAM](https://github.com/RNA-ROB/amalgam). From that repo, `scripts/Build_Transcriptome.py`, `assets/human.refTSS_v4.1.hg38.bed.gz`, `assets/atlas.clusters.2.0.GRCh38.bed.gz`, and the `.tbi` index files for the `.bed.gz` files will be used

Sorted `.sam` or `.bam` files are needed. The example `*_filtered.sam` files are already sorted. Sorted alignment files can be created with a command like: `samtools sort -o sorted.bam unsorted.sam`

First run `stringtie` to get a `.gtf` for each input file:
```
stringtie -o example/gs689_1.gtf example/gs689_1_filtered.sam
stringtie -o example/gs689_2.gtf example/gs689_2_filtered.sam
stringtie -o example/gs689_3.gtf example/gs689_3_filtered.sam
stringtie -o example/pc3e_1.gtf example/pc3e_1_filtered.sam
stringtie -o example/pc3e_2.gtf example/pc3e_2_filtered.sam
stringtie -o example/pc3e_3.gtf example/pc3e_3_filtered.sam
```

Combine those `.gtf` files with the GENCODE annotation by creating `example/gtf_list.txt` with:
```
example/gencode.v43.annotation_filtered.gtf
example/gs689_1.gtf
example/gs689_2.gtf
example/gs689_3.gtf
example/pc3e_1.gtf
example/pc3e_2.gtf
example/pc3e_3.gtf
```

And then run:
```
gffcompare -i example/gtf_list.txt -T -o example/gffcompare
```

The command to create the GTF of high-confidence transcripts is:
```
python amalgam/scripts/Build_Transcriptome.py -i example/gffcompare -g example/gencode.v43.annotation_filtered.gtf -f example/GRCh38.primary_assembly.genome_filtered.fa -x amalgam/assets/human.refTSS_v4.1.hg38.bed.gz -y amalgam/assets/atlas.clusters.2.0.GRCh38.bed.gz -o example/combined.gtf
```

The `gene_name` and `Ensembl_canonical` tag are copied from the GENCODE GTF attributes to the new GTF:
```
rmats-long copy_gtf_attributes.py --gencode-gtf example/gencode.v43.annotation_filtered.gtf --other-gtf example/combined.gtf --out-gtf example/combined_with_attributes.gtf
```

#### Isoform Analysis Example

First create a directory of sorted annotation files based on the `combined_with_attributes.gtf` from [Creating an input GTF](#creating-an-input-gtf) (or a reference `.gtf`):
```
rmats-long organize_gene_info_by_chr.py --gtf ./example/combined_with_attributes.gtf --out-dir ./example_out/annotation
```

Next the read alignments can be processed. The input files can be `.sam` or `.bam` files. The example files are from `minimap2`, but the alignments could have been produced by another tool. Each file is processed with a command like:
```
rmats-long simplify_alignment_info.py --in-file ./example/gs689_1_filtered.sam --out-tsv ./example_out/gs689_1_simplified.tsv
```

The commands for the other files are:
```
rmats-long simplify_alignment_info.py --in-file ./example/gs689_2_filtered.sam --out-tsv ./example_out/gs689_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/gs689_3_filtered.sam --out-tsv ./example_out/gs689_3_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_1_filtered.sam --out-tsv ./example_out/pc3e_1_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_2_filtered.sam --out-tsv ./example_out/pc3e_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_3_filtered.sam --out-tsv ./example_out/pc3e_3_simplified.tsv
```

The simplified alignments can then be used to create a directory of sorted alignment files. The script requires a `.tsv` file listing each simplified alignment file along with its sample name. Create a file, `./example_out/samples.tsv`, with:
```
gs689_1	./example_out/gs689_1_simplified.tsv
gs689_2	./example_out/gs689_2_simplified.tsv
gs689_3	./example_out/gs689_3_simplified.tsv
pc3e_1	./example_out/pc3e_1_simplified.tsv
pc3e_2	./example_out/pc3e_2_simplified.tsv
pc3e_3	./example_out/pc3e_3_simplified.tsv
```

Then run:
```
rmats-long organize_alignment_info_by_gene_and_chr.py --gtf-dir ./example_out/annotation --out-dir ./example_out/alignments --samples-tsv ./example_out/samples.tsv
```

The alignments need to be checked against the annotation to determine the compatible isoforms for each alignment. In order to use the same isoform compatibility code as the ASM workflow, [detect_splicing_events.py](#detect_splicing_eventspy) is run with `--output-full-gene-asm` to format the isoforms for each gene:
```
rmats-long detect_splicing_events.py --gtf-dir ./example_out/annotation --align-dir ./example_out/alignments --out-dir ./example_out/events --output-full-gene-asm
```

Next determine the compatible isoforms for each read with:
```
rmats-long count_reads_for_asms.py --event-dir ./example_out/events --gtf-dir ./example_out/annotation --align-dir ./example_out/alignments --out-dir ./example_out/asm_counts
```

Finally [rmats_long.py](#rmats_longpy) is run to determine the significant isoforms and produce the final output files. It requires two sample groups to be defined as in `group_1.txt`:
```
pc3e_1,pc3e_2,pc3e_3
```

and `group_2.txt`:
```
gs689_1,gs689_2,gs689_3
```

`--no-splice-graph-plot` is used since the splice graph plot at the gene level can have hundreds of isoforms. Here is the main command:
```
rmats-long rmats_long.py --gtf-dir ./example_out/annotation --align-dir ./example_out/alignments --event-dir ./example_out/events --asm-counts-dir ./example_out/asm_counts --gencode-gtf ./example/combined_with_attributes.gtf --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --out-dir ./example_out/rmats_long --plot-file-type .png --no-splice-graph-plot
```

`rmats_long.py` will run other commands. For this example it first runs:
```
rmats-long detect_differential_isoforms.py --out-dir ./example_out/rmats_long --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --adj-pvalue 0.05 --delta-proportion 0.05 --num-threads 1 --min-isoform-reads 1 --min-cpm-per-asm 0 --sample-read-total-tsv ./example_out/alignments/sample_read_totals.tsv --limit-asm-to-top-n-isoforms 50 --average-reads-per-group 10 --gene-cpm-tsv ./example_out/alignments/sample_gene_cpm.tsv --asm-proportion-of-gene 0.05 --asm-counts-dir ./example_out/asm_counts
```

Along with other status messages, that command should print: `found 17 isoforms from 5 ASMs from 5 genes with adj_pvalue <= 0.05 and abs(delta_isoform_proportion) >= 0.05 and average reads per group >= 10.0 and ASM CPM >= 5.0% of gene CPM`. One significant row from `./example_out/rmats_long/differential_isoforms_filtered.tsv` is:
```
asm_id	gene_id	gene_name	isoform_id	lr	df	pvalue	adj_pvalue	pc3e_1_proportion	pc3e_2_proportion	pc3e_3_proportion	gs689_1_proportion	gs689_2_proportion	gs689_3_proportion	group_1_average_proportion	group_2_average_proportion	delta_isoform_proportion	pc3e_1_count	pc3e_2_count	pc3e_3_count	gs689_1_count	gs689_2_count	gs689_3_count	pc3e_1_cpm	pc3e_2_cpm	pc3e_3_cpm	gs689_1_cpm	gs689_2_cpm	gs689_3_cpm
0_3	ENSG00000198561.16	CTNND1	ENST00000358694.10	222.3	1	2.842e-50	8.384e-49	0	0	0	0.2553	0.3632	0.2641	0	0.2942	-0.2942	0	0	0	55.45	88.51	106.5	0	0	0	6.562e+04	9.621e+04	6.916e+04
```

Next it will run a command similar to what is below using some temporary files:
```
rmats-long visualize_isoforms.py --gene-id ENSG00000198561.16 --abundance ./example_out/rmats_long/rmats_long_tmp/0_3_abun.tsv --updated-gtf ./example_out/rmats_long/rmats_long_tmp/0_3.gtf --diff-transcripts ./example_out/rmats_long/rmats_long_tmp/diff.tsv --out-dir ./example_out/rmats_long/results_by_gene/ENSG00000198561.16 --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --plot-file-type .png --intron-scaling 1 --max-transcripts 7 --gene-name CTNND1 --is-asm --graph-file ./example_out/events/graph_0.txt --asm-id 0_3 --out-transcript-colors ./example_out/rmats_long/rmats_long_tmp/colors.tsv
```

And produce `./example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_abundance.png`:
![0_3_abundance](example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_abundance.png)

And `./example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_structure.png`:
![0_3_structure](example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_structure.png)

The plots show that ENST00000358694.10 is abundant in GS689 but not in PC3E. ENST00000529986.5 has reads in both groups, but is more abundant in PC3E. ENST00000358694.10 is the most significant isoform for this gene and ENST00000529986.5 is the most significant isoform that has a delta proportion in the opposite direction of ENST00000358694.10.

The differences between those two selected isoforms are determined with:
```
rmats-long classify_isoform_differences.py --updated-gtf ./example_out/rmats_long/rmats_long_tmp/0_3.gtf --out-tsv ./example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_isoform_differences_ENST00000358694.10_to_ENST00000529986.5.tsv --main-transcript-id ENST00000358694.10 --second-transcript-id ENST00000529986.5
```

`./example_out/rmats_long/results_by_gene/ENSG00000198561.16/0_3_isoform_differences_ENST00000358694.10_to_ENST00000529986.5.tsv` shows that the difference is due to consecutive exons skipping. The difference is classified as `COMPLEX` and the coordinates are provided. This can also be seen in the isoform structure plot:
```
transcript1	transcript2	event	coordinates
ENST00000358694.10	ENST00000529986.5	COMPLEX	chr11:57762046:57762046:+;chr11:57794010:57794010:+;chr11:57762046:57762046:+;chr11:57789037:57789155:+;chr11:57791385:57791673:+;chr11:57794010:57794010:+
```

Similar commands are run for the other significant genes. A summary is written to `./example_out/rmats_long/summary.txt`:
```
## [...]/python rmats_long.py --gtf-dir ./example_out/annotation [...]
## YYYY-mm-ddTHH:MM:SS
## source code commit: [...]
## significant differential isoform usage
total significant isoforms: 17
total genes with significant isoforms: 5
total ASMs with significant isoforms: 5
adjusted pvalue threshold: 0.05
delta isoform proportion threshold: 0.05
## alternative splicing classifications between isoform pairs
total classified isoform pairs: 5
exon skipping: 2
alternative 5'-splice site: 0
alternative 3'-splice site: 0
mutually exclusive exons: 0
intron retention: 0
alternative first exon: 1
alternative last exon: 0
complex: 1
combinatorial: 1
alternative endpoints: 0
## Number of isoforms per ASM
total ASMs with 2 isoforms: 0
total ASMs with 3 isoforms: 0
total ASMs with 4 isoforms: 1
total ASMs with 5 isoforms: 0
total ASMs with 6 isoforms: 1
[...]
```

#### From Abundance Example

If isoform abundance by sample is already available then [rmats_long.py](#rmats_longpy) can be run without quantifying the isoform abundance again. The example data includes isoform definitions and abundance values from [ESPRESSO](https://github.com/Xinglab/espresso). The sample names from the abundance file need to be split into groups as done in `group_1.txt`:
```
pc3e_1,pc3e_2,pc3e_3
```

and `group_2.txt`:
```
gs689_1,gs689_2,gs689_3
```

Here is the main command:
```
rmats-long rmats_long.py --abundance ./example/samples_N2_R0_abundance.esp --updated-gtf ./example/samples_N2_R0_updated.gtf --gencode-gtf ./example/gencode.v43.annotation_filtered.gtf --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --out-dir ./example_out_from_abun --plot-file-type .png
```

`rmats_long.py` will run other commands. For this example it first runs:
```
rmats-long detect_differential_isoforms.py --out-dir ./example_out_from_abun --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --adj-pvalue 0.05 --delta-proportion 0.05 --num-threads 1 --min-isoform-reads 1 --min-cpm-per-asm 0 --sample-read-total-tsv ./example_out_from_abun/sample_read_totals.tsv --limit-asm-to-top-n-isoforms 50 --average-reads-per-group 10 --gene-cpm-tsv ./example_out_from_abun/sample_gene_cpm.tsv --asm-proportion-of-gene 0.05 --abundance ./example/samples_N2_R0_abundance.esp
```

Along with other status messages, that command should print: `found 8 isoforms from 3 genes with adj_pvalue <= 0.05 and abs(delta_isoform_proportion) >= 0.05 and average reads per group >= 10.0 and ASM CPM >= 5.0% of gene CPM`. One significant row from `./example_out_from_abun/differential_transcripts_filtered.tsv` is:
```
gene_id	feature_id	lr	df	pvalue	adj_pvalue	pc3e_1_proportion	pc3e_2_proportion	pc3e_3_proportion	gs689_1_proportion	gs689_2_proportion	gs689_3_proportion	group_1_average_proportion	group_2_average_proportion	delta_isoform_proportion	pc3e_1_count	pc3e_2_count	pc3e_3_count	gs689_1_count	gs689_2_count	gs689_3_count	pc3e_1_cpm	pc3e_2_cpm	pc3e_3_cpm	gs689_1_cpm	gs689_2_cpm	gs689_3_cpm
ENSG00000204580.14	ENST00000418800.6	101.1	1	8.619e-24	3.62e-23	0.2363	0.3629	0.2497	0	0	0.01065	0.283	0.003551	0.2794	16.54	26.49	22	0	0	1.31	8.148e+04	1.311e+05	1.106e+05	0	0	2319
```

Next it will run a command similar to what is below using some temporary files:
```
rmats-long visualize_isoforms.py --gene-id ENSG00000204580.14 --abundance ./example_out_from_abun/rmats_long_tmp/gene_abundance.esp --updated-gtf ./example_out_from_abun/rmats_long_tmp/gene_updated.gtf --diff-transcripts ./example_out_from_abun/rmats_long_tmp/gene_diff_transcripts.tsv --out-dir ./example_out_from_abun/results_by_gene/ENSG00000204580.14 --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --plot-file-type .png --intron-scaling 1 --max-transcripts 7 --gencode-gtf ./example_out_from_abun/rmats_long_tmp/gene_gencode.gtf
```

And produce `./example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_abundance.png`:
![ENSG00000204580.14_abundance](example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_abundance.png)

And `./example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_structure.png`:
![ENSG00000204580.14_structure](example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_structure.png)

The plots show that ENST00000418800.6 is abundant in PC3E but only has a few reads in GS689. ENST00000376568.8 has reads in both groups, but is more abundant in GS689. ENST00000418800.6 is the most significant isoform for this gene and ENST00000376568.8 is the most significant isoform that has a delta proportion in the opposite direction of ENST00000418800.6.

The differences between those two selected isoforms are determined with:
```
rmats-long classify_isoform_differences.py --updated-gtf ./example_out_from_abun/rmats_long_tmp/gene_updated.gtf --out-tsv ./example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_isoform_differences_ENST00000418800.6_to_ENST00000376568.8.tsv --main-transcript-id ENST00000418800.6 --second-transcript-id ENST00000376568.8 --gencode-gtf ./example_out_from_abun/rmats_long_tmp/gene_gencode.gtf
```

`./example_out_from_abun/results_by_gene/ENSG00000204580.14/ENSG00000204580.14_isoform_differences_ENST00000418800.6_to_ENST00000376568.8.tsv` shows that those two isoforms differ by an alternative first exon and also a skipped exon. This can also be seen in the isoform structure plot:
```
transcript1	transcript2	event	coordinates
ENST00000418800.6	ENST00000376568.8	AFE	chr6:30882613:30882983:+;chr6:30884519:30884710:+
ENST00000418800.6	ENST00000376568.8	SE	chr6:30895404:30895514:+
```

Similar commands are run for the other significant genes. A summary is written to `./example_out_from_abun/summary.txt`:
```
## [...]/python rmats_long.py --abundance ./example/samples_N2_R0_abundance.esp [...]
## YYYY-mm-ddTHH:MM:SS
## source code commit: [...]
## significant differential isoform usage
total significant isoforms: 8
total genes with significant isoforms: 3
adjusted pvalue threshold: 0.05
delta isoform proportion threshold: 0.05
## alternative splicing classifications between isoform pairs
total classified isoform pairs: 3
exon skipping: 0
alternative 5'-splice site: 0
alternative 3'-splice site: 0
mutually exclusive exons: 0
intron retention: 0
alternative first exon: 0
alternative last exon: 0
complex: 1
combinatorial: 2
alternative endpoints: 0
```

#### ASM Analysis Example

The workflow will search for ASMs within each gene. An ASM is a set of isoforms where all isoforms have the same start node and the same end node in the splice graph while no other node is shared by all isoforms. The isoforms do not need to be full length transcripts. If `--output-strict-only` is used, then each ASM is required not to overlap any splice graph edges which are not used in the ASM.

First create a directory of sorted annotation files based on the `combined_with_attributes.gtf` from [Creating an input GTF](#creating-an-input-gtf) (or a reference `.gtf`):
```
rmats-long organize_gene_info_by_chr.py --gtf ./example/combined_with_attributes.gtf --out-dir ./example_out_asm/annotation
```

Next the read alignments can be processed. The input files can be `.sam` or `.bam` files. The example files are from `minimap2`, but the alignments could have been produced by another tool. Each file is processed with a command like:
```
rmats-long simplify_alignment_info.py --in-file ./example/gs689_1_filtered.sam --out-tsv ./example_out_asm/gs689_1_simplified.tsv
```

The commands for the other files are:
```
rmats-long simplify_alignment_info.py --in-file ./example/gs689_2_filtered.sam --out-tsv ./example_out_asm/gs689_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/gs689_3_filtered.sam --out-tsv ./example_out_asm/gs689_3_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_1_filtered.sam --out-tsv ./example_out_asm/pc3e_1_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_2_filtered.sam --out-tsv ./example_out_asm/pc3e_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_3_filtered.sam --out-tsv ./example_out_asm/pc3e_3_simplified.tsv
```

The simplified alignments can then be used to create a directory of sorted alignment files. The script requires a `.tsv` file listing each simplified alignment file along with its sample name. Create a file, `./example_out_asm/samples.tsv`, with:
```
gs689_1	./example_out_asm/gs689_1_simplified.tsv
gs689_2	./example_out_asm/gs689_2_simplified.tsv
gs689_3	./example_out_asm/gs689_3_simplified.tsv
pc3e_1	./example_out_asm/pc3e_1_simplified.tsv
pc3e_2	./example_out_asm/pc3e_2_simplified.tsv
pc3e_3	./example_out_asm/pc3e_3_simplified.tsv
```

Then run:
```
rmats-long organize_alignment_info_by_gene_and_chr.py --gtf-dir ./example_out_asm/annotation --out-dir ./example_out_asm/alignments --samples-tsv ./example_out_asm/samples.tsv
```

With the `annotation/` and `alignments/` directories ready, the ASMs present in the data can be detected with:
```
rmats-long detect_splicing_events.py --gtf-dir ./example_out_asm/annotation --align-dir ./example_out_asm/alignments --out-dir ./example_out_asm/events
```

A `.gtf` file with the ASM definitions can be created with:
```
rmats-long create_gtf_from_asm_definitions.py --event-dir ./example_out_asm/events --out-gtf ./example_out_asm/asm.gtf
```

Next the alignments are checked against the ASM definitions to determine the compatible isoforms:
```
rmats-long count_reads_for_asms.py --event-dir ./example_out_asm/events --gtf-dir ./example_out_asm/annotation --align-dir ./example_out_asm/alignments --out-dir ./example_out_asm/asm_counts
```

Finally [rmats_long.py](#rmats_longpy) is run to determine the significant ASM isoforms and produce the final output files. It requires two sample groups to be defined as in `group_1.txt`:
```
pc3e_1,pc3e_2,pc3e_3
```

and `group_2.txt`:
```
gs689_1,gs689_2,gs689_3
```

Here is the main command:
```
rmats-long rmats_long.py --gtf-dir ./example_out_asm/annotation --align-dir ./example_out_asm/alignments --event-dir ./example_out_asm/events --asm-counts-dir ./example_out_asm/asm_counts --gencode-gtf ./example/combined_with_attributes.gtf --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --out-dir ./example_out_asm/rmats_long --plot-file-type .png
```

`rmats_long.py` will run other commands. For this example it first runs:
```
rmats-long detect_differential_isoforms.py --out-dir ./example_out_asm/rmats_long --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --adj-pvalue 0.05 --delta-proportion 0.05 --num-threads 1 --min-isoform-reads 1 --min-cpm-per-asm 0 --sample-read-total-tsv ./example_out_asm/alignments/sample_read_totals.tsv --limit-asm-to-top-n-isoforms 50 --average-reads-per-group 10 --gene-cpm-tsv ./example_out_asm/alignments/sample_gene_cpm.tsv --asm-proportion-of-gene 0.05 --asm-counts-dir ./example_out_asm/asm_counts
```

Along with other status messages, that command should print: `found 66 isoforms from 19 ASMs from 3 genes with adj_pvalue <= 0.05 and abs(delta_isoform_proportion) >= 0.05 and average reads per group >= 10.0 and ASM CPM >= 5.0% of gene CPM`. One significant row from `./example_out_asm/rmats_long/differential_isoforms_filtered.tsv` is:
```
asm_id	gene_id	gene_name	isoform_id	lr	df	pvalue	adj_pvalue	pc3e_1_proportion	pc3e_2_proportion	pc3e_3_proportion	gs689_1_proportion	gs689_2_proportion	gs689_3_proportion	group_1_average_proportion	group_2_average_proportion	delta_isoform_proportion	pc3e_1_count	pc3e_2_count	pc3e_3_count	gs689_1_count	gs689_2_count	gs689_3_count	pc3e_1_cpm	pc3e_2_cpm	pc3e_3_cpm	gs689_1_cpm	gs689_2_cpm	gs689_3_cpm
0_9	ENSG00000198561.16	CTNND1	0_9_3	242.9	1	9.22e-55	1.899e-52	1	1	1	0.25	0.2117	0.305	1	0.2556	0.7444	36.97	33.97	53.92	18.75	18.21	45.44	6.213e+04	7.686e+04	9.733e+04	2.219e+04	1.979e+04	2.951e+04
```

Next it will run a command similar to what is below using some temporary files:
```
rmats-long visualize_isoforms.py --gene-id ENSG00000198561.16 --abundance ./example_out_asm/rmats_long/rmats_long_tmp/0_9_abun.tsv --updated-gtf ./example_out_asm/rmats_long/rmats_long_tmp/0_9.gtf --diff-transcripts ./example_out_asm/rmats_long/rmats_long_tmp/diff.tsv --out-dir ./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16 --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --plot-file-type .png --intron-scaling 1 --max-transcripts 7 --gene-name CTNND1 --is-asm --graph-file ./example_out_asm/events/graph_0.txt --asm-id 0_9 --out-transcript-colors ./example_out_asm/rmats_long/rmats_long_tmp/colors.tsv
```

And produce `./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_abundance.png`:
![0_9_abundance](example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_abundance.png)

And `./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_structure.png`:
![0_9_structure](example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_structure.png)

The plots show that `0_9_3` which skips both exons is the only isoform in PC3E samples. In GS689 samples the most abundant isoform is `0_9_0` which includes both exons. `0_9_3` is the most significant isoform for this gene and `0_9_0` is the most significant isoform that has a delta proportion in the opposite direction of `0_9_3`.

Another [visualize_isoforms.py](#visualize_isoformspy) command will be run to produce `./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_in_gene.png` which shows where the ASM region is located within the full-length isoforms of the gene:
![0_9_in_gene](example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_in_gene.png)

The differences between those two significant isoforms are determined with:
```
rmats-long classify_isoform_differences.py --updated-gtf ./example_out_asm/rmats_long/rmats_long_tmp/0_9.gtf --out-tsv ./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_isoform_differences_0_9_3_to_0_9_0.tsv --main-transcript-id 0_9_3 --second-transcript-id 0_9_0
```

`./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_isoform_differences_0_9_3_to_0_9_0.tsv` shows that this is classified as `COMPLEX` (because of the consecutive exon skipping) and lists the exon coordinates:
```
transcript1	transcript2	event	coordinates
0_9_3	0_9_0	COMPLEX	chr11:57762046:57762046:+;chr11:57794010:57794010:+;chr11:57762046:57762046:+;chr11:57789037:57789155:+;chr11:57791385:57791673:+;chr11:57794010:57794010:+
```

Other significant ASMs detected in ENSG00000198561.16 have the same difference between their selected significant isoforms. `example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/duplicate_asms.tsv` has one line for each set of ASMs that have the same isoform differences. The line for this consecutive exon skipping is:
```
0_9	0_13	0_8	0_7	0_12	0_11	0_15	0_16
```

`0_9` is the first ASM in that row because it has the fewest isoforms. Only `0_9` from that row will be used when producing the values in `summary.txt`

A diagram of the splice graph will also be produced with:
```
rmats-long plot_splice_graph.py --event-dir ./example_out_asm/events --chr chr11 --gene-id ENSG00000198561.16 --asm-id 0_9 --out-file ./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_graph.png
```
![0_9_graph](example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_graph.png)

Another splice graph will be produced with:
```
rmats-long plot_simple_splice_graph.py --out-file ./example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_simple_graph.png --gtf ./example_out_asm/rmats_long/rmats_long_tmp/0_9.gtf --intron-scaling 1 --transcript-colors ./example_out_asm/rmats_long/rmats_long_tmp/colors.tsv --color-by-isoform --ids-from-color-file --show-isoform-endpoint-symbols
```
![0_9_simple_graph](example_out_asm/rmats_long/results_by_gene/ENSG00000198561.16/0_9_simple_graph.png)

Similar commands are run for the other significant ASMs. A summary is written to `./example_out_asm/rmats_long/summary.txt`:
```
## [...]/python rmats_long.py --gtf-dir ./example_out_asm/annotation [...]
## YYYY-mm-ddTHH:MM:SS
## source code commit: [...]
## significant differential isoform usage
total significant isoforms: 22
total genes with significant isoforms: 3
total ASMs with significant isoforms: 9
adjusted pvalue threshold: 0.05
delta isoform proportion threshold: 0.05
## alternative splicing classifications between isoform pairs
total classified isoform pairs: 9
exon skipping: 4
alternative 5'-splice site: 0
alternative 3'-splice site: 0
mutually exclusive exons: 0
intron retention: 0
alternative first exon: 2
alternative last exon: 0
complex: 3
combinatorial: 0
alternative endpoints: 0
## Number of isoforms per ASM
total ASMs with 2 isoforms: 2
total ASMs with 3 isoforms: 3
total ASMs with 4 isoforms: 1
total ASMs with 5 isoforms: 2
total ASMs with 6 isoforms: 0
total ASMs with 7 isoforms: 0
total ASMs with 8 isoforms: 1
```

An analysis of ASM coupling can be run with these two commands:
```
rmats-long check_asm_coupling.py --event-dir ./example_out_asm/events --asm-counts-dir ./example_out_asm/asm_counts --count-tsv ./example_out_asm/rmats_long/count.tsv --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --out-tsv ./example_out_asm/distal_coupling_counts.tsv --isoform-out-tsv ./example_out_asm/distal_coupling_main_isoforms.tsv
```

```
rmats-long check_asm_coupling.R ./example_out_asm/distal_coupling_counts.tsv ./example_out_asm/distal_coupling_stats.tsv
```

#### Basic Events Example

The workflow can identify basic splicing events instead of more complex ASMs. The procedure is the same as [ASM Analysis Example](#asm-analysis-example) except that `--output-basic-events` is used with [detect_splicing_events.py](#detect_splicing_eventspy). If using the snakemake, basic events can be detected with `quantify_basic_events: true` set in the config

Create a file, `./example_out_basic/samples.tsv`, with:
```
gs689_1	./example_out_basic/gs689_1_simplified.tsv
gs689_2	./example_out_basic/gs689_2_simplified.tsv
gs689_3	./example_out_basic/gs689_3_simplified.tsv
pc3e_1	./example_out_basic/pc3e_1_simplified.tsv
pc3e_2	./example_out_basic/pc3e_2_simplified.tsv
pc3e_3	./example_out_basic/pc3e_3_simplified.tsv
```

Run these commands
```
rmats-long organize_gene_info_by_chr.py --gtf ./example/combined_with_attributes.gtf --out-dir ./example_out_basic/annotation
```
```
rmats-long simplify_alignment_info.py --in-file ./example/gs689_1_filtered.sam --out-tsv ./example_out_basic/gs689_1_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/gs689_2_filtered.sam --out-tsv ./example_out_basic/gs689_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/gs689_3_filtered.sam --out-tsv ./example_out_basic/gs689_3_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_1_filtered.sam --out-tsv ./example_out_basic/pc3e_1_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_2_filtered.sam --out-tsv ./example_out_basic/pc3e_2_simplified.tsv
rmats-long simplify_alignment_info.py --in-file ./example/pc3e_3_filtered.sam --out-tsv ./example_out_basic/pc3e_3_simplified.tsv
```
```
rmats-long organize_alignment_info_by_gene_and_chr.py --gtf-dir ./example_out_basic/annotation --out-dir ./example_out_basic/alignments --samples-tsv ./example_out_basic/samples.tsv
```
```
rmats-long detect_splicing_events.py --gtf-dir ./example_out_basic/annotation --align-dir ./example_out_basic/alignments --out-dir ./example_out_basic/events --output-basic-events
```
```
rmats-long create_gtf_from_asm_definitions.py --event-dir ./example_out_basic/events --out-gtf ./example_out_basic/basic_events.gtf
```
```
rmats-long count_reads_for_asms.py --event-dir ./example_out_basic/events --gtf-dir ./example_out_basic/annotation --align-dir ./example_out_basic/alignments --out-dir ./example_out_basic/asm_counts
```
```
rmats-long rmats_long.py --gtf-dir ./example_out_basic/annotation --align-dir ./example_out_basic/alignments --event-dir ./example_out_basic/events --asm-counts-dir ./example_out_basic/asm_counts --gencode-gtf ./example/combined_with_attributes.gtf --group-1 ./example/group_1.txt --group-2 ./example/group_2.txt --group-1-name PC3E --group-2-name GS689 --out-dir ./example_out_basic/rmats_long --plot-file-type .png
```

Along with other status messages, the `rmats_long.py` command should print: `found 14 isoforms from 7 ASMs from 2 genes with adj_pvalue <= 0.05 and abs(delta_isoform_proportion) >= 0.05 and average reads per group >= 10.0 and ASM CPM >= 5.0% of gene CPM`. One significant row from `./example_out_basic/rmats_long/differential_isoforms_filtered.tsv` is:
```
asm_id	gene_id	gene_name	isoform_id	lr	df	pvalue	adj_pvalue	pc3e_1_proportion	pc3e_2_proportion	pc3e_3_proportion	gs689_1_proportion	gs689_2_proportion	gs689_3_proportion	group_1_average_proportion	group_2_average_proportion	delta_isoform_proportion	pc3e_1_count	pc3e_2_count	pc3e_3_count	gs689_1_count	gs689_2_count	gs689_3_count	pc3e_1_cpm	pc3e_2_cpm	pc3e_3_cpm	gs689_1_cpm	gs689_2_cpm	gs689_3_cpm
0_9	ENSG00000198561.16	CTNND1	0_9_0	89.5	1	3.065e-21	4.597e-20	0.02703	0.0294	0.01962	0.55	0.5385	0.411	0.02535	0.4998	-0.4745	1.027	1.029	1.079	23.65	21.54	33.29	1726	2328	1948	2.799e+04	2.341e+04	2.162e+04
```

The plots for that event are:
`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_abundance.png`:
![0_9_abundance](example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_abundance.png)

`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_structure.png`:
![0_9_structure](example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_structure.png)

`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_in_gene.png`:
![0_9_in_gene](example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_in_gene.png)

`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_graph.png`
![0_9_graph](example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_graph.png)

`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_simple_graph.png`
![0_9_simple_graph](example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_simple_graph.png)

`./example_out_basic/rmats_long/results_by_gene/ENSG00000198561.16/0_9_isoform_differences_0_9_0_to_0_9_1.tsv` shows that this is classified as exon skipping:
```
transcript1	transcript2	event	coordinates
0_9_0	0_9_1	SE	chr11:57791385:57791673:+
```

A summary is written to `./example_out_basic/rmats_long/summary.txt`:
```
## [...]/python rmats_long.py --gtf-dir ./example_out_basic/annotation [...]
## YYYY-mm-ddTHH:MM:SS
## source code commit: [...]
## significant differential isoform usage
total significant isoforms: 12
total genes with significant isoforms: 2
total ASMs with significant isoforms: 6
adjusted pvalue threshold: 0.05
delta isoform proportion threshold: 0.05
## alternative splicing classifications between isoform pairs
total classified isoform pairs: 6
exon skipping: 4
alternative 5'-splice site: 0
alternative 3'-splice site: 0
mutually exclusive exons: 1
intron retention: 1
alternative first exon: 0
alternative last exon: 0
complex: 0
combinatorial: 0
alternative endpoints: 0
## Number of isoforms per ASM
total ASMs with 2 isoforms: 6
```

### Individual Scripts

#### rmats_long.py

```
rmats-long rmats_long.py -h

usage: rmats_long.py [-h] [--abundance ABUNDANCE] [--updated-gtf UPDATED_GTF]
                     [--use-drimseq] [--event-dir EVENT_DIR]
                     [--asm-counts-dir ASM_COUNTS_DIR] [--align-dir ALIGN_DIR]
                     [--gtf-dir GTF_DIR] --group-1 GROUP_1 --group-2 GROUP_2
                     --out-dir OUT_DIR [--gencode-gtf GENCODE_GTF]
                     [--group-1-name GROUP_1_NAME]
                     [--group-2-name GROUP_2_NAME] [--num-threads NUM_THREADS]
                     [--process-top-n PROCESS_TOP_N]
                     [--process-selected PROCESS_SELECTED]
                     [--plot-file-type {.pdf,.png,all}]
                     [--intron-scaling INTRON_SCALING]
                     [--max-transcripts MAX_TRANSCRIPTS]
                     [--diff-transcripts DIFF_TRANSCRIPTS]
                     [--adj-pvalue ADJ_PVALUE] [--use-unadjusted-pvalue]
                     [--delta-proportion DELTA_PROPORTION]
                     [--compare-all-within-gene] [--covar-tsv COVAR_TSV]
                     [--min-cpm-per-asm MIN_CPM_PER_ASM]
                     [--no-splice-graph-plot]
                     [--min-isoform-reads MIN_ISOFORM_READS]
                     [--limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS]
                     [--average-reads-per-group AVERAGE_READS_PER_GROUP]
                     [--average-cpm-per-group AVERAGE_CPM_PER_GROUP]
                     [--min-cpm-per-group MIN_CPM_PER_GROUP]
                     [--asm-proportion-of-gene ASM_PROPORTION_OF_GENE]

Identify significant splicing changes and produce plots

options:
  -h, --help            show this help message and exit

Gene isoforms:
  Use either Gene isoforms or ASM isoforms

  --abundance ABUNDANCE
                        The path to the abundance.esp file from ESPRESSO
  --updated-gtf UPDATED_GTF
                        The path to the updated.gtf file from ESPRESSO
  --use-drimseq         Use DRIMSeq instead of run_stat_model.py

ASM isoforms:
  Use either ASM isoforms or Gene isoforms

  --event-dir EVENT_DIR
                        The output directory from detect_splicing_events.py
  --asm-counts-dir ASM_COUNTS_DIR
                        The output directory from count_reads_for_asms.py
  --align-dir ALIGN_DIR
                        The output directory from
                        organize_alignment_info_by_gene_and_chr.py
  --gtf-dir GTF_DIR     The output directory from organize_gene_info_by_chr.py

Required:
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list. The sample names
                        should match with the ESPRESSO abundance column names
                        or --asm-counts-dir names.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2
  --out-dir OUT_DIR     The path to use as the output directory

Optional:
  --gencode-gtf GENCODE_GTF
                        The path to a gencode annotation.gtf file. Will be
                        used to identify the Ensembl canonical isoform and the
                        gene name
  --group-1-name GROUP_1_NAME
                        A name for group 1 (default group 1)
  --group-2-name GROUP_2_NAME
                        A name for group 2 (default group 2)
  --num-threads NUM_THREADS
                        The number of threads to use (default 1)
  --process-top-n PROCESS_TOP_N
                        Generate plots and classify isoform differences for
                        the top "n" significant genes. By default all
                        significant genes are processed
  --process-selected PROCESS_SELECTED
                        A comma separated list of gene IDs or ASM IDs to
                        generate plots and classify isoform differences for
  --plot-file-type {.pdf,.png,all}
                        The file type for output plots (default .png))
  --intron-scaling INTRON_SCALING
                        The factor to use to reduce intron length in the plot.
                        A value of 2 would reduce introns to 1/2 of the
                        original plot length (default 1)
  --max-transcripts MAX_TRANSCRIPTS
                        How many transcripts to plot individually. The
                        remaining transcripts in the gene will be grouped
                        together (default 7)
  --diff-transcripts DIFF_TRANSCRIPTS
                        The path to the differential transcript results. If
                        given then skip the differential isoform calculation.
  --adj-pvalue ADJ_PVALUE
                        The cutoff for adjusted p-value (default 0.05)
  --use-unadjusted-pvalue
                        Use pvalue instead of adj_pvalue for the cutoff
  --delta-proportion DELTA_PROPORTION
                        The cutoff for delta isoform proportion (default 0.05)
  --compare-all-within-gene
                        Compare the most significant isoform to all other
                        isoforms in the gene. By default, the most significant
                        isoform is only compared to the most significant
                        isoform with a delta proportion in the opposite
                        direction.
  --covar-tsv COVAR_TSV
                        A .tsv with 1 line per sample. The first line has the
                        column names. The first column is sample_id. Each
                        additional column is a covariate.
  --min-cpm-per-asm MIN_CPM_PER_ASM
                        Only consider ASMs where at least 1 sample has at
                        least this CPM of reads assigned to the ASM. (default
                        0)
  --no-splice-graph-plot
                        Do not run plot_splice_graph.py
  --min-isoform-reads MIN_ISOFORM_READS
                        Only consider isoforms with at least this many reads
                        (default 1)
  --limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS
                        Only consider the top N isoforms with the highest
                        total proportion across samples for each ASM (default
                        50)
  --average-reads-per-group AVERAGE_READS_PER_GROUP
                        For each sample group require the average read count
                        to be at least this value in order to be significant
                        (default 10)
  --average-cpm-per-group AVERAGE_CPM_PER_GROUP
                        For each sample group require the average CPM to be at
                        least this value in order to be significant
  --min-cpm-per-group MIN_CPM_PER_GROUP
                        For each sample group require the min CPM to be at
                        least this value in order to be significant
  --asm-proportion-of-gene ASM_PROPORTION_OF_GENE
                        Require the ASM CPM to be at least this proportion of
                        the gene CPM in at least 1 sample in order to be
                        significant (default 0.05)
```

#### detect_differential_isoforms.py

[src/rmats_long/detect_differential_isoforms.py](src/rmats_long/detect_differential_isoforms.py) detects differential isoform usage. The samples need to be separated into `--group-1` and `--group-2` input files as comma separated lists.

The main output file has these columns:
* `gene_id`
* `isoform_id`: isoform ID
* `lr`: likelihood ratio statistic
* `df`: degrees of freedom
* `pvalue`
* `adj_pvalue`: adjusted p-value
* `{sample_name}_proportion`: proportion of this isoform among all isoforms in this gene (1 column per sample)
* `group_1_average_proportion`
* `group_2_average_proportion`
* `delta_isoform_proportion`: `group_1_average_proportion - group_2_average_proportion`

`differential_isoforms_filtered.tsv` contains only the rows meeting the significance cutoffs.

A summary of the number of isoforms and genes passing the default filters will be printed to stdout. The counts using different filters can be printed using [count_significant_isoforms.py](#count_significant_isoformspy).

```
rmats-long detect_differential_isoforms.py -h

usage: detect_differential_isoforms.py [-h] [--abundance ABUNDANCE]
                                       [--asm-counts-dir ASM_COUNTS_DIR]
                                       --out-dir OUT_DIR --group-1 GROUP_1
                                       --group-2 GROUP_2
                                       [--num-threads NUM_THREADS]
                                       [--adj-pvalue ADJ_PVALUE]
                                       [--use-unadjusted-pvalue]
                                       [--delta-proportion DELTA_PROPORTION]
                                       [--covar-tsv COVAR_TSV] [--use-drimseq]
                                       [--min-cpm-per-asm MIN_CPM_PER_ASM]
                                       --sample-read-total-tsv
                                       SAMPLE_READ_TOTAL_TSV
                                       [--min-isoform-reads MIN_ISOFORM_READS]
                                       [--limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS]
                                       [--average-reads-per-group AVERAGE_READS_PER_GROUP]
                                       [--average-cpm-per-group AVERAGE_CPM_PER_GROUP]
                                       [--min-cpm-per-group MIN_CPM_PER_GROUP]
                                       [--asm-proportion-of-gene ASM_PROPORTION_OF_GENE]
                                       --gene-cpm-tsv GENE_CPM_TSV

Detect differential isoform expression using a multinomial model

options:
  -h, --help            show this help message and exit
  --abundance ABUNDANCE
                        The path to the abundance.esp file from ESPRESSO
  --asm-counts-dir ASM_COUNTS_DIR
                        The output directory from count_reads_for_asms.py
  --out-dir OUT_DIR     The path to use as the output directory
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list. The sample names
                        should match with the ESPRESSO abundance column names.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2.
  --num-threads NUM_THREADS
                        The number of threads to use (default: 1)
  --adj-pvalue ADJ_PVALUE
                        The cutoff for adjusted p-value (default: 0.05)
  --use-unadjusted-pvalue
                        Use pvalue instead of adj_pvalue for the cutoff
  --delta-proportion DELTA_PROPORTION
                        The cutoff for delta isoform proportion (default:
                        0.05)
  --covar-tsv COVAR_TSV
                        A .tsv with 1 line per sample. The first line has the
                        column names. The first column is sample_id. Each
                        additional column is a covariate.
  --use-drimseq         Use DRIMSeq instead of run_stat_model.py
  --min-cpm-per-asm MIN_CPM_PER_ASM
                        Only consider ASMs where at least 1 sample has at
                        least this CPM of reads assigned to the ASM. (default:
                        0)
  --sample-read-total-tsv SAMPLE_READ_TOTAL_TSV
                        A .tsv file with two columns: sample and total. The
                        1st line is the header
  --min-isoform-reads MIN_ISOFORM_READS
                        Only consider isoforms with at least this many reads
                        (default: 1)
  --limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS
                        Only consider the top N isoforms with the highest
                        total proportion across samples for each ASM (default:
                        50)
  --average-reads-per-group AVERAGE_READS_PER_GROUP
                        For each sample group require the average read count
                        to be at least this value in order to be significant
                        (default: 10)
  --average-cpm-per-group AVERAGE_CPM_PER_GROUP
                        For each sample group require the average CPM to be at
                        least this value in order to be significant
  --min-cpm-per-group MIN_CPM_PER_GROUP
                        For each sample group require the min CPM to be at
                        least this value in order to be significant
  --asm-proportion-of-gene ASM_PROPORTION_OF_GENE
                        Require the ASM CPM to be at least this proportion of
                        the gene CPM in at least 1 sample in order to be
                        significant (default: 0.05)
  --gene-cpm-tsv GENE_CPM_TSV
                        A .tsv file with gene_id as the 1st header and then an
                        additional header for each sample name
```

#### count_significant_isoforms.py

```
rmats-long count_significant_isoforms.py -h

usage: count_significant_isoforms.py [-h] --diff-transcripts DIFF_TRANSCRIPTS
                                     --diff-asms DIFF_ASMS --out-tsv OUT_TSV
                                     [--adj-pvalue ADJ_PVALUE]
                                     [--use-unadjusted-pvalue]
                                     [--delta-proportion DELTA_PROPORTION]
                                     [--is-asm]
                                     [--average-reads-per-group AVERAGE_READS_PER_GROUP]
                                     [--average-cpm-per-group AVERAGE_CPM_PER_GROUP]
                                     [--min-cpm-per-group MIN_CPM_PER_GROUP]
                                     [--asm-proportion-of-gene ASM_PROPORTION_OF_GENE]
                                     --group-1 GROUP_1 --group-2 GROUP_2

Count isoforms that meet the cutoff values

options:
  -h, --help            show this help message and exit
  --diff-transcripts DIFF_TRANSCRIPTS
                        The path to the differential transcript results
  --diff-asms DIFF_ASMS
                        The path to the differential asms results
  --out-tsv OUT_TSV     The path to write transcripts that meet the cutoff
                        values
  --adj-pvalue ADJ_PVALUE
                        The cutoff for adjusted p-value (default: 0.05)
  --use-unadjusted-pvalue
                        Use pvalue instead of adj_pvalue for the cutoff
  --delta-proportion DELTA_PROPORTION
                        The cutoff for delta isoform proportion (default:
                        0.05)
  --is-asm              Use if running with ASM output
  --average-reads-per-group AVERAGE_READS_PER_GROUP
                        For each sample group require the average read count
                        to be at least this value in order to be significant
                        (default: 10)
  --average-cpm-per-group AVERAGE_CPM_PER_GROUP
                        For each sample group require the average CPM to be at
                        least this value in order to be significant
  --min-cpm-per-group MIN_CPM_PER_GROUP
                        For each sample group require the min CPM to be at
                        least this value in order to be significant
  --asm-proportion-of-gene ASM_PROPORTION_OF_GENE
                        Require the ASM CPM to be at least this proportion of
                        the gene CPM in at least 1 sample in order to be
                        significant (default: 0.05)
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2
```

#### visualize_isoforms.py

[src/rmats_long/visualize_isoforms.py](src/rmats_long/visualize_isoforms.py) creates plots showing the isoform abundance and structure. The `--gene-id` can be selected from the [differential isoform test](#detect_differential_isoformspy). The `--abundance` file could come from the ESPRESSO output or be created by `rmats_long.py`. By default, the most abundant isoforms for the gene will be plotted. Specific isoforms can be plotted with `--main-transcript-ids` or isoforms can be determined automatically if `--diff-transcripts` or `--gencode-gtf` are given. The most significant isoform and another significant isoform with opposite `delta_isoform_proportion` will be chosen from `--diff-transcripts` and the Ensembl canonical transcript will be selected based on a tag in the `--gencode-gtf`.

```
rmats-long visualize_isoforms.py -h

usage: visualize_isoforms.py [-h] --gene-id GENE_ID [--gene-name GENE_NAME]
                             [--asm-id ASM_ID] [--abundance ABUNDANCE]
                             --updated-gtf UPDATED_GTF
                             [--gencode-gtf GENCODE_GTF] [--asm-gtf ASM_GTF]
                             [--diff-transcripts DIFF_TRANSCRIPTS] --out-dir
                             OUT_DIR [--plot-file-type {.pdf,.png,all}]
                             [--main-transcript-ids MAIN_TRANSCRIPT_IDS]
                             [--max-transcripts MAX_TRANSCRIPTS]
                             [--intron-scaling INTRON_SCALING]
                             [--group-1 GROUP_1] [--group-2 GROUP_2]
                             [--group-1-name GROUP_1_NAME]
                             [--group-2-name GROUP_2_NAME] [--is-asm]
                             [--start-coord START_COORD]
                             [--end-coord END_COORD] [--graph-file GRAPH_FILE]
                             [--out-transcript-colors OUT_TRANSCRIPT_COLORS]

Visualize the structure and abundance of isoforms

options:
  -h, --help            show this help message and exit
  --gene-id GENE_ID     The gene_id to visualize
  --gene-name GENE_NAME
                        The name for the gene (used as plot title). If not
                        given then the gene_name from --gencode-gtf will be
                        used. If no other name is found then --gene-id is used
                        as a default
  --asm-id ASM_ID       The asm_id to use with --asm-gtf
  --abundance ABUNDANCE
                        The path to the abundance.esp file from ESPRESSO
  --updated-gtf UPDATED_GTF
                        The path to the updated.gtf file from ESPRESSO
  --gencode-gtf GENCODE_GTF
                        The path to a gencode annotation.gtf file. Can be used
                        to identify the gene_name and Ensembl canonical
                        isoform
  --asm-gtf ASM_GTF     The path to a .gtf file with ASM transcripts. Can be
                        used to select which gene transcripts to plot
  --diff-transcripts DIFF_TRANSCRIPTS
                        The path to the differential transcript results. Can
                        be used to determine --main-transcript-ids
  --out-dir OUT_DIR     The path to use as the output directory
  --plot-file-type {.pdf,.png,all}
                        The file type for output plots (default .png))
  --main-transcript-ids MAIN_TRANSCRIPT_IDS
                        A comma separated list of transcript IDs to plot as
                        the main transcripts. If not given then the most
                        significant isoform from --diff-transcripts, a second
                        significant isoform with a delta proportion in the
                        opposite direction, and the Ensembl canonical isoform
                        from --gencode-gtf will be used if possible
  --max-transcripts MAX_TRANSCRIPTS
                        How many transcripts to plot individually. The
                        remaining transcripts in the gene will be grouped
                        together (default 7)
  --intron-scaling INTRON_SCALING
                        The factor to use to reduce intron length in the plot.
                        A value of 2 would reduce introns to 1/2 of the
                        original plot length (default 1)
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list. The sample names
                        should match with the ESPRESSO abundance column names.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2.
  --group-1-name GROUP_1_NAME
                        A name for group 1 (default group 1)
  --group-2-name GROUP_2_NAME
                        A name for group 2 (default group 2)
  --is-asm              Use ASM data
  --start-coord START_COORD
                        Indicate the start of a region at this coordinate
  --end-coord END_COORD
                        Indicate the end of a region at this coordinate
  --graph-file GRAPH_FILE
                        The path to graph_{chr}.txt which has the splice graph
                        details for this gene
  --out-transcript-colors OUT_TRANSCRIPT_COLORS
                        Where to write a .tsv of colors used for transcripts
```

#### classify_isoform_differences.py

[src/rmats_long/classify_isoform_differences.py](src/rmats_long/classify_isoform_differences.py) compares the structures of isoforms within a gene by calling [FindAltTSEvents.py](#findalttseventspy) with a "main" isoform and a second isoform (or all other isoforms in the gene by default).

```
rmats-long classify_isoform_differences.py -h

usage: classify_isoform_differences.py [-h] --main-transcript-id
                                       MAIN_TRANSCRIPT_ID --updated-gtf
                                       UPDATED_GTF [--gencode-gtf GENCODE_GTF]
                                       --out-tsv OUT_TSV
                                       [--second-transcript-id SECOND_TRANSCRIPT_ID]

Compare the structures of isoforms within a gene

options:
  -h, --help            show this help message and exit
  --main-transcript-id MAIN_TRANSCRIPT_ID
                        The transcript_id of the main isoform in the .gtf file
  --updated-gtf UPDATED_GTF
                        The path to the updated.gtf file from ESPRESSO
  --gencode-gtf GENCODE_GTF
                        The path to a gencode annotation.gtf file. Can be used
                        to compare against isoforms not detected by ESPRESSO
  --out-tsv OUT_TSV     The path of the output file
  --second-transcript-id SECOND_TRANSCRIPT_ID
                        If given, only compare the main transcript to this
                        transcript
```

#### FindAltTSEvents.py

[src/rmats_long/FindAltTSEvents.py](src/rmats_long/FindAltTSEvents.py) compares the structures of any two transcript isoforms. Local differences in transcript structure are classified into 7 basic alternative splicing categories:

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
rmats-long FindAltTSEvents.py -h

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

#### organize_gene_info_by_chr.py

[src/rmats_long/organize_gene_info_by_chr.py](src/rmats_long/organize_gene_info_by_chr.py) creates a directory of sorted files which allows other steps to be more efficient. `chr_name_id_mapping.tsv` identifies which chromosome corresponds to each file.

```
rmats-long organize_gene_info_by_chr.py -h

usage: organize_gene_info_by_chr.py [-h] --gtf GTF --out-dir OUT_DIR

Create 1 file per chr with transcript info by gene

options:
  -h, --help         show this help message and exit
  --gtf GTF          The path to a gtf with transcript info
  --out-dir OUT_DIR  The directory to create and where chr files will be
                     written
```

#### simplify_alignment_info.py

```
rmats-long simplify_alignment_info.py -h

usage: simplify_alignment_info.py [-h] --in-file IN_FILE --out-tsv OUT_TSV
                                  [--sort-buffer-size SORT_BUFFER_SIZE]

Read a SAM or BAM file and output the splice junctions for each read

options:
  -h, --help            show this help message and exit
  --in-file IN_FILE     The path to a .sam or .bam file
  --out-tsv OUT_TSV     The path to write the output reads
  --sort-buffer-size SORT_BUFFER_SIZE
                        Used for the --buffer-size argument of sort. Default:
                        2G
```

#### organize_alignment_info_by_gene_and_chr.py

Similar to [organize_gene_info_by_chr.py](#organize_gene_info_by_chrpy) [src/rmats_long/organize_alignment_info_by_gene_and_chr.py](src/rmats_long/organize_alignment_info_by_gene_and_chr.py) creates a directory of sorted files which allows other steps to be more efficient. `sample_read_totals.tsv` records the total number of reads for each sample.

```
rmats-long organize_alignment_info_by_gene_and_chr.py -h

usage: organize_alignment_info_by_gene_and_chr.py [-h] --gtf-dir GTF_DIR
                                                  --out-dir OUT_DIR
                                                  --samples-tsv SAMPLES_TSV
                                                  [--sort-buffer-size SORT_BUFFER_SIZE]
                                                  [--num-threads NUM_THREADS]

Create 1 file per chr with read info by gene

options:
  -h, --help            show this help message and exit
  --gtf-dir GTF_DIR     The output directory from organize_gene_info_by_chr.py
  --out-dir OUT_DIR     The directory to create and where new files will be
                        written
  --samples-tsv SAMPLES_TSV
                        The path to a file where each line has 2 tab separated
                        columns: sample name, then a path to an output file
                        from simplify_alignment_info.py. A sample name can
                        have multiple lines if it has multiple input files.
  --sort-buffer-size SORT_BUFFER_SIZE
                        Used for the --buffer-size argument of sort. Default:
                        2G
  --num-threads NUM_THREADS
                        The number of threads to use (default 1)
```

#### detect_splicing_events.py

[src/rmats_long/detect_splicing_events.py](src/rmats_long/detect_splicing_events.py) finds ASMs using a `.gtf` file. The transcripts for each gene are used to build a splice graph which is filtered if `--align-dir` and `--min-reads-per-edge` are provided. Then a search is performed for regions in the graph with multiple possible paths. In order to avoid computational issues with large and complex genes, the search is limited using `--max-nodes-in-event`.

```
rmats-long detect_splicing_events.py -h

usage: detect_splicing_events.py [-h] --gtf-dir GTF_DIR
                                 [--align-dir ALIGN_DIR] --out-dir OUT_DIR
                                 [--max-nodes-in-event MAX_NODES_IN_EVENT]
                                 [--max-paths-in-event MAX_PATHS_IN_EVENT]
                                 [--num-threads NUM_THREADS]
                                 [--min-reads-per-edge MIN_READS_PER_EDGE]
                                 [--output-full-gene-asm]
                                 [--output-basic-events]
                                 [--simplify-gene-isoform-endpoints]
                                 [--filter-gene-isoforms-by-edge]
                                 [--output-strict-only] [--novel-junctions]

Detect alternative splicing events from a set of isoforms

options:
  -h, --help            show this help message and exit
  --gtf-dir GTF_DIR     The output directory from organize_gene_info_by_chr.py
  --align-dir ALIGN_DIR
                        The output directory from
                        organize_alignment_info_by_gene_and_chr.py
  --out-dir OUT_DIR     The directory to create and where chr files will be
                        written
  --max-nodes-in-event MAX_NODES_IN_EVENT
                        Only look for events between nodes (splice sites) in
                        the splice graph that are at most --max-nodes-in-event
                        apart. Default: 50
  --max-paths-in-event MAX_PATHS_IN_EVENT
                        Only output events with at most --max-paths-in-event.
                        Default: 50
  --num-threads NUM_THREADS
                        how many threads to use. Default: 1
  --min-reads-per-edge MIN_READS_PER_EDGE
                        Only include an edge in the splice graph if there are
                        at least this many supporting reads. Default: 5
  --output-full-gene-asm
                        Output only ASMs that cover the entire gene (from
                        transcript start to end)
  --output-basic-events
                        Output only ASMs for basic events (SE, A5SS, A3SS,
                        MXE, RI)
  --simplify-gene-isoform-endpoints
                        Combine gene isoforms where the only difference is the
                        transcripts start and/or end
  --filter-gene-isoforms-by-edge
                        With --output-full-gene-asm, require each isoform to
                        have --min-reads-per-edge for each splice junction
  --output-strict-only  Only output events where is_strict is True
  --novel-junctions     Consider splice junctions found in --align-dir but not
                        --gtf-dir
```

#### count_reads_for_asms.py

```
rmats-long count_reads_for_asms.py -h

usage: count_reads_for_asms.py [-h] --event-dir EVENT_DIR --align-dir
                               ALIGN_DIR --gtf-dir GTF_DIR --out-dir OUT_DIR
                               [--num-threads NUM_THREADS]
                               [--sort-buffer-size SORT_BUFFER_SIZE]

Determine read counts for isoforms in ASMs

options:
  -h, --help            show this help message and exit
  --event-dir EVENT_DIR
                        The output directory from detect_splicing_events.py
  --align-dir ALIGN_DIR
                        The output directory from
                        organize_alignment_info_by_gene_and_chr.py
  --gtf-dir GTF_DIR     The output directory from organize_gene_info_by_chr.py
  --out-dir OUT_DIR     The directory to write ASM read counts in output files
                        by chr
  --num-threads NUM_THREADS
                        how many threads to use
  --sort-buffer-size SORT_BUFFER_SIZE
                        Used for the --buffer-size argument of sort. Default:
                        2G
```

#### plot_splice_graph.py

[src/rmats_long/plot_splice_graph.py](src/rmats_long/plot_splice_graph.py) creates a diagram of the splice graph for an ASM. The edges in the graphs will show the number of reads aligned to that junction or exon if [detect_splicing_events.py](#detect_splicing_eventspy) was run with `--align-dir`.

```
rmats-long plot_splice_graph.py -h

usage: plot_splice_graph.py [-h] --event-dir EVENT_DIR [--chr CHR]
                            [--gene-id GENE_ID] [--asm-id ASM_ID] --out-file
                            OUT_FILE [--edge-label-lines]
                            [--min-edge-weight MIN_EDGE_WEIGHT]
                            [--show-extra-nodes]

Plot the splice graph for a gene or ASM

options:
  -h, --help            show this help message and exit
  --event-dir EVENT_DIR
                        The output directory from detect_splicing_events.py
  --chr CHR             The chromosome name with the gene or ASM
  --gene-id GENE_ID     The gene ID to plot
  --asm-id ASM_ID       The ASM ID to plot
  --out-file OUT_FILE   The file for the output plot
  --edge-label-lines    draw a line connecting each edge label to its edge
  --min-edge-weight MIN_EDGE_WEIGHT
                        only draw edges with at least this weight
  --show-extra-nodes    include nodes other than the main annotated
                        coordinates
```

#### plot_simple_splice_graph.py

```
rmats-long plot_simple_splice_graph.py -h

usage: plot_simple_splice_graph.py [-h] --gtf GTF --out-file OUT_FILE
                                   [--gene-id GENE_ID] [--asm-id ASM_ID]
                                   [--transcript-ids TRANSCRIPT_IDS]
                                   [--intron-scaling INTRON_SCALING]
                                   [--equal-spacing]
                                   [--junction-counts JUNCTION_COUNTS]
                                   [--color-by-isoform]
                                   [--color-first-transcript-only]
                                   [--transcript-colors TRANSCRIPT_COLORS]
                                   [--ids-from-color-file]
                                   [--show-isoform-endpoint-symbols]
                                   [--font-size FONT_SIZE]
                                   [--plot-width PLOT_WIDTH]
                                   [--plot-height PLOT_HEIGHT]
                                   [--plot-dpi PLOT_DPI]

Output a diagram of a splicing graph

options:
  -h, --help            show this help message and exit
  --gtf GTF             The path to a .gtf file with isoform definitions
  --out-file OUT_FILE   The file for the output plot
  --gene-id GENE_ID     Only plot isoforms with this "gene_id" attribute
  --asm-id ASM_ID       Only plot isoforms with this "asm_id" attribute
  --transcript-ids TRANSCRIPT_IDS
                        A comma separated list of isoform IDs. Only plot
                        isoforms with one of these "transcript_id" attributes.
  --intron-scaling INTRON_SCALING
                        The factor to use to reduce intron length in the plot.
                        A value of 2 would reduce introns to 1/2 of the
                        original plot length. (default 1)
  --equal-spacing       Plot each splice site a fixed distance from the
                        previous splice site
  --junction-counts JUNCTION_COUNTS
                        A tab separated file with headers: ["chr", "start",
                        "end", "count"]. Any junction plotted between "start"
                        and "end" will have "count" displayed.
  --color-by-isoform    Use a color for each isoform and apply that color to
                        exons and junctions only used by that isoform
  --color-first-transcript-only
                        Color all exons and junctions for the 1st value from
                        --transcript-ids
  --transcript-colors TRANSCRIPT_COLORS
                        A tab separated file with headers: ["transcript_id",
                        "color"].The color column is an RGB hex string
                        (#RRGGBB)
  --ids-from-color-file
                        set --transcripts-ids from the transcript_id column of
                        --transcript-colors
  --show-isoform-endpoint-symbols
                        Add a symbol at each isoform's start and a symbol at
                        each isoform's end
  --font-size FONT_SIZE
                        The font size for the plot (default 8)
  --plot-width PLOT_WIDTH
                        The plot width in inches (default 12)
  --plot-height PLOT_HEIGHT
                        The plot height in inches (default 6)
  --plot-dpi PLOT_DPI   The plot resolution in dots per inch (default 300)
```

#### create_gtf_from_asm_definitions.py

```
rmats-long create_gtf_from_asm_definitions.py -h

usage: create_gtf_from_asm_definitions.py [-h] --event-dir EVENT_DIR --out-gtf
                                          OUT_GTF

Create a .gtf file with ASM definitions

options:
  -h, --help            show this help message and exit
  --event-dir EVENT_DIR
                        The output directory from detect_splicing_events.py
  --out-gtf OUT_GTF     The output .gtf file to create
```

#### run_stat_model.py

```
rmats-long run_stat_model.py -h

usage: run_stat_model.py [-h] [--counts-dir COUNTS_DIR]
                         [--abundance ABUNDANCE] --out-dir OUT_DIR --group-1
                         GROUP_1 --group-2 GROUP_2 [--num-threads NUM_THREADS]
                         [--group-1-name GROUP_1_NAME]
                         [--group-2-name GROUP_2_NAME]
                         [--min-isoform-reads MIN_ISOFORM_READS]
                         [--limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS]
                         [--min-cpm-per-asm MIN_CPM_PER_ASM]
                         --sample-read-total-tsv SAMPLE_READ_TOTAL_TSV
                         [--em-tolerance EM_TOLERANCE]
                         [--em-max-iter EM_MAX_ITER]
                         [--random-seed RANDOM_SEED]
                         [--progress-every-n PROGRESS_EVERY_N]
                         [--sort-buffer-size SORT_BUFFER_SIZE]
                         [--covar-tsv COVAR_TSV]

Run a statistical model on isoform counts

options:
  -h, --help            show this help message and exit

Read compatibility:
  Use either Read compatibility or Transcript abundance

  --counts-dir COUNTS_DIR
                        The input directory with read isoform compatibility
                        files

Transcript abundance:
  Use either Transcript abundance or Read compatibility

  --abundance ABUNDANCE
                        The path to the transcript abundance by sample

Required:
  --out-dir OUT_DIR     The directory to write output to
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2

Optional:
  --num-threads NUM_THREADS
                        How many threads to use (default 1)
  --group-1-name GROUP_1_NAME
                        A name for group 1 (default group_1)
  --group-2-name GROUP_2_NAME
                        A name for group 2 (default group_2)
  --min-isoform-reads MIN_ISOFORM_READS
                        Only consider isoforms with at least this many reads
                        (default 1)
  --limit-asm-to-top-n-isoforms LIMIT_ASM_TO_TOP_N_ISOFORMS
                        Only consider the top N isoforms with the highest
                        total proportion across samples for each ASM (default:
                        50)
  --min-cpm-per-asm MIN_CPM_PER_ASM
                        Only consider ASMs where at least 1 sample has at
                        least this CPM of reads assigned to the ASM. (default
                        0)
  --sample-read-total-tsv SAMPLE_READ_TOTAL_TSV
                        A .tsv file with two columns: sample and total. The
                        1st line is the header
  --em-tolerance EM_TOLERANCE
                        Stop performing EM iterations when the change is at or
                        below this value (default 0.001)
  --em-max-iter EM_MAX_ITER
                        Perform at most this many EM iterations (default 100)
  --random-seed RANDOM_SEED
                        Passed to R base::set.seed() (default 123)
  --progress-every-n PROGRESS_EVERY_N
                        Print a status message after a certain number of ASMs
                        (default 100)
  --sort-buffer-size SORT_BUFFER_SIZE
                        Used for the --buffer-size argument of sort. Default:
                        2G
  --covar-tsv COVAR_TSV
                        A .tsv with 1 line per sample. The first line has the
                        column names. The first column is sample_id. Each
                        additional column is a covariate.
```

#### check_asm_coupling.py

```
rmats-long check_asm_coupling.py -h

usage: check_asm_coupling.py [-h] --event-dir EVENT_DIR --asm-counts-dir
                             ASM_COUNTS_DIR --count-tsv COUNT_TSV
                             [--group-1 GROUP_1] [--group-2 GROUP_2]
                             [--group-1-name GROUP_1_NAME]
                             [--group-2-name GROUP_2_NAME]
                             [--group-tsv GROUP_TSV] --out-tsv OUT_TSV
                             --isoform-out-tsv ISOFORM_OUT_TSV

Check pairs of ASMs for splicing coupling

options:
  -h, --help            show this help message and exit
  --event-dir EVENT_DIR
                        The output directory from detect_splicing_events.py
  --asm-counts-dir ASM_COUNTS_DIR
                        The output directory from count_reads_for_asms.py
  --count-tsv COUNT_TSV
                        The count.tsv from rmats_long.py
  --group-1 GROUP_1     The path to a file listing the sample names for group
                        1. The file should have a single line with the sample
                        names as a comma separated list.
  --group-2 GROUP_2     The path to a file listing the sample names for group
                        2
  --group-1-name GROUP_1_NAME
                        A name for --group-1 (default group 1)
  --group-2-name GROUP_2_NAME
                        A name for --group-2 (default group 2)
  --group-tsv GROUP_TSV
                        A .tsv with the headers "sample" and "group". Can be
                        used instead of --group-1 and --group-2 to define more
                        than two groups
  --out-tsv OUT_TSV     Where to write coupling results
  --isoform-out-tsv ISOFORM_OUT_TSV
                        Where to write the main isoform determined for each
                        ASM
```
