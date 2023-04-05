library('BiocParallel')
library('DRIMSeq')
library('ggplot2')

args <- base::commandArgs(trailingOnly=TRUE)
espresso_abundance_path <- args[1]
out_dir_path <- args[2]
num_threads <- base::as.integer(args[3])
group_1_file_path <- args[4]
group_2_file_path <- args[5]

string_is_all_whitespace <- function(string) {
    return(base::grepl('^\\s*$', string)[1])
}

read_sample_file <- function(file_path) {
    connection <- base::file(file_path)
    base::open(connection)
    lines <- base::readLines(connection)
    base::close(connection)
    samples_string <- lines[1]
    num_lines <- base::length(lines)
    if (num_lines > 1) {
        for (line_i in 2:num_lines) {
            if (!string_is_all_whitespace(lines[line_i])) {
                stop(base::paste0('unexpected extra line in ', file_path))
            }
        }
    }
    sample_names <- base::strsplit(samples_string, ',')[[1]]
    return(sample_names)
}

make_sample_group_data_frame <- function(sample_names, group_1_file_path,
                                         group_2_file_path) {
    group_1_sample_names <- read_sample_file(group_1_file_path)
    group_2_sample_names <- read_sample_file(group_2_file_path)
    num_samples <- base::length(sample_names)
    num_names_with_group <- (base::length(group_1_sample_names) +
                             base::length(group_2_sample_names))
    if (num_samples != num_names_with_group) {
        stop(base::paste0('number of samples in group files (',
                          num_names_with_group,
                          ') should match number in abundance file (',
                          num_samples, ')'))
    }
    group_names <- base::vector(mode='character',
                                length=num_samples)
    for (sample_i in 1:base::length(sample_names)) {
        sample_name <- sample_names[sample_i]
        if (sample_name %in% group_1_sample_names) {
            group_names[sample_i] <- 'group_1'
        } else if (sample_name %in% group_2_sample_names) {
            group_names[sample_i] <- 'group_2'
        } else {
            stop(base::paste0(sample_name, ' not matched to group 1 or 2'))
        }
    }
    df <- base::data.frame(sample_id=sample_names,
                           group=group_names)
    return(df)
}

read_tsv_with_rows_as_first_column <- function(path, first_column_name) {
    ## check.names=FALSE prevents '-' (dash) from being replaced with '.' in
    ## the sample names.
    data <- utils::read.table(path, header=TRUE, sep='\t', quote='',
                              row.names=1, check.names=FALSE)
    num_columns <- base::ncol(data)
    ## read.table sometimes keeps the first column name as a regular column
    ## instead of as the header for the row names.
    ## In that case, remove the final column of NA values and shift the
    ## column names by 1.
    if (base::colnames(data)[1] == first_column_name) {
        num_columns <- base::ncol(data)
        base::colnames(data) <- base::colnames(data)[2:num_columns]
        data[num_columns] <- NULL
    }
    return(data)
}

read_abundance <- function(espresso_abundance_path) {
    ## columns: transcript_ID, transcript_name, gene_ID, sample_name1, sample_name2, ...
    data <- read_tsv_with_rows_as_first_column(espresso_abundance_path,
                                               'transcript_ID')
    all_column_names <- base::colnames(data)
    sample_names <- all_column_names[3:base::length(all_column_names)]
    return(base::list(data=data, sample_names=sample_names))
}

make_counts_data_frame <- function(abundance) {
    ## rows are transcripts
    ## gene_id column (filter out NA genes)
    ## feature_id column with transcript IDs
    ## column for each sample with counts
    not_na_gene_bools <- !base::is.na(abundance$data$gene_ID)
    gene_ids <- abundance$data$gene_ID[not_na_gene_bools]
    transcript_ids <- base::rownames(abundance$data)[not_na_gene_bools]
    df <- base::data.frame(gene_id=gene_ids, feature_id=transcript_ids)
    for (sample in abundance$sample_names) {
        df[[sample]] <- abundance$data[[sample]][not_na_gene_bools]
    }
    return(df)
}

## From the vignette:
## [the default params give] a very relaxed filtering, where transcripts with
## zero expression in all the samples and genes with only one transcript
## remain[ing] are removed.
filter_data <- function(drim_data, sample_df) {
    return(DRIMSeq::dmFilter(drim_data))
}

plot_transcripts_per_gene <- function(drim_data, out_dir_path) {
    plot <- DRIMSeq::plotData(drim_data)
    plot_file_path <- base::paste0(out_dir_path, '/transcripts_per_gene.png')
    ggplot2::ggsave(plot=plot, plot_file_path)
}

plot_precision <- function(drim_data, out_dir_path) {
    plot <- DRIMSeq::plotPrecision(drim_data)
    plot_file_path <- base::paste0(out_dir_path,
                                   '/precision_by_gene_expression.png')
    ggplot2::ggsave(plot=plot, plot_file_path)
}

plot_pvalues <- function(drim_data, out_dir_path) {
    plot <- DRIMSeq::plotPValues(drim_data)
    plot_file_path <- base::paste0(out_dir_path, '/gene_pvalues.png')
    ggplot2::ggsave(plot=plot, plot_file_path)

    plot <- DRIMSeq::plotPValues(drim_data, level='feature')
    plot_file_path <- base::paste0(out_dir_path, '/transcript_pvalues.png')
    ggplot2::ggsave(plot=plot, plot_file_path)
}

write_output <- function(drim_data, out_dir_path) {
    gene_results <- DRIMSeq::results(drim_data)
    out_path <- base::paste0(out_dir_path, '/differential_genes.tsv')
    utils::write.table(gene_results, file=out_path, quote=FALSE, sep='\t',
                       row.names=FALSE, col.names=TRUE)

    transcript_results <- DRIMSeq::results(drim_data, level='feature')
    out_path <- base::paste0(out_dir_path, '/differential_transcripts.tsv')
    utils::write.table(transcript_results, file=out_path, quote=FALSE, sep='\t',
                       row.names=FALSE, col.names=TRUE)
}

create_out_dir <- function(dir_path) {
    if (!base::file.exists(dir_path)) {
        base::cat(base::paste0('creating: ', dir_path, '\n'))
        base::dir.create(dir_path)
    }
}

main <- function() {
    create_out_dir(out_dir_path)
    bpparam <- BiocParallel::MulticoreParam(workers=num_threads)
    base::set.seed(123)  ## DRIMSeq recommends setting the seed

    ## Read and filter input data
    abundance <- read_abundance(espresso_abundance_path)
    sample_df <- make_sample_group_data_frame(abundance$sample_names,
                                              group_1_file_path,
                                              group_2_file_path)
    counts_df <- make_counts_data_frame(abundance)
    drim_data <- DRIMSeq::dmDSdata(counts_df, sample_df)
    drim_data <- filter_data(drim_data, sample_df)
    plot_transcripts_per_gene(drim_data, out_dir_path)

    ## Calculate results
    design_full <- stats::model.matrix(~ group,
                                       data=DRIMSeq::samples(drim_data))
    null_model <- stats::model.matrix(~ 1, data=DRIMSeq::samples(drim_data))
    base::cat('\nAbout to run DRIMSeq::dmPrecision(...) which could take a while...\n')
    drim_data <- DRIMSeq::dmPrecision(drim_data, design=design_full,
                                      BPPARAM=bpparam)
    plot_precision(drim_data, out_dir_path)

    drim_data <- DRIMSeq::dmFit(drim_data, design=design_full,
                                BPPARAM=bpparam)
    drim_data <- DRIMSeq::dmTest(drim_data, design=null_model,
                                 BPPARAM=bpparam)

    ## Output final results
    plot_pvalues(drim_data, out_dir_path)
    write_output(drim_data, out_dir_path)
}

main()
