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

split_multi_gene_rows <- function(df) {
    splits <- base::strsplit(df$gene_id, ',')
    split_lengths <- base::sapply(splits, function(split) {base::length(split)})
    highest_gene_count <- base::max(split_lengths)
    if (highest_gene_count == 1) {
        renamed <- base::data.frame(orig=base::c(), renamed=base::c())
        return(base::list(counts_df=df, renamed_transcripts=renamed))
    }

    multi_gene_bools <- split_lengths > 1
    multi_gene_rows <- df[multi_gene_bools, ]
    multi_gene_splits <- splits[multi_gene_bools]
    without_multi_gene_rows <- df[!multi_gene_bools, ]
    all_colnames <- base::colnames(multi_gene_rows)
    non_gene_colnames <- all_colnames[all_colnames != 'gene_id']
    non_id_colnames <- non_gene_colnames[non_gene_colnames != 'feature_id']
    orig_transcript_names <- base::c()
    new_transcript_names <- base::c()
    new_single_gene_rows <- NULL
    for (multi_row_i in 1:base::nrow(multi_gene_rows)) {
        split <- multi_gene_splits[[multi_row_i]]
        orig_row <- multi_gene_rows[multi_row_i, ]
        rename_i <- 0
        for (gene_id in split) {
            new_transcript_name <- base::paste0(orig_row$feature_id, '_', rename_i)
            orig_transcript_names <- base::append(orig_transcript_names, orig_row$feature_id)
            new_transcript_names <- base::append(new_transcript_names, new_transcript_name)
            rename_i <- rename_i + 1
            new_df <- base::data.frame(gene_id=gene_id, feature_id=new_transcript_name)
            for (colname in non_id_colnames) {
                new_df[[colname]] <- orig_row[[colname]]
            }
            if (base::is.null(new_single_gene_rows)) {
                new_single_gene_rows <- new_df
            } else {
                new_single_gene_rows <- base::rbind(new_single_gene_rows, new_df)
            }
        }
    }

    result_df <- base::rbind(without_multi_gene_rows, new_single_gene_rows)
    renamed_transcripts <- base::data.frame(orig=orig_transcript_names, renamed=new_transcript_names)
    return(base::list(counts_df=result_df, renamed_transcripts=renamed_transcripts))
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

    return(split_multi_gene_rows(df))
}

## DRIMSeq behavior depends on the order of transcript IDs within a gene.
## If the last transcript ID for a gene has no counts for either sample group,
## then DRIMSeq will give NA results for that gene.
## Sort so that transcript IDs that have counts for both sample groups are last.
sort_counts_df <- function(counts_df, sample_df) {
    groups <- base::unique(sample_df$group)
    groups <- base::sort(groups)
    sort_columns <- list(gene_id=counts_df$gene_id)
    for (group in groups) {
        samples <- sample_df$sample_id[sample_df$group == group]
        group_sums <- counts_df[, samples]
        if (base::length(samples) > 1) {
            group_sums <- base::rowSums(group_sums)
        }
        sort_columns[[group]] <- group_sums > 0
    }
    sort_columns$feature_id <- counts_df$feature_id
    ## sort by: gene_id, has_counts_group_1, has_counts_group_2, transcript_ID
    sort_order <- base::order(sort_columns[[1]], sort_columns[[2]],
                              sort_columns[[3]], sort_columns[[4]])
    sorted_counts_df <- counts_df[sort_order, ]
    return(sorted_counts_df)
}

## From the vignette:
## [the default params give] a very relaxed filtering, where transcripts with
## zero expression in all the samples and genes with only one transcript
## remain[ing] are removed.
##
## Then from the documentation for dmFilter:
## 'min_samps_feature_expr' and 'min_samps_feature_prop'
## defines the minimal number of samples where features are required
## to be expressed at the minimal levels of counts 'min_feature_expr'
## or proportions 'min_feature_prop'. In differential transcript/exon
## usage analysis, we suggest using 'min_samps_feature_expr' and
## 'min_samps_feature_prop' equal to the minimal number of replicates
## in any of the conditions.
filter_data <- function(drim_data, sample_df) {
    groups <- base::unique(sample_df$group)
    samples_per_group <- base::sapply(groups, function(group) {
        return(base::sum(sample_df$group == group))
    })
    min_samples_per_group <- base::min(samples_per_group)
    return(DRIMSeq::dmFilter(drim_data,
                             min_samps_feature_expr=min_samples_per_group,
                             min_feature_expr=1))
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

write_output <- function(drim_data, renamed_transcripts, out_dir_path) {
    gene_results <- DRIMSeq::results(drim_data)
    out_path <- base::paste0(out_dir_path, '/differential_genes.tsv')
    utils::write.table(gene_results, file=out_path, quote=FALSE, sep='\t',
                       row.names=FALSE, col.names=TRUE)

    transcript_results <- DRIMSeq::results(drim_data, level='feature')
    fixed_feature_ids <- base::sapply(transcript_results$feature_id, function(feature_id) {
        match_i <- renamed_transcripts$renamed == feature_id
        match <- renamed_transcripts$orig[match_i]
        if (base::length(match) == 1) {
            return(match)
        }
        return(feature_id)
    })
    transcript_results$feature_id <- fixed_feature_ids
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
    counts_and_renamed_transcripts <- make_counts_data_frame(abundance)
    counts_df <- counts_and_renamed_transcripts$counts_df
    renamed_transcripts <- counts_and_renamed_transcripts$renamed_transcripts
    sorted_counts_df <- sort_counts_df(counts_df, sample_df)
    drim_data <- DRIMSeq::dmDSdata(sorted_counts_df, sample_df)
    drim_data <- filter_data(drim_data, sample_df)
    plot_transcripts_per_gene(drim_data, out_dir_path)

    ## Calculate results
    design_full <- stats::model.matrix(~ group,
                                       data=DRIMSeq::samples(drim_data))
    null_model <- stats::model.matrix(~ 1, data=DRIMSeq::samples(drim_data))
    base::cat('\nAbout to run DRIMSeq::dmPrecision(...) which could take a while...\n')
    ## From the DRIMSeq documentation:
    ## add_uniform: Whether to add a small fractional count to zeros, (adding
    ##           a uniform random variable between 0 and 0.1). This option
    ##           allows for the fitting of genewise precision and coefficients
    ##           for genes with two features having all zero for one group, or
    ##           the last feature having all zero for one group.
    ##
    ## prec_moderation: Precision moderation method. One can choose to shrink
    ##      the precision estimates toward the common precision
    ##      ("common") or toward the (precision versus mean expression)
    ##      trend ("trended")
    ##
    ## If there are only 2 samples then both 'trended' (the default) and 'common'
    ## result in an error. prec_moderation='none' avoids the error.
    if (base::nrow(sample_df) == 2) {
       prec_moderation <- 'none'
    } else {
        prec_moderation <- 'trended'
    }
    drim_data <- DRIMSeq::dmPrecision(drim_data, design=design_full,
                                      add_uniform=TRUE,
                                      prec_moderation=prec_moderation,
                                      BPPARAM=bpparam)
    plot_precision(drim_data, out_dir_path)

    drim_data <- DRIMSeq::dmFit(drim_data, design=design_full,
                                add_uniform=TRUE,
                                BPPARAM=bpparam)
    drim_data <- DRIMSeq::dmTest(drim_data, design=null_model,
                                 BPPARAM=bpparam)

    ## Output final results
    plot_pvalues(drim_data, out_dir_path)
    write_output(drim_data, renamed_transcripts, out_dir_path)
}

main()
