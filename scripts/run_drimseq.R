library('BiocParallel')
library('DRIMSeq')
library('ggplot2')

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
    group_sample_names <- base::c(group_1_sample_names, group_2_sample_names)
    num_names_with_group <- base::length(group_sample_names)
    extra_samples <- base::vector(mode='character')
    missing_samples <- base::c(group_1_sample_names, group_2_sample_names)
    sample_order <- base::vector(mode='character',
                                 length=num_names_with_group)
    group_names <- base::vector(mode='character',
                                length=num_names_with_group)
    next_i <- 1
    for (sample_name in sample_names) {
        missing_samples <- missing_samples[missing_samples != sample_name]
        group_name <- ''
        if (sample_name %in% group_1_sample_names) {
            group_name <- 'group_1'
        } else if (sample_name %in% group_2_sample_names) {
            group_name <- 'group_2'
        } else {
            extra_samples <- base::append(extra_samples, sample_name)
            next
        }

        sample_order[next_i] <- sample_name
        group_names[next_i] <- group_name
        next_i <- next_i + 1
    }

    df <- base::data.frame(sample_id=sample_order,
                           group=group_names)
    if (base::length(missing_samples) != 0) {
        comma_samples <- base::paste0(missing_samples, collapse=',')
        stop(base::paste0('Samples from group files not found in data: ', comma_samples))
    }
    if (base::length(extra_samples) != 0) {
        comma_samples <- base::paste0(extra_samples, collapse=',')
        base::cat('Skipping samples not assigned to a group: ', comma_samples, '\n')
    }

    return(df)
}

add_covar_columns <- function(covar_df, sample_df) {
    indices <- base::match(sample_df$sample_id, covar_df$sample_id)
    covars <- base::colnames(covar_df)
    covars <- covars[covars != 'sample_id']
    for (name in covars) {
        sample_df[[name]] <- covar_df[[name]][indices]
    }

    return(sample_df)
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

write_output_with_rename <- function(drim_data, renamed_transcripts, out_dir_path) {
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

run_drimseq <- function (counts_df, sample_df, renamed_transcripts,
                         num_threads, covar_df, out_dir_path) {
    create_out_dir(out_dir_path)
    bpparam <- BiocParallel::MulticoreParam(workers=num_threads)
    base::set.seed(123)  ## DRIMSeq recommends setting the seed

    sample_df <- add_covar_columns(covar_df, sample_df)
    sorted_counts_df <- sort_counts_df(counts_df, sample_df)
    drim_data <- DRIMSeq::dmDSdata(sorted_counts_df, sample_df)
    drim_data <- filter_data(drim_data, sample_df)
    plot_transcripts_per_gene(drim_data, out_dir_path)

    ## Calculate results
    full_formula <- ~ group
    null_formula <- ~ 1
    if (base::nrow(covar_df) > 0) {
        covars <- base::colnames(covar_df)
        covars <- covars[covars != 'sample_id']
        full_formula <- stats::reformulate(base::c('group', covars))
        null_formula <- stats::reformulate(covars)
    }
    design_full <- stats::model.matrix(full_formula,
                                       data=DRIMSeq::samples(drim_data))
    null_model <- stats::model.matrix(null_formula,
                                      data=DRIMSeq::samples(drim_data))
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
    if ((base::length(renamed_transcripts) == 1)
        && (base::is.na(renamed_transcripts))) {
        write_output(drim_data, out_dir_path)
    } else {
        write_output_with_rename(drim_data, renamed_transcripts, out_dir_path)
    }
}
