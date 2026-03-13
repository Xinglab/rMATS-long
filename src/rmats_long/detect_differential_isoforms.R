library('this.path')

base::source(this.path::here('run_drimseq.R'))

args <- base::commandArgs(trailingOnly=TRUE)
num_args = base::length(args)
espresso_abundance_path <- args[1]
out_dir_path <- args[2]
num_threads <- base::as.integer(args[3])
group_1_file_path <- args[4]
group_2_file_path <- args[5]
covar_path <- ''
if (num_args == 6) {
    covar_path <- args[6]
}

read_tsv_with_rows_as_first_column <- function(path, first_column_name) {
    ## check.names=FALSE prevents '-' (dash) from being replaced with '.' in
    ## the sample names.
    data <- utils::read.table(path, header=TRUE, sep='\t', quote='',
                              row.names=1, check.names=FALSE)
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

read_covar_tsv <- function(covar_path) {
    if (covar_path == '') {
        return(base::data.frame())
    }

    df <- utils::read.table(covar_path, header=TRUE, sep='\t', quote='',
                            check.names=FALSE)
    return(df)
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

make_counts_data_frame <- function(abundance, used_samples) {
    ## rows are transcripts
    ## gene_id column (filter out NA genes)
    ## feature_id column with transcript IDs
    ## column for each sample with counts
    not_na_gene_bools <- !base::is.na(abundance$data$gene_ID)
    gene_ids <- abundance$data$gene_ID[not_na_gene_bools]
    transcript_ids <- base::rownames(abundance$data)[not_na_gene_bools]
    df <- base::data.frame(gene_id=gene_ids, feature_id=transcript_ids)
    for (sample in used_samples) {
        df[[sample]] <- abundance$data[[sample]][not_na_gene_bools]
    }

    return(split_multi_gene_rows(df))
}

main <- function() {
    ## Read and filter input data
    abundance <- read_abundance(espresso_abundance_path)
    sample_df <- make_sample_group_data_frame(abundance$sample_names,
                                              group_1_file_path,
                                              group_2_file_path)
    counts_and_renamed_transcripts <- make_counts_data_frame(abundance, sample_df$sample_id)
    counts_df <- counts_and_renamed_transcripts$counts_df
    renamed_transcripts <- counts_and_renamed_transcripts$renamed_transcripts
    covar_df <- read_covar_tsv(covar_path)

    run_drimseq(counts_df, sample_df, renamed_transcripts, num_threads,
                covar_df, out_dir_path)
}

main()
