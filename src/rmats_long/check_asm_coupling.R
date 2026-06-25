## Rscript check_asm_coupling.R coupling_out.tsv coupling_stat_out.tsv [min_reads_per_group]
library(emmeans)
library(lme4)

args <- base::commandArgs(trailingOnly=TRUE)
in_path <- args[1]
out_path <- args[2]

if (base::length(args) >= 3) {
  min_reads_per_group_str <- args[3]
  min_reads_per_group <- base::as.integer(min_reads_per_group_str)
} else {
  min_reads_per_group <- 10
}

read_input_file <- function(path) {
  df <- utils::read.table(path, header=TRUE, sep='\t', check.names=FALSE)
  return(df)
}

write_to_out_file <- function(df, path) {
  utils::write.table(df, path, quote=FALSE, sep='\t', row.names=FALSE)
}

get_all_group_pvalue <- function(model) {
  joint_result <- tryCatch(
      emmeans::joint_tests(model),
      error=function(e) return(NULL))
  if (base::is.null(joint_result)) {
    return(NA)
  }

  all_group_pvalue <- joint_result$p.value[
      joint_result[['model term']] == 'group:is_high1']

  return(all_group_pvalue)
}

round_by_largest_remainder <- function(values) {
  floored <- base::floor(values)
  remainders <- values - floored
  extra <- base::sum(remainders)
  if (extra < 1) {
   return(floored)
  }

  rounded <- floored
  ordered <- base::order(remainders, decreasing=TRUE)
  for (extra_i in 1:extra) {
    index <- ordered[extra_i]
    rounded[index] <- rounded[index] + 1
  }

  return(rounded)
}

check_if_enough_group_reads <- function(min_reads_per_group, glmer_df,
                                        unique_groups) {
  for (group in unique_groups) {
    group_select <- glmer_df$group == group
    high2_sum <- base::sum(glmer_df$high2[group_select])
    low2_sum <- base::sum(glmer_df$low2[group_select])
    total <- high2_sum + low2_sum
    if (total < min_reads_per_group) {
      return(FALSE)
    }
  }

  return(TRUE)
}

run_test <- function(min_reads_per_group, asm_id_1, asm_id_2, asm_pair_df) {
  num_rows <- base::nrow(asm_pair_df)
  new_num_rows <- num_rows * 2
  sample_ids <- base::vector(mode='character', length=new_num_rows)
  groups <- base::vector(mode='character', length=new_num_rows)
  is_high1_values <- base::vector(mode='logical', length=new_num_rows)
  high2_values <- base::vector(mode='numeric', length=new_num_rows)
  low2_values <- base::vector(mode='numeric', length=new_num_rows)
  new_row_i <- 1
  for (row_i in 1:num_rows) {
     new_row_i_2 <- new_row_i + 1
     sample_id <- asm_pair_df$sample_id[row_i]
     sample_ids[new_row_i] <- sample_id
     sample_ids[new_row_i_2] <- sample_id

     group <- asm_pair_df$group[row_i]
     groups[new_row_i] <- group
     groups[new_row_i_2] <- group

     is_high1_values[new_row_i] <- TRUE
     is_high1_values[new_row_i_2] <- FALSE

     high_high <- asm_pair_df$high1_high2[row_i]
     low_high <- asm_pair_df$low1_high2[row_i]
     high_low <- asm_pair_df$high1_low2[row_i]
     low_low <- asm_pair_df$low1_low2[row_i]
     ## Add a small amount to avoid extreme values from the model
     to_add <- 1
     rounded <- round_by_largest_remainder(
         base::c(high_high + to_add,
                 low_high + to_add,
                 high_low + to_add,
                 low_low + to_add))

     high2_values[new_row_i] <- rounded[1]
     high2_values[new_row_i_2] <- rounded[2]
     low2_values[new_row_i] <- rounded[3]
     low2_values[new_row_i_2] <- rounded[4]

     new_row_i <- new_row_i + 2
  }

  unique_groups <- base::unique(groups)
  num_groups <- base::length(unique_groups)
  error_df <- base::data.frame(
      asm_id_1=asm_id_1, asm_id_2=asm_id_2,
      group=unique_groups, odds_ratio=NA, lower_95_ci=NA, upper_95_ci=NA,
      one_group_pvalue=NA, all_group_pvalue=NA)

  glmer_df <- base::data.frame(sample_id=sample_ids, group=groups,
                               is_high1=is_high1_values, high2=high2_values,
                               low2=low2_values)
  has_enough_reads <- check_if_enough_group_reads(min_reads_per_group, glmer_df,
                                                  unique_groups)
  if (!has_enough_reads) {
    return(error_df)
  }

  formula <- cbind(high2, low2) ~ group + group:is_high1 + (1|sample_id)
  control <- lme4::glmerControl(optimizer='bobyqa')
  family <- stats::binomial
  ## suppressMessages to ignore: "boundary (singular) fit: see help('isSingular')"
  model <- tryCatch(
    base::suppressMessages(
      lme4::glmer(formula, data=glmer_df, family=family, control=control)),
    error=function(e) return(NULL))

  if (base::is.null(model)) {
    return(error_df)
  }

  em_obj <- tryCatch(
    emmeans::emmeans(model, revpairwise ~ is_high1 | group, type='response'),
    error=function(e) return(NULL))

  if (base::is.null(em_obj)) {
    return(error_df)
  }

  contrast_ci <- stats::confint(em_obj$contrasts, level=0.95)
  contrast_summary <- base::summary(em_obj$contrasts)
  all_group_pvalue <- get_all_group_pvalue(model)

  ratios <- base::vector(mode='numeric', length=num_groups)
  lower_95_cis <- base::vector(mode='numeric', length=num_groups)
  upper_95_cis <- base::vector(mode='numeric', length=num_groups)
  one_group_pvalues <- base::vector(mode='numeric', length=num_groups)
  for (group_i in 1:num_groups) {
    group <- unique_groups[group_i]
    conf_row <- contrast_ci[contrast_ci$group == group, ]
    ratios[group_i] <- conf_row$odds.ratio
    lower_95_cis[group_i] <- conf_row$asymp.LCL
    upper_95_cis[group_i] <- conf_row$asymp.UCL
    one_group_pvalues[group_i] <- contrast_summary$p.value[
        contrast_summary$group == group]
  }

  return(base::data.frame(
    asm_id_1=asm_id_1, asm_id_2=asm_id_2,
    group=unique_groups, odds_ratio=ratios, lower_95_ci=lower_95_cis,
    upper_95_ci=upper_95_cis, one_group_pvalue=one_group_pvalues,
    all_group_pvalue=all_group_pvalue))
}

run_test_on_each_asm_pair <- function(min_reads_per_group, df) {
  num_rows <- base::nrow(df)
  if (num_rows == 0) {
    return(base::data.frame(
      asm_id_1=NA, asm_id_2=NA, group=NA, odds_ratio=NA, lower_95_ci=NA,
      upper_95_ci=NA, one_group_pvalue=NA, one_group_fdr=NA,
      all_group_pvalue=NA, all_group_fdr=NA))
  }

  results <- base::list()
  all_group_pvalues <- base::vector(mode='numeric')
  group_start_i <- 1
  group_asm_1 <- df$asm_id_1[1]
  group_asm_2 <- df$asm_id_2[1]
  for (row_i in 1:num_rows) {
    found_asm_1 <- df$asm_id_1[row_i]
    found_asm_2 <- df$asm_id_2[row_i]
    if ((found_asm_1 == group_asm_1) && (found_asm_2 == group_asm_2)) {
      next
    }

    pair_df <- df[group_start_i:(row_i-1), ]
    result <- run_test(min_reads_per_group, group_asm_1, group_asm_2, pair_df)
    ## Need to wrap result in a list() to append as a single value
    results <- base::append(results, base::list(result))
    all_group_pvalues <- base::append(all_group_pvalues,
                                      result$all_group_pvalue[1])
    group_start_i <- row_i
    group_asm_1 <- found_asm_1
    group_asm_2 <- found_asm_2
  }

  pair_df <- df[group_start_i:num_rows, ]
  result <- run_test(min_reads_per_group, group_asm_1, group_asm_2, pair_df)
  results <- base::append(results, base::list(result))
  all_group_pvalues <- base::append(all_group_pvalues,
                                    result$all_group_pvalue[1])
  all_group_fdr <- stats::p.adjust(all_group_pvalues, method='BH')

  num_results <- base::length(results)
  for (result_i in 1:num_results) {
    results[[result_i]]$all_group_fdr <- all_group_fdr[result_i]
  }
  out_df <- base::do.call(base::rbind, results)

  out_df$one_group_fdr <- stats::p.adjust(out_df$one_group_pvalue, method='BH')
  ## reorder pvalue columns
  all_group_pvalues <- out_df$all_group_pvalue
  all_group_fdrs <- out_df$all_group_fdr
  out_df$all_group_pvalue <- NULL
  out_df$all_group_fdr <- NULL
  out_df$all_group_pvalue <- all_group_pvalues
  out_df$all_group_fdr <- all_group_fdrs
  return(out_df)
}

main <- function() {
  emmeans::emm_options(msg.nesting=FALSE)
  in_df <- read_input_file(in_path)
  out_df <- run_test_on_each_asm_pair(min_reads_per_group, in_df)
  write_to_out_file(out_df, out_path)
}

main()
