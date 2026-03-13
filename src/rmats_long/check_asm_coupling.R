## Rscript check_asm_coupling.R coupling_out.tsv coupling_updated_out.tsv
args <- base::commandArgs(trailingOnly=TRUE)
in_path <- args[1]
out_path <- args[2]

read_input_file <- function(path) {
  df <- utils::read.table(path, header=TRUE, sep='\t', check.names=FALSE)
  return(df)
}

write_to_out_file <- function(df, path) {
  utils::write.table(df, path, quote=FALSE, sep='\t', row.names=FALSE)
}

run_fisher_exact <- function(main_main, main_other, other_main, other_other) {
  values <- base::c(main_main, main_other, other_main, other_other)
  table <- base::matrix(data=values, nrow=2, ncol=2)
  result <- stats::fisher.test(table)
  p_value <- result$p.value
  odds_ratio <- result$estimate
  return(base::list(p_value=p_value, odds_ratio=odds_ratio))
}

calculate_new_columns <- function(df) {
  num_rows <- base::nrow(df)
  num_useds <- base::vector(mode='numeric', length=num_rows)
  num_not_useds <- base::vector(mode='numeric', length=num_rows)
  psi_1s <- base::vector(mode='numeric', length=num_rows)
  psi_2s <- base::vector(mode='numeric', length=num_rows)
  odds_ratios <- base::vector(mode='numeric', length=num_rows)
  p_values <- base::vector(mode='numeric', length=num_rows)
  for (i in 1:num_rows) {
    mm <- df$main_main[i]
    mo <- df$main_other[i]
    om <- df$other_main[i]
    oo <- df$other_other[i]
    result <- run_fisher_exact(mm, mo, om, oo)
    p_values[i] = result$p_value
    odds_ratios[i] = result$odds_ratio
    num_used <- mm + mo + om + oo
    num_not_used <- df$unclear[i] + df$only_1[i] + df$only_2[i]
    num_useds[i] <- num_used
    num_not_useds[i] <- num_not_used
    psi_1s[i] <- (mm + mo) / num_used
    psi_2s[i] <- (mm + om) / num_used
  }

  df$num_used <- num_useds
  df$num_not_used <- num_not_useds
  df$psi_1 <- psi_1s
  df$psi_2 <- psi_2s
  df$odds_ratio <- odds_ratios
  df$p_value <- p_values
  return(df)
}


main <- function() {
  df <- read_input_file(in_path)
  df <- calculate_new_columns(df)
  keep_rows <- ((df$num_used > 3*df$num_not_used)
                & (df$num_used >= 10)
                & (df$psi_1 < 0.9)
                & (df$psi_1 > 0.1)
                & (df$psi_2 < 0.9)
                & (df$psi_2 > 0.1))
  df <- df[keep_rows, ]
  df$fdr <- stats::p.adjust(df$p_value, method='BH')
  write_to_out_file(df, out_path)
}

main()
