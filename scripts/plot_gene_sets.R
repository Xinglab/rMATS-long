library('ggVennDiagram')
library('ggplot2')

## Rscript plot_gene_sets.R genes_1.tsv genes_2.tsv 'Group 1 Name' 'Group 2 Name' out_file.png
##
## Each .tsv file needs to have 'gene_id' as one of the headers
args <- base::commandArgs(trailingOnly=TRUE)
set_1_path <- args[1]
set_2_path <- args[2]
set_1_name <- args[3]
set_2_name <- args[4]
out_path <- args[5]

get_colors <- function() {
  return(base::list(red='#FF0000',
                    white='#FFFFFF'))
}

read_gene_set_from_file <- function(path) {
  df <- utils::read.table(path, header=TRUE, sep='\t', check.names=FALSE)
  genes <- df$gene_id
  unique_genes <- base::unique(genes)
  return(unique_genes)
}

create_ggplot_theme <- function() {
  return(
    ggplot2::theme(panel.background=ggplot2::element_blank(),
                   axis.line=ggplot2::element_blank(),
                   axis.text=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   plot.title=ggplot2::element_blank(),
                   legend.text=ggplot2::element_blank())
  )
}

make_plot <- function(set_1, set_2, name_1, name_2) {
  colors <- get_colors()
  data <- base::list()
  data[[name_1]] <- set_1
  data[[name_2]] <- set_2

  venn <- ggVennDiagram::Venn(data)
  ## '201' is two horizontally overlapping circles
  processed <- ggVennDiagram::process_data(venn, shape_id='201')
  ## 'both' gives counts and percentage
  ## label_alpha=0 removes the label background
  plot <- ggVennDiagram::plot_venn(processed, label='both', label_alpha=0) +
      create_ggplot_theme() +
      ggplot2::scale_fill_gradient(low=colors$white, high=colors$red,
                                   limits=base::c(0, NA)) +
      ggplot2::labs(fill='Genes')

  return(plot)
}

make_empty_plot <- function(set_1_size, set_2_size, name_1, name_2) {
  if ((set_1_size > 0) && (set_2_size > 0)) {
    stop('make_empty_plot called with 2 non-empty sets')
  }

  colors <- get_colors()
  data <- base::list()
  ## Create sets so that:
  ##   * all regions (unique to 1, unique to 2, and overlap) have 1 element
  ##   * the non-zero set size (if any) has an extra element
  ## The regions with 1 value will be colored white.
  ## Any region with more than 1 value will be colored red.
  if ((set_1_size == 0) && (set_2_size == 0)) {
    set_1 <- base::c(1, 2)
    set_2 <- base::c(2, 3)
  } else if (set_1_size == 0) {
    set_1 <- base::c(1, 2)
    set_2 <- base::c(2, 3, 4)
  } else {
    set_1 <- base::c(1, 2, 3)
    set_2 <- base::c(3, 4)
  }
  data[[name_1]] <- set_1
  data[[name_2]] <- set_2

  venn <- ggVennDiagram::Venn(data)
  ## '201' is two horizontally overlapping circles
  processed <- ggVennDiagram::process_data(venn, shape_id='201')
  ## Create labels based on set_1_size and set_2_size
  label_data <- processed$regionLabel
  name_1_x <- label_data$X[label_data$name == name_1]
  name_1_y <- label_data$Y[label_data$name == name_1]
  name_2_x <- label_data$X[label_data$name == name_2]
  name_2_y <- label_data$Y[label_data$name == name_2]
  overlap_index <- (label_data$name != name_1) & (label_data$name != name_2)
  overlap_x <- label_data$X[overlap_index]
  overlap_y <- label_data$Y[overlap_index]
  labels <- base::data.frame(x=base::c(name_1_x, name_2_x, overlap_x),
                             y=base::c(name_1_y, name_2_y, overlap_y),
                             label=base::c(base::paste(set_1_size),
                                           base::paste(set_2_size),
                                           '0'))

  plot <- ggVennDiagram::plot_venn(processed, label='none') +
      create_ggplot_theme() +
      ggplot2::scale_fill_gradient(low=colors$white, high=colors$red,
                                   limits=base::c(1, 2)) +
      ggplot2::labs(fill='Genes') +
      ggplot2::geom_text(mapping=ggplot2::aes(x=x, y=y, label=label),
                         data=labels)

  return(plot)
}

save_plot <- function(plot, out_path) {
  ggplot2::ggsave(out_path, plot=plot, dpi=300, width=10, height=8)
}

main <- function() {
  gene_set_1 <- read_gene_set_from_file(set_1_path)
  gene_set_2 <- read_gene_set_from_file(set_2_path)
  set_1_size <- base::length(gene_set_1)
  set_2_size <- base::length(gene_set_2)
  if ((set_1_size == 0) || (set_2_size == 0)) {
    plot <- make_empty_plot(set_1_size, set_2_size, set_1_name, set_2_name)
  } else {
    plot <- make_plot(gene_set_1, gene_set_2, set_1_name, set_2_name)
  }

  save_plot(plot, out_path)
}

main()
