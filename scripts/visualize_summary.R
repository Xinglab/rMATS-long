## Based on code from Yongjun Li
library(ggplot2)
library(ggrepel)

args <- base::commandArgs(trailingOnly=TRUE)
summary_path <- args[1]
out_dir <- args[2]

read_summary <- function(summary_path) {
  summary_lines <- base::readLines(summary_path)

  event_keys <- base::c('exon skipping', "alternative 5'-splice site",
                        "alternative 3'-splice site", 'mutually exclusive exons',
                        'intron retention', 'alternative first exon',
                        'alternative last exon', 'complex', 'combinatorial')
  categories <- base::vector(mode='character')
  counts <- base::vector(mode='integer')
  for (line in summary_lines) {
    line_length <- base::length(line)
    if (line_length == 0 || base::substr(line, 1, 1) == '#') {
      next
    }
    parts <- base::strsplit(line, ': ', fixed=TRUE)[[1]]
    if (base::length(parts) != 2) {
      next
    }
    key <- parts[1]
    value <- parts[2]
    if (key %in% event_keys) {
      categories <- base::append(categories, key)
      count <- base::as.integer(value)
      counts <- base::append(counts, count)
    }
  }

  df <- base::data.frame(category=categories, count=counts)
  return(df)
}

create_ggplot_theme <- function() {
  return(
    ggplot2::theme(panel.background=ggplot2::element_blank(),
                   axis.line=ggplot2::element_blank(),
                   axis.text=ggplot2::element_blank(),
                   axis.ticks=ggplot2::element_blank(),
                   plot.title=ggplot2::element_blank(),
                   legend.text=ggplot2::element_text(size=20))
  )
}

make_event_plot <- function(df) {
  colors <- base::c('#1E79B3',  ## blue
                    '#A5CDE4',  ## pale blue
                    '#B4D789',  ## pale green
                    '#E42A89',  ## pink
                    '#E4AB23',  ## orange
                    '#199C77',  ## blue green
                    '#F06046',  ## light red
                    '#7670B2',  ## purple
                    '#686868')  ## grey

  total_count <- base::sum(df$count)
  df$percent <- base::round(df$count / total_count, digits=3)
  cumulative <- base::cumsum(df$percent)
  label_pos <- base::vector(mode='numeric')
  for (i in 1:base::length(df$percent)) {
    offset <- 1
    if (i > 1) {
      offset <- 1 - cumulative[i - 1]
    }
    offset <- offset - (df$percent[i] / 2)
    label_pos <- base::append(label_pos, offset)
  }
  df$label_pos <- label_pos

  base::rownames(df) <- df$category
  df$category <- base::paste0(df$category, ' (', df$count, ')')
  df$category <- base::factor(df$category, levels=df$category)

  slice_labels <- base::sapply(df$percent, function(percent) {
    if (percent == 0) {
      return('')
    }
    return(base::paste0(percent * 100, '%'))
  })

  plot <- ggplot2::ggplot(data=df, ggplot2::aes(x=1, y=percent,
                                                fill=category)) +
    ggplot2::geom_bar(width=1, stat='identity') +
    ggplot2::coord_polar('y') +
    ggrepel::geom_text_repel(ggplot2::aes(x=1.1, y=label_pos, label=slice_labels),
                             nudge_x=0.15,
                             size=6) +
    create_ggplot_theme() +
    ggplot2::labs(x='', y='', fill='', title='') +
    ggplot2::guides(fill=ggplot2::guide_legend()) +
    ggplot2::scale_fill_manual(values=colors)

  return(plot)
}

save_plot <- function(plot, out_dir) {
  out_path_pdf <- base::paste0(out_dir, '/summary_plot.pdf')
  ggplot2::ggsave(out_path_pdf, plot=plot, dpi=300, width=10, height=8)

  out_path_png <- base::paste0(out_dir, '/summary_plot.png')
  ggplot2::ggsave(out_path_png, plot=plot, dpi=300, width=10, height=8)
}

main <- function() {
  df <- read_summary(summary_path)
  plot <- make_event_plot(df)
  save_plot(plot, out_dir)
}

main()
