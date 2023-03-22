## Based on code from Yang Xu
library(cowplot)
library(ggplot2)

args <- base::commandArgs(trailingOnly=TRUE)
cpm_and_proportion_path <- args[1]
max_transcripts_str <- args[2]
transcript_color_string <- args[3]
group_color_string <- args[4]
out_plot_path <- args[5]

parse_color_strings <- function(color_string) {
    strings <- base::strsplit(color_string, ',')[[1]]
    return(strings)
}

read_cpm_and_proportion_file <- function(file_path) {
    data <- utils::read.table(file_path, header=TRUE, sep='\t')
    return(data)
}

plot_proportion <- function(data, group_colors, transcript_colors,
                            sample_colors, font_size) {
    ## Convert to percent (100% instead of 1.0).
    data$proportion <- data$proportion * 100
    plot <- ggplot2::ggplot(data=data,
                            ggplot2::aes(x=sample, y=proportion,
                                         fill=transcript, group=transcript)) +
        ## The geom_point is not plotted, but is used to get a legend.
        ## (y=NA_integer_) is needed instead of NA which has logical type.
        ## The legend will use the point size and shape (15=solid square).
        ## The group_colors are used in the theme to color the x-axis labels.
        ggplot2::geom_point(aes(y=NA_integer_, color=group), size=6, shape=15) +
        ggplot2::scale_color_manual(name='Group', values=group_colors) +
        ggplot2::geom_bar(stat='identity',
                          position=ggplot2::position_stack(vjust=1), width=0.9,
                          color='black', linewidth=0.5) +
        ggplot2::scale_fill_manual(name='Transcript ID',
                                   values=transcript_colors) +
        create_prop_theme(font_size, sample_colors) +
        ggplot2::labs(x='', y='Isoform proportion (%)', title='')
    return(plot)
}

plot_cpm <- function(data, group_colors, transcript_colors, font_size) {
    plot <- ggplot2::ggplot(data=data,
                            ggplot2::aes(x=sample, y=cpm, fill=transcript,
                                         group=transcript)) +
        ## The CPM plot has the bars stacked in reverse order.
        ggplot2::geom_bar(stat='identity',
                          position=ggplot2::position_stack(vjust=1,
                                                           reverse=TRUE),
                          width=0.9, color='black', linewidth=0.5) +
        ggplot2::scale_fill_manual(values=transcript_colors) +
        create_cpm_theme(font_size) +
        ggplot2::labs(x='', y='Abundance (CPM)', title='')
    return(plot)
}

combine_plots <- function(prop_plot, cpm_plot, font_size) {
    prop_with_legend_plot <- prop_plot +
        guides(fill=ggplot2::guide_legend(nrow=2, title.position='top',
                                          label.theme=ggplot2::element_text(size=font_size-1.5),
                                          keywidth=1, keyheight=1, order=1)) +
        guides(color=ggplot2::guide_legend(nrow=2, title.position='top',
                                           label.theme=ggplot2::element_text(size=font_size-1.5),
                                           keywidth=1, keyheight=1, order=2)) +
        ggplot2::theme(legend.position='bottom')

    legend_plot <- cowplot::get_legend(prop_with_legend_plot)
    box_plot <- cowplot::plot_grid(cpm_plot, prop_plot, align='v', ncol=1,
                                   rel_heights=c(1, 1.8))
    combined_plot <- cowplot::plot_grid(box_plot, legend_plot, ncol=1,
                                        align='v', axis='t',
                                        rel_heights=c(1, 0.2))
    return(combined_plot)
}

save_plot <- function(num_samples, combined_plot, out_plot_path) {
    out_width <- base::max(0.25*num_samples + 1, 8)
    ggplot2::ggsave(out_plot_path, plot=combined_plot, width=out_width,
                    height=5.8, dpi=300)
}

create_ggplot_theme <- function(font_size) {
    return(
        ggplot2::theme(panel.grid.major=ggplot2::element_line(color='white',
                                                              linewidth=0.01),
                       panel.grid.minor=ggplot2::element_line(color='white',
                                                              linewidth=0.01),
                       axis.title.y=ggplot2::element_text(size=font_size),
                       axis.text.y=ggplot2::element_text(size=font_size+2),
                       plot.title=ggplot2::element_blank(),
                       legend.position='none',
                       panel.border=ggplot2::element_rect(fill=NA,
                                                          linewidth=0.3,
                                                          color='black'),
                       panel.background=ggplot2::element_blank(),
                       legend.box='horizontal')
    )
}

create_prop_theme <- function(font_size, sample_colors) {
    return(
        create_ggplot_theme(font_size) +
        ggplot2::theme(axis.title.x=ggplot2::element_text(size=font_size),
                       axis.text.x=ggplot2::element_text(size=font_size+1,
                                                         color=sample_colors,
                                                         angle=90, vjust=0.5,
                                                         hjust=1),
                       legend.title=ggplot2::element_text(size=font_size-1,
                                                          hjust=0.5),
                       legend.text=ggplot2::element_text(size=font_size-1.5),
                       legend.direction='vertical',
                       plot.margin=unit(c(0.5, 5, 0.5, 5), 'pt'))
    )
}

create_cpm_theme <- function(font_size) {
    return(
        create_ggplot_theme(font_size) +
        ggplot2::theme(axis.title.x=ggplot2::element_blank(),
                       axis.ticks.x=ggplot2::element_blank(),
                       axis.text.x=ggplot2::element_blank(),
                       legend.title=ggplot2::element_text(size=font_size-1),
                       legend.text=ggplot2::element_text(size=font_size-2),
                       legend.direction='horizontal',
                       plot.margin=unit(c(5,5,0.5,5), 'pt'))
    )
}

get_unique_sample_colors <- function(unique_samples, unique_groups,
                                     group_colors, data) {
    unique_sample_colors <- base::sapply(unique_samples, function(sample) {
        first_i_for_sample <- base::which(data$sample == sample)[1]
        group <- data$group[first_i_for_sample]
        group_i <- base::which(unique_groups == group)
        return(group_colors[group_i])
    })
    return(unique_sample_colors)
}


rename_transcripts_and_define_factors <- function(data, max_transcripts,
                                                  unique_samples, unique_groups) {
    ## The input file had the transcripts in sorted order
    sorted_transcripts <- base::unique(data$transcript)
    num_transcripts <- base::length(sorted_transcripts)
    ## Combine additional transcripts into a single bar labeled 'Others"
    named_transcripts <- sorted_transcripts[1:max_transcripts]
    other_transcripts <- sorted_transcripts[(max_transcripts + 1):num_transcripts]
    data$transcript <- base::lapply(data$transcript, function(transcript) {
        if (transcript %in% other_transcripts) {
            return('Others')
        }
        return(transcript)
    })
    new_df <- data[data$transcript != 'Others', ]
    for (sample in unique_samples) {
        sample_indices <- base::which(data$sample == sample &
                                      data$transcript == 'Others')
        cpm_sum <- base::sum(data$cpm[sample_indices])
        prop_sum <- base::sum(data$proportion[sample_indices])
        first_sample_row <- data[sample_indices[1], ]
        first_sample_row$cpm <- cpm_sum
        first_sample_row$proportion <- prop_sum
        new_df <- base::rbind(new_df, first_sample_row)
    }

    ## Define sorting for the plot with factors
    new_df$sample <- base::factor(new_df$sample, levels=unique_samples)
    named_and_other <- base::append(named_transcripts, 'Others')
    new_df$transcript <- base::factor(new_df$transcript, levels=named_and_other)
    new_df$group <- base::factor(new_df$group, levels=unique_groups)
    return(new_df)
}


main <- function() {
    max_transcripts <- base::as.integer(max_transcripts_str)
    transcript_colors <- parse_color_strings(transcript_color_string)
    group_colors <- parse_color_strings(group_color_string)

    data <- read_cpm_and_proportion_file(cpm_and_proportion_path)
    unique_groups <- base::unique(data$group)
    unique_samples <- base::unique(data$sample)
    num_samples <- base::length(unique_samples)
    unique_sample_colors <- get_unique_sample_colors(unique_samples,
                                                     unique_groups,
                                                     group_colors, data)
    data <- rename_transcripts_and_define_factors(data, max_transcripts,
                                                  unique_samples, unique_groups)
    font_size <- 10
    prop_plot <- plot_proportion(data, group_colors, transcript_colors,
                                 unique_sample_colors, font_size)
    cpm_plot <- plot_cpm(data, group_colors, transcript_colors, font_size)
    combined <- combine_plots(prop_plot, cpm_plot, font_size)
    save_plot(num_samples, combined, out_plot_path)
}

main()
