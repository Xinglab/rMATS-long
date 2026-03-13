library('shiny')

SHINY_DIR <- base::getwd()
NOTIFICATION_SECONDS=10

get_last_part_of_path <- function(dir_path) {
  parts <- base::strsplit(dir_path, '/')[[1]]
  num_parts <- base::length(parts)
  last_part <- parts[num_parts]
  return(last_part)
}

get_run_mode <- function(main_dir) {
  last_part <- get_last_part_of_path(main_dir)
  parts <- base::strsplit(last_part, '_')[[1]]
  num_parts <- base::length(parts)
  if (num_parts < 3) {
    return('')
  }

  if (parts[1] == 'rmats' && parts[2] == 'long') {
    run_mode <- base::paste(parts[3:num_parts], collapse='_')
    return(run_mode)
  }

  return('')
}


get_parent_dir <- function(main_dir) {
  last_part <- get_last_part_of_path(main_dir)
  len_main_dir <- base::nchar(main_dir)
  len_last_part <- base::nchar(last_part)
  parent_dir <- base::substr(main_dir, 1, (len_main_dir - len_last_part) - 1)
  return(parent_dir)
}

get_work_dir <- function(rel_path, abs_path) {
  rel_parts <- base::strsplit(rel_path, '/')[[1]]
  abs_parts <- base::strsplit(abs_path, '/')[[1]]
  rel_len <- base::length(rel_parts)
  abs_len <- base::length(abs_parts)
  min_len <- min(rel_len, abs_len)
  i <- 0
  while (i < min_len) {
    rel_part <- rel_parts[rel_len - i]
    abs_part <- abs_parts[abs_len - i]
    if (rel_part == abs_part) {
      i <- i + 1
      next
    }

    work_dir <- base::paste(abs_parts[1:(abs_len - i)], collapse='/')
    return(work_dir)
  }

  if (i < abs_len) {
    work_dir <- base::paste(abs_parts[1:(abs_len - i)], collapse='/')
    return(work_dir)
  }

  return(FALSE)
}

convert_rel_to_abs <- function(rel_path, work_dir) {
  path <- base::paste0(work_dir, '/', rel_path)
  return(base::normalizePath(path))
}

find_relative_path <- function(target, main_dir, parent_dir) {
  if (target == '') {
    return('')
  }

  target_is_abs <- base::substr(target, 1, 1) == '/'
  main_is_abs <- base::substr(main_dir, 1, 1) == '/'
  if (target_is_abs || main_is_abs) {
    return(target)
  }

  found_parent <- get_parent_dir(main_dir)
  work_dir <- get_work_dir(found_parent, parent_dir)
  if (work_dir == FALSE) {
    return(target)
  }
  converted <- convert_rel_to_abs(target, work_dir)
  return(converted)
}


check_for_arguments_in_summary <- function(parent_dir, summary_txt) {
  result <- base::list(group_1='', group_2='', group_1_name='', group_2_name='',
                       out_dir='', abundance='', updated_gtf='', gencode_gtf='')

  lines <- base::readLines(summary_txt, n=1)
  parts <- strsplit(lines[1], ' ')[[1]]
  prev_part <- ''
  for (part in parts) {
    if (prev_part == '--group-1') {
      result$group_1 <- part
    }
    if (prev_part == '--group-2') {
      result$group_2 <- part
    }
    if (prev_part == '--group-1-name') {
      result$group_1_name <- part
    }
    if (prev_part == '--group-2-name') {
      result$group_2_name <- part
    }
    if (prev_part == '--out-dir') {
      result$out_dir <- part
    }
    if (prev_part == '--abundance') {
      result$abundance <- part
    }
    if (prev_part == '--updated-gtf') {
      result$updated_gtf <- part
    }
    if (prev_part == '--gencode-gtf') {
      result$gencode_gtf <- part
    }

    prev_part <- part
  }

  result$group_1 <- find_relative_path(result$group_1, result$out_dir,
                                       parent_dir)
  result$group_2 <- find_relative_path(result$group_2, result$out_dir,
                                       parent_dir)
  result$abundance <- find_relative_path(result$abundance, result$out_dir,
                                         parent_dir)
  result$updated_gtf <- find_relative_path(result$updated_gtf, result$out_dir,
                                           parent_dir)
  result$gencode_gtf <- find_relative_path(result$gencode_gtf, result$out_dir,
                                           parent_dir)

  return(result)
}

get_file_paths <- function(main_dir) {
  result <- base::list(main_dir='', diff_asms='', diff_isos='', summary_plot='',
                       summary_txt='', align_dir='', event_dir='', count_tsv='',
                       combined_gtf='', sample_total_tsv='', group_1='',
                       group_2='', group_1_name='', group_2_name='', abundance='',
                       updated_gtf='', gencode_gtf='')
  if (main_dir == '') {
    return(result)
  }

  result$main_dir <- base::normalizePath(main_dir)
  parent_dir <- get_parent_dir(result$main_dir)
  result$diff_asms <- base::paste0(result$main_dir, '/differential_asms.tsv')
  result$diff_isos <- base::paste0(result$main_dir, '/differential_isoforms.tsv')
  result$summary_plot <- base::paste0(result$main_dir, '/summary_plot.png')
  result$summary_txt <- base::paste0(result$main_dir, '/summary.txt')
  result$align_dir <- base::paste0(parent_dir, '/alignments')
  result$event_dir <- base::paste0(parent_dir, '/events')
  result$count_tsv <- base::paste0(result$main_dir, '/count.tsv')
  result$combined_gtf <- base::paste0(parent_dir, '/combined.gtf')
  result$sample_total_tsv <- base::paste0(result$align_dir,
                                          '/sample_read_totals.tsv')
  result$group_1 <- base::paste0(parent_dir, '/group_1.txt')
  result$group_2 <- base::paste0(parent_dir, '/group_2.txt')
  summary_args <- check_for_arguments_in_summary(parent_dir, result$summary_txt)
  result$abundance <- summary_args$abundance
  result$updated_gtf <- summary_args$updated_gtf
  result$gencode_gtf <- summary_args$gencode_gtf
  if (!base::file.exists(result$group_1)) {
    result$group_1 <- summary_args$group_1
  }
  if (!base::file.exists(result$group_2)) {
    result$group_2 <- summary_args$group_2
  }

  if (!base::file.exists(result$diff_asms)) {
    result$diff_asms <- base::paste0(result$main_dir, '/differential_genes.tsv')
    result$diff_isos <- base::paste0(result$main_dir,
                                     '/differential_transcripts.tsv')
  }

  if (!base::file.exists(result$event_dir)) {
    run_mode <- get_run_mode(result$main_dir)
    result$event_dir <- base::paste0(parent_dir, '/', run_mode, '_events')
  }

  keys <- base::names(result)
  for (key in keys) {
    value <- result[[key]]
    if (!base::file.exists(value)) {
      result[[key]] <- ''
    }
  }

  result$group_1_name <- 'group-1'
  if (summary_args$group_1_name != '') {
    result$group_1_name <- summary_args$group_1_name
  }
  result$group_2_name <- 'group-2'
  if (summary_args$group_2_name != '') {
    result$group_2_name <- summary_args$group_2_name
  }

  return(result)
}

read_text_file <- function(path) {
  lines <- base::readLines(path)
  text <- base::paste(lines, collapse='\n')
  return(text)
}

get_table_from_file <- function(path) {
  parsed <- utils::read.table(path, header=TRUE, sep='\t')
  return(parsed)
}

select_at_most_n_rows_of_table <- function(table, max_rows) {
  selected <- base::list()
  names <- base::names(table)
  for (name in names) {
    selected[[name]] <- table[[name]][1:max_rows]
  }

  return(selected)
}

sort_table_by_column <- function(table, column, decreasing) {
  names <- base::names(table)
  found_column <- base::sum(names == column) == 1
  if (!found_column) {
    return(table)
  }

  sort_i <- base::order(table[[column]], decreasing=decreasing)
  for (name in names) {
    table[[name]] <- table[[name]][sort_i]
  }

  return(table)
}

get_default_table <- function() {
  return(base::list(header=base::c('...')))
}

get_default_image <- function() {
  return (base::paste0(SHINY_DIR, '/clear_1px.png'))
}

print_and_run_command <- function(command, args, error_message) {
  base::cat(command, args, '\n')
  return_code <- base::system2(command, args)
  if (return_code != 0) {
    error <- base::paste0(error_message, ': ', return_code)
  } else {
    error <- FALSE
  }

  return(error)
}

print_and_run_command_with_stdout_redirect <- function(command, args, out_path,
                                                       error_message) {
  base::cat(command, args, '\n')
  return_code <- base::system2(command, args, stdout=out_path)
  if (return_code != 0) {
    error <- base::paste0(error_message, ': ', return_code)
  } else {
    error <- FALSE
  }

  return(error)
}

generate_summary_plot <- function(summary_txt, out_dir) {
  command <- 'rmats-long'
  script <- 'visualize_summary.R'
  args <- base::c(script, summary_txt, out_dir)
  error <- print_and_run_command(command, args,
                                 'command to generate summary plot failed')
  out_path <- base::paste0(out_dir, '/summary_plot.png')
  return(base::list(error=error, out_path=out_path))
}

create_abundance_for_asm <- function(asm_id, event_dir, sample_total_tsv,
                                     count_tsv, gene_out_dir) {
  out_name <- base::paste0('abun_', asm_id, '.esp')
  out_path <- base::paste0(gene_out_dir, '/', out_name)
  command <- 'rmats-long'
  script <- 'create_abundance_file_for_asm.py'
  args <- base::c(script, '--asm-id', asm_id, '--event-dir', event_dir,
                  '--sample-totals', sample_total_tsv, '--count-tsv', count_tsv,
                  '--out-file', out_path)
  error <- print_and_run_command(command, args,
                                 'command to get ASM abundance failed')
  return(base::list(error=error, out_path=out_path))
}

get_chr_i_from_asm_id <- function(asm_id) {
  parts <- base::strsplit(asm_id, '_')[[1]]
  return(parts[1])
}

create_gtf_for_asms <- function(event_dir, out_dir) {
  out_path <- base::paste0(out_dir, '/asm_as_gene.gtf')
  command <- 'rmats-long'
  script <- 'create_gtf_from_asm_definitions.py'
  args <- base::c(script, '--event-dir', event_dir, '--out-gtf', out_path)
  error <- print_and_run_command(command, args,
                                 'command to create ASM gtf failed')
  if (error != FALSE) {
    return(base::list(error=error, out_path=out_path))
  }

  command <- 'sed'
  sed_script <- '\'s/gene_id "[^"]*"; asm_id "\\([^"]*\\)";/gene_id "\\1"; asm_id "\\1";/\''
  args <- base::c('-i', sed_script, out_path)
  error <- print_and_run_command(command, args,
                                 'command to update gene_ids in ASM gtf failed')
  return(base::list(error=error, out_path=out_path))
}

update_differential_isoforms_headers <- function(diff_isos, out_dir) {
  out_path <- base::paste0(out_dir, '/diff_isos_modified_headers.tsv')
  command <- 'sed'
  sed_script <- '\'s/asm_id\tgene_id\tisoform_id/gene_id\torig_gene_id\tfeature_id/\''
  args <- base::c(sed_script, diff_isos)
  error <- print_and_run_command_with_stdout_redirect(
      command, args, out_path,
      'command to modifiy differential_isoforms.tsv headers failed')
  return(base::list(error=error, out_path=out_path))
}

get_gene_out_dir <- function(main_dir, gene_id) {
  gene_out_dir <- base::paste0(main_dir, '/results_by_gene/', gene_id)
  return(gene_out_dir)
}

get_graph_file_for_asm <- function(asm_id, event_dir) {
  chr_i <- get_chr_i_from_asm_id(asm_id)
  graph_file <- base::paste0(event_dir, '/graph_', chr_i, '.txt')
  return(graph_file)
}

generate_abundance_and_structure_plots_from_abundance <- function(
    gene_id, main_dir, abundance, updated_gtf, gencode_gtf, diff_isos,
    group_1, group_2, group_1_name, group_2_name, plot_file_type, intron_scaling,
    max_transcripts) {
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  dir.create(gene_out_dir, recursive=TRUE)
  abun_plot_path <- base::paste0(gene_out_dir, '/', gene_id, '_abundance.png')
  struct_plot_path <- base::paste0(gene_out_dir, '/', gene_id, '_structure.png')
  struct_sj_plot_path <- base::paste0(gene_out_dir, '/', gene_id,
                                      '_structure_sj_counts.png')
  colors_out_path <- base::paste0(gene_out_dir, '/', gene_id, '_colors.tsv')
  command <- 'rmats-long'
  script <- 'visualize_isoforms.py'
  args <- base::c(script, '--gene-id', gene_id, '--abundance',
                  abundance, '--updated-gtf', updated_gtf,
                  '--gencode-gtf', gencode_gtf,
                  '--diff-transcripts', diff_isos, '--out-dir',
                  gene_out_dir, '--group-1', group_1, '--group-2', group_2,
                  '--group-1-name', group_1_name, '--group-2-name', group_2_name,
                  '--plot-file-type', plot_file_type, '--intron-scaling',
                  intron_scaling, '--max-transcripts', max_transcripts,
                  '--out-transcript-colors', colors_out_path)
  error <- print_and_run_command(command, args,
                                 'command to get ASM abundance failed')
  return(base::list(error=error, abun_plot_path=abun_plot_path,
                    struct_plot_path=struct_plot_path,
                    struct_sj_plot_path=struct_sj_plot_path,
                    colors_out_path=colors_out_path,
                    asm_gtf_out_path=updated_gtf,
                    modified_diff_isos=diff_isos))
}

generate_abundance_and_structure_plots <- function(
    gene_id, asm_id, main_dir, event_dir, sample_total_tsv, count_tsv, diff_isos,
    group_1, group_2, group_1_name, group_2_name, plot_file_type, intron_scaling,
    max_transcripts) {
  graph_file <- get_graph_file_for_asm(asm_id, event_dir)
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  dir.create(gene_out_dir, recursive=TRUE)
  abun_plot_path <- base::paste0(gene_out_dir, '/', asm_id, '_abundance.png')
  struct_plot_path <- base::paste0(gene_out_dir, '/', asm_id, '_structure.png')
  struct_sj_plot_path <- base::paste0(gene_out_dir, '/', asm_id,
                                      '_structure_sj_counts.png')
  colors_out_path <- base::paste0(gene_out_dir, '/', asm_id, '_colors.tsv')
  abun_result <- create_abundance_for_asm(asm_id, event_dir, sample_total_tsv,
                                          count_tsv, gene_out_dir)
  if (abun_result$error != FALSE) {
    return(base::list(error=abun_result$error, abun_plot_path=abun_plot_path,
                      struct_plot_path=struct_plot_path,
                      struct_sj_plot_path=struct_sj_plot_path,
                      colors_out_path=colors_out_path, asm_gtf_out_path='',
                      modified_diff_isos=''))
  }

  gtf_result <- create_gtf_for_asms(event_dir, main_dir)
  if (gtf_result$error != FALSE) {
    return(base::list(error=gtf_result$error, abun_plot_path=abun_plot_path,
                      struct_plot_path=struct_plot_path,
                      struct_sj_plot_path=struct_sj_plot_path,
                      colors_out_path=colors_out_path,
                      asm_gtf_out_path=gtf_result$out_path,
                      modified_diff_isos=''))
  }

  diff_result <- update_differential_isoforms_headers(diff_isos, main_dir)
  if (diff_result$error != FALSE) {
    return(base::list(error=diff_result$error, abun_plot_path=abun_plot_path,
                      struct_plot_path=struct_plot_path,
                      struct_sj_plot_path=struct_sj_plot_path,
                      colors_out_path=colors_out_path,
                      asm_gtf_out_path=gtf_result$out_path,
                      modified_diff_isos=''))
  }

  command <- 'rmats-long'
  script <- 'visualize_isoforms.py'
  args <- base::c(script, '--gene-id', gene_id, '--abundance',
                  abun_result$out_path, '--updated-gtf', gtf_result$out_path,
                  '--diff-transcripts', diff_result$out_path, '--out-dir',
                  gene_out_dir, '--group-1', group_1, '--group-2', group_2,
                  '--group-1-name', group_1_name, '--group-2-name', group_2_name,
                  '--plot-file-type', plot_file_type, '--intron-scaling',
                  intron_scaling, '--max-transcripts', max_transcripts,
                  '--is-asm', '--graph-file', graph_file, '--asm-id', asm_id,
                  '--out-transcript-colors', colors_out_path)
  error <- print_and_run_command(command, args,
                                 'command to get ASM abundance failed')
  return(base::list(error=error, abun_plot_path=abun_plot_path,
                    struct_plot_path=struct_plot_path,
                    struct_sj_plot_path=struct_sj_plot_path,
                    colors_out_path=colors_out_path,
                    asm_gtf_out_path=gtf_result$out_path,
                    modified_diff_isos=diff_result$out_path))
}

file_move <- function(from, to) {
  ok <- base::file.copy(from, to, overwrite=TRUE)
  if (!ok) {
    return(base::paste0('failed to move ', from, ' to ', to))
  }

  ok <- base::file.remove(from)
  if (!ok) {
    return(base::paste0('failed to remove during move ', from, ' to ', to))
  }

  return(FALSE)
}

read_table_with_possible_empty_columns <- function(path) {
  lines <- base::readLines(path)
  headers <- NULL
  num_headers <- 0
  values <- NULL
  for (line in lines) {
    parts <- strsplit(line, '\t')[[1]]
    if (base::is.null(headers)) {
      headers <- parts
      num_headers <- base::length(headers)
      next
    }

    if (base::is.null(values)) {
      values <- base::list()
      for (i in 1:num_headers) {
        header <- headers[i]
        values[[header]] <- base::c(parts[i])
      }

      next
    }

    for (i in 1:num_headers) {
      header <- headers[i]
      values[[header]] <- base::append(values[[header]], parts[i])
    }
  }

  return(values)
}

read_asm_from_event <- function(asm_id, event_dir) {
  chr_i <- get_chr_i_from_asm_id(asm_id)
  event_file <- base::paste0(event_dir, '/chr_id_', chr_i, '.tsv')
  tmp_file <- base::paste0(event_dir, 'tmp_grep.txt')
  command <- 'grep'
  pattern <- base::paste0('\'^\\(asm_id\\|', asm_id, '\\)\t\'')
  args <- base::c(pattern, event_file)
  error <- print_and_run_command_with_stdout_redirect(
      command, args, tmp_file, 'grep for the asm ID failed')
  if (error != FALSE) {
    strand <- ''
    start <- ''
    end <- ''
  } else {
    parsed <- read_table_with_possible_empty_columns(tmp_file)
    base::file.remove(tmp_file)
    strand <- parsed$strand
    start <- parsed$start
    end <- parsed$end
  }

  return(base::list(error=error, strand=strand, start=start, end=end))
}

generate_in_gene_plot <- function(
    gene_id, asm_id, main_dir, event_dir, modified_diff_isos, asm_gtf,
    combined_gtf, plot_file_type, intron_scaling, max_transcripts) {
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  tmp_out_dir <- base::paste0(gene_out_dir, '/', 'tmp_in_gene')
  in_gene_plot_path <- base::paste0(gene_out_dir, '/', asm_id, '_in_gene.png')
  in_gene_sj_plot_path <- base::paste0(gene_out_dir, '/', asm_id,
                                       '_in_gene_sj_counts.png')
  in_gene_tmp_path <- base::paste0(tmp_out_dir, '/', gene_id, '_structure.png')
  in_gene_sj_tmp_path <- base::paste0(tmp_out_dir, '/', gene_id,
                                      '_structure_sj_counts.png')
  graph_file <- get_graph_file_for_asm(asm_id, event_dir)
  event_result <- read_asm_from_event(asm_id, event_dir)
  if (event_result$error != FALSE) {
    return(base::list(error=event_result$error,
                      in_gene_plot_path=in_gene_plot_path,
                      in_gene_sj_plot_path=in_gene_sj_plot_path))
  }

  if (event_result$strand == '-') {
    start_coord <- event_result$end
    end_coord <- event_result$start
  } else {
    start_coord <- event_result$start
    end_coord <- event_result$end
  }

  command <- 'rmats-long'
  script <- 'visualize_isoforms.py'
  args <- base::c(script, '--gene-id', gene_id, '--asm-id', asm_id, '--asm-gtf',
                  asm_gtf, '--gencode-gtf', combined_gtf, '--updated-gtf',
                  combined_gtf, '--diff-transcripts', modified_diff_isos,
                  '--out-dir', tmp_out_dir, '--plot-file-type', plot_file_type,
                  '--intron-scaling', intron_scaling, '--max-transcripts',
                  max_transcripts, '--graph-file', graph_file)
  if (start_coord != 'source' && start_coord != 'sink') {
    args <- base::append(args, base::c('--start-coord', start_coord))
  }
  if (end_coord != 'source' && end_coord != 'sink') {
    args <- base::append(args, base::c('--end-coord', end_coord))
  }

  error <- print_and_run_command(command, args,
                                 'command to plot ASM region in gene failed')
  if (error == FALSE) {
    error <- file_move(in_gene_tmp_path, in_gene_plot_path)
    if (error == FALSE) {
      error <- file_move(in_gene_sj_tmp_path, in_gene_sj_plot_path)
    }
  }

  return(base::list(error=error, in_gene_plot_path=in_gene_plot_path,
                    in_gene_sj_plot_path=in_gene_sj_plot_path))
}

generate_simple_graph_plot <- function(gene_id, asm_id, main_dir, asm_gtf,
                                       isoform_colors, intron_scaling,
                                       equal_spacing, color_first) {
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  out_path <- base::paste0(gene_out_dir, '/', asm_id, '_simple_graph.png')
  command <- 'rmats-long'
  script <- 'plot_simple_splice_graph.py'
  args <- base::c(script, '--out-file', out_path, '--gtf', asm_gtf,
                  '--intron-scaling', intron_scaling, '--transcript-colors',
                  isoform_colors, '--color-by-isoform', '--ids-from-color-file',
                  '--show-isoform-endpoint-symbols')
  if (equal_spacing) {
    args <- base::append(args, '--equal-spacing')
  }
  if (color_first) {
    args <- base::append(args, '--color-first-transcript-only')
  }

  error <- print_and_run_command(command, args,
                                 'command to create simple splice graph failed')
  return(base::list(error=error, out_path=out_path))
}

generate_graph_plot <- function(gene_id, asm_id, main_dir, event_dir) {
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  out_path <- base::paste0(gene_out_dir, '/', asm_id, '_graph.png')
  command <- 'rmats-long'
  script <- 'plot_splice_graph.py'
  args <- base::c(script, '--event-dir', event_dir, '--gene-id', gene_id,
                  '--asm-id', asm_id, '--out-file', out_path)
  error <- print_and_run_command(command, args,
                                 'command to create splice graph failed')
  return(base::list(error=error, out_path=out_path))
}

get_main_and_second_id_from_colors_tsv <- function(colors_tsv) {
  error <- FALSE
  main_id <- ''
  second_id <- ''
  if (!base::file.exists(colors_tsv)) {
    error <- base::paste0('colors_tsv (', colors_tsv, ') not found')
    return(base::list(error=error, main_id=main_id, second_id=second_id))
  }

  parsed <- utils::read.table(colors_tsv, header=TRUE, sep='\t')
  num_transcripts <- base::length(parsed$transcript_id)
  if (num_transcripts == 0) {
    error <- base::paste0('no transcripts found in ', colors_tsv)
  } else if (num_transcripts == 1) {
    main_id <- parsed$transcript_id[1]
  } else {
    main_id <- parsed$transcript_id[1]
    second_id <- parsed$transcript_id[2]
  }
  return(base::list(error=error, main_id=main_id, second_id=second_id))
}

generate_diff_table <- function(gene_id, asm_id, main_dir, asm_gtf,
                                isoform_colors) {
  gene_out_dir <- get_gene_out_dir(main_dir, gene_id)
  id_result <- get_main_and_second_id_from_colors_tsv(isoform_colors)
  if (id_result$error != FALSE) {
    return(base::list(error=id_result$error, out_path=''))
  }

  command <- 'rmats-long'
  script <- 'classify_isoform_differences.py'
  if (id_result$second_id == '') {
    out_path <- base::paste0(gene_out_dir, '/', asm_id,
                             '_isoform_differences_from_', id_result$main_id,
                             '.tsv')
    args <- base::c(script, '--updated-gtf', asm_gtf, '--out-tsv', out_path,
                    '--main-transcript-id', id_result$main_id)
  } else {
    out_path <- base::paste0(gene_out_dir, '/', asm_id, '_isoform_differences_',
                             id_result$main_id, '_to_', id_result$second_id,
                             '.tsv')
    args <- base::c(script, '--updated-gtf', asm_gtf, '--out-tsv', out_path,
                    '--main-transcript-id', id_result$main_id,
                    '--second-transcript-id', id_result$second_id)
  }
  error_text <- 'command to create splicing differences table failed'
  error <- print_and_run_command(command, args, error_text)

  return(base::list(error=error, out_path=out_path))
}

server <- function(input, output, session) {
  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Setting main directory...', session=session,
                              duration=NULL, id='setting_main')
      file_paths <- get_file_paths(input$main_directory)
      if (file_paths$diff_asms == '') {
        shiny::showNotification(
          base::paste0('No ASM file found in ', input$main_directory),
          session=session, duration=NOTIFICATION_SECONDS, type='error')
        shiny::removeNotification(id='setting_main', session=session)
        return()
      }

      shiny::updateTextInput(session=session, 'diff_asms',
                             value=file_paths$diff_asms)
      shiny::updateTextInput(session=session, 'diff_isos',
                             value=file_paths$diff_isos)
      shiny::updateTextInput(session=session, 'summary_plot',
                             value=file_paths$summary_plot)
      shiny::updateTextInput(session=session, 'summary_txt',
                             value=file_paths$summary_txt)
      shiny::updateTextInput(session=session, 'align_dir',
                             value=file_paths$align_dir)
      shiny::updateTextInput(session=session, 'event_dir',
                             value=file_paths$event_dir)
      shiny::updateTextInput(session=session, 'count_tsv',
                             value=file_paths$count_tsv)
      shiny::updateTextInput(session=session, 'combined_gtf',
                             value=file_paths$combined_gtf)
      shiny::updateTextInput(session=session, 'sample_total_tsv',
                             value=file_paths$sample_total_tsv)
      shiny::updateTextInput(session=session, 'group_1',
                             value=file_paths$group_1)
      shiny::updateTextInput(session=session, 'group_2',
                             value=file_paths$group_2)
      shiny::updateTextInput(session=session, 'group_1_name',
                             value=file_paths$group_1_name)
      shiny::updateTextInput(session=session, 'group_2_name',
                             value=file_paths$group_2_name)
      shiny::updateTextInput(session=session, 'abundance',
                             value=file_paths$abundance)
      shiny::updateTextInput(session=session, 'updated_gtf',
                             value=file_paths$updated_gtf)
      shiny::updateTextInput(session=session, 'gencode_gtf',
                             value=file_paths$gencode_gtf)

      default_image_path <- get_default_image()
      output$abundance_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='700',
             alt='default image'),
        deleteFile=FALSE)
      output$structure_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='333',
             alt='default image'),
        deleteFile=FALSE)
      output$structure_sj_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='333',
             alt='default image'),
        deleteFile=FALSE)
      output$in_gene_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='333',
             alt='default image'),
        deleteFile=FALSE)
      output$in_gene_sj_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='333',
             alt='default image'),
        deleteFile=FALSE)
      output$simple_graph_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='500',
             alt='default image'),
        deleteFile=FALSE)
      output$graph_image <- shiny::renderImage(
        base::list(src=default_image_path, width='1000', height='333',
             alt='default image'),
        deleteFile=FALSE)

      output$isoform_differences_table <- shiny::renderTable(get_default_table())

      if (base::file.exists(file_paths$summary_plot)) {
          output$summary_pie_image <- shiny::renderImage(
            base::list(src=file_paths$summary_plot, width='1000', height='800',
                 alt='summary_pie_plot'),
            deleteFile=FALSE)
      } else {
          output$summary_pie_image <- shiny::renderImage(
            base::list(src=default_image_path, width='1000', height='800',
                 alt='default image'),
            deleteFile=FALSE)
      }

      if (base::file.exists(file_paths$summary_txt)) {
        text_value <- read_text_file(file_paths$summary_txt)
        output$summary_text <- shiny::renderText(text_value, sep='')
      } else {
        output$summary_text <- shiny::renderText('')
      }

      if (base::file.exists(file_paths$diff_asms)) {
        asm_table <- get_table_from_file(file_paths$diff_asms)
        limited_table <- select_at_most_n_rows_of_table(asm_table,
                                                        input$asm_max_rows)
        output$asm_stat_table <- shiny::renderTable(limited_table)
        updateSelectInput(session=session, 'asm_sort_col_1',
                          choices=base::names(asm_table))
        updateSelectInput(session=session, 'asm_sort_col_2',
                          choices=base::names(asm_table))
      } else {
        output$asm_stat_table <- shiny::renderTable(get_default_table())
      }

      if (base::file.exists(file_paths$diff_isos)) {
        iso_table <- get_table_from_file(file_paths$diff_isos)
        limited_table <- select_at_most_n_rows_of_table(iso_table,
                                                        input$iso_max_rows)
        output$isoform_stat_table <- shiny::renderTable(limited_table)
        updateSelectInput(session=session, 'iso_sort_col_1',
                          choices=base::names(iso_table))
        updateSelectInput(session=session, 'iso_sort_col_2',
                          choices=base::names(iso_table))
      } else {
        output$isoform_stat_table <- shiny::renderTable(get_default_table())
      }

      bslib::nav_select('summary_navset', selected='summary_pie_tab',
                        session=session)
      bslib::nav_select('asm_navset', selected='asm_abundance_tab',
                        session=session)

      shiny::removeNotification(id='setting_main', session=session)
    }),
    input$set_main_directory
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Sorting ASM results...', session=session,
                              duration=NULL, id='sorting_asm')

      if (base::file.exists(input$diff_asms)) {
        table <- get_table_from_file(input$diff_asms)
        sort_col_1 <- input$asm_sort_col_1
        sort_col_2 <- input$asm_sort_col_2
        sort_col_1_dec <- input$asm_sort_col_1_decreasing
        sort_col_2_dec <- input$asm_sort_col_2_decreasing
        table <- sort_table_by_column(table, sort_col_2, sort_col_2_dec)
        table <- sort_table_by_column(table, sort_col_1, sort_col_1_dec)
        limited_table <- select_at_most_n_rows_of_table(table,
                                                        input$asm_max_rows)
        output$asm_stat_table <- shiny::renderTable(limited_table)
      } else {
        shiny::showNotification(
          'No ASM result file found',
          session=session, duration=NOTIFICATION_SECONDS, type='error')
      }
      shiny::removeNotification(id='sorting_asm', session=session)
    }),
    input$sort_asm_results
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Sorting isoform results...', session=session,
                              duration=NULL, id='sorting_iso')

      if (base::file.exists(input$diff_isos)) {
        table <- get_table_from_file(input$diff_isos)
        sort_col_1 <- input$iso_sort_col_1
        sort_col_2 <- input$iso_sort_col_2
        sort_col_1_dec <- input$iso_sort_col_1_decreasing
        sort_col_2_dec <- input$iso_sort_col_2_decreasing
        table <- sort_table_by_column(table, sort_col_2, sort_col_2_dec)
        table <- sort_table_by_column(table, sort_col_1, sort_col_1_dec)
        limited_table <- select_at_most_n_rows_of_table(table,
                                                        input$iso_max_rows)
        output$isoform_stat_table <- shiny::renderTable(limited_table)
      } else {
        shiny::showNotification(
          'No isoform result file found',
          session=session, duration=NOTIFICATION_SECONDS, type='error')
      }
      shiny::removeNotification(id='sorting_iso', session=session)
    }),
    input$sort_iso_results
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Finding gene ID...', session=session,
                              duration=NULL, id='finding_gene')

      found_gene_id <- ''
      if (base::file.exists(input$diff_isos)) {
        table <- get_table_from_file(input$diff_isos)
        matches <- table$asm_id == input$asm_id
        gene_ids <- table$gene_id[matches]
        if (base::length(gene_ids) == 0) {
          shiny::showNotification(
            base::paste0('"', input$asm_id, '" not found in isoform result file'),
            session=session, duration=NOTIFICATION_SECONDS, type='error')
        } else {
          found_gene_id <- gene_ids[1]
        }
      } else {
        shiny::showNotification(
          'No isoform result file found',
          session=session, duration=NOTIFICATION_SECONDS, type='error')
      }
      shiny::updateTextInput(session=session, 'gene_id',
                             value=found_gene_id)
      shiny::removeNotification(id='finding_gene', session=session)
    }),
    input$find_gene_id_for_asm
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating summary plot...', session=session,
                              duration=NULL, id='creating_summary')
      result <- generate_summary_plot(input$summary_txt, input$main_directory)
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        shiny::updateTextInput(session=session, 'summary_plot',
                               value=result$out_path)
        output$summary_pie_image <- shiny::renderImage(
          base::list(src=result$out_path, width='1000', height='800',
               alt='summary_pie_plot'),
          deleteFile=FALSE)
      }
      shiny::removeNotification(id='creating_summary', session=session)
    }),
    input$plot_summary_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating abundance plot...', session=session,
                              duration=NULL, id='creating_abundance')
      ## TODO use input$abundance_plot_file_type
      plot_file_type <- '.png'
      if (input$abundance != '') {
        result <- generate_abundance_and_structure_plots_from_abundance(
            input$gene_id, input$main_directory, input$abundance, input$updated_gtf,
            input$gencode_gtf, input$diff_isos, input$group_1,
            input$group_2, input$group_1_name, input$group_2_name,
            plot_file_type, input$abundance_intron_scaling,
            input$abundance_max_transcripts)
      } else {
        result <- generate_abundance_and_structure_plots(
            input$gene_id, input$asm_id, input$main_directory, input$event_dir,
            input$sample_total_tsv, input$count_tsv, input$diff_isos, input$group_1,
            input$group_2, input$group_1_name, input$group_2_name,
            plot_file_type, input$abundance_intron_scaling,
            input$abundance_max_transcripts)
      }
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        shiny::updateTextInput(session=session, 'isoform_colors',
                               value=result$colors_out_path)
        shiny::updateTextInput(session=session, 'asm_gtf',
                               value=result$asm_gtf_out_path)
        shiny::updateTextInput(session=session, 'modified_diff_isos',
                               value=result$modified_diff_isos)
        output$abundance_image <- shiny::renderImage(
          base::list(src=result$abun_plot_path, width='1000', height='700',
               alt='abundance_plot'),
          deleteFile=FALSE)
        output$structure_image <- shiny::renderImage(
          base::list(src=result$struct_plot_path, width='1000', height='333',
               alt='structure_plot'),
          deleteFile=FALSE)
        output$structure_sj_image <- shiny::renderImage(
          base::list(src=result$struct_sj_plot_path, width='1000', height='333',
               alt='structure_sj_plot'),
          deleteFile=FALSE)
      }
      shiny::removeNotification(id='creating_abundance', session=session)
    }),
    input$plot_abundance_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating in-gene plot...', session=session,
                              duration=NULL, id='creating_in_gene')
      if (base::file.exists(input$combined_gtf)) {
        in_gene_gtf <- input$combined_gtf
      } else {
        in_gene_gtf <- input$gencode_gtf
      }
      ## TODO use input$in_gene_plot_file_type
      plot_file_type <- '.png'
      result <- generate_in_gene_plot(
          input$gene_id, input$asm_id, input$main_directory, input$event_dir,
          input$modified_diff_isos, input$asm_gtf, in_gene_gtf,
          plot_file_type, input$in_gene_intron_scaling,
          input$in_gene_max_transcripts)
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        output$in_gene_image <- shiny::renderImage(
          base::list(src=result$in_gene_plot_path, width='1000', height='333',
                     alt='in_gene_plot'),
          deleteFile=FALSE)
        output$in_gene_sj_image <- shiny::renderImage(
          base::list(src=result$in_gene_sj_plot_path, width='1000', height='333',
                     alt='in_gene_sj_plot'),
          deleteFile=FALSE)
      }
      shiny::removeNotification(id='creating_in_gene', session=session)
    }),
    input$plot_in_gene_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating splicing differences table...',
                              session=session, duration=NULL,
                              id='creating_diffs')
      result <- generate_diff_table(input$gene_id, input$asm_id,
                                    input$main_directory, input$asm_gtf,
                                    input$isoform_colors)
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        output$isoform_differences_table <- shiny::renderTable(
            get_table_from_file(result$out_path))
      }
      shiny::removeNotification(id='creating_diffs', session=session)
    }),
    input$create_diff_table_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating simple splice graph plot...',
                              session=session, duration=NULL,
                              id='creating_simple_graph')
      result <- generate_simple_graph_plot(
          input$gene_id, input$asm_id, input$main_directory, input$asm_gtf,
          input$isoform_colors, input$simple_intron_scaling, input$equal_spacing,
          input$color_first)
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        output$simple_graph_image <- shiny::renderImage(
          base::list(src=result$out_path, width='1000', height='500',
                     alt='splice_graph_plot'),
          deleteFile=FALSE)
      }
      shiny::removeNotification(id='creating_simple_graph', session=session)
    }),
    input$plot_simple_graph_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('Creating splice graph plot...', session=session,
                              duration=NULL, id='creating_graph')
      result <- generate_graph_plot(input$gene_id, input$asm_id,
                                    input$main_directory, input$event_dir)
      if (result$error != FALSE) {
        shiny::showNotification(
          result$error, session=session, duration=NOTIFICATION_SECONDS,
          type='error')
      } else {
        output$graph_image <- shiny::renderImage(
          base::list(src=result$out_path, width='1000', height='333',
                     alt='splice_graph_plot'),
          deleteFile=FALSE)
      }
      shiny::removeNotification(id='creating_graph', session=session)
    }),
    input$plot_graph_button
  )
}
