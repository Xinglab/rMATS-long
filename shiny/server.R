library('shiny')

## These values can be set to override default locations
BASE_PATH <- base::normalizePath('..')
SHINY_SCRIPT_PATH <- base::paste0(BASE_PATH, '/shiny')
RMATS_LONG_SCRIPT_PATH <- base::paste0(BASE_PATH, '/scripts')
DATA_DIR <- base::paste0(SHINY_SCRIPT_PATH, '/data')
OUT_DIR <- base::paste0(SHINY_SCRIPT_PATH, '/output')
FORCE_775 <- FALSE  ## workaround for a filesystem issue
PYTHON_PATH <- 'python'
## Setting these to a conda environment can allow using conda
## installed dependencies without activating the conda environment.
## Added to PATH when creating plots (/path/to/conda_env/bin)
ENV_PATH <- ''
## Added to LD_LIBRARY_PATH when creating plots (/path/to/conda_env/lib)
ENV_LD_LIBRARY_PATH <- ''
## Used for R_LIBS_{USER,SITE} when creating plots (/path/to/conda_env/lib/R/library)
R_LIBS_PATH <- ''

base::cat(base::paste0('Using files at ', BASE_PATH, '\n'))

get_dataset_names <- function() {
  names <- base::list.files(DATA_DIR)
  return(names)
}

get_output_dir <- function(datasets) {
  dataset <- datasets[[1]]
  time_str <- base::format(base::Sys.time(), '%Y%m%dT%H%M%S')
  out_dir <- base::paste0(OUT_DIR, '/', dataset, '_', time_str)
  if (FORCE_775) {
    base::dir.create(out_dir, recursive=TRUE, mode='0775')
    base::Sys.chmod(out_dir, mode='0775', use_umask=FALSE)
  } else {
    base::dir.create(out_dir, recursive=TRUE)
  }
  return(out_dir)
}

parse_sample_names <- function(abundance_path) {
  if (abundance_path == '') {
    return(base::c(''))
  }

  header <- base::readLines(abundance_path, n=1)
  columns <- base::strsplit(header, '\t')[[1]]
  num_columns <- base::length(columns)
  sample_names <- columns[4:num_columns]
  return(sample_names)
}

parse_group_sample_names <- function(group_path) {
  if (group_path == '') {
    return(base::c(''))
  }

  header <- base::readLines(group_path, n=1)
  sample_names <- base::strsplit(header, ',')[[1]]
  return(sample_names)
}

write_group_samples <- function(samples, path) {
  con <- base::file(path, open='w')
  comma_samples <- base::paste(samples, collapse=',')
  base::write(comma_samples, file=con)
  base::close(con)
}

format_sample_names <- function(sample_names) {
  comma_samples <- base::paste(sample_names, collapse=', ')
  return(base::paste0('sample names: ', comma_samples))
}

get_file_paths <- function(datasets) {
  dataset <- datasets[[1]]
  if (dataset == '') {
    return(base::list(abun='', updated='', group_1='', group_2='', ref='',
                      diff_trans='', diff_genes='', summary_txt='',
                      summary_png='', summary_pdf=''))
  }

  sub_dir <- base::paste0(DATA_DIR, '/', dataset)
  abun <- base::paste0(sub_dir, '/abundance.esp')
  updated <- base::paste0(sub_dir, '/updated.gtf')
  group_1 <- base::paste0(sub_dir, '/group_1.txt')
  group_2 <- base::paste0(sub_dir, '/group_2.txt')
  ref <- base::paste0(sub_dir, '/reference.gtf')
  diff_trans <- base::paste0(sub_dir, '/differential_transcripts.tsv')
  diff_genes <- base::paste0(sub_dir, '/differential_genes.tsv')
  summary_txt <- base::paste0(sub_dir, '/summary.txt')
  summary_png <- base::paste0(sub_dir, '/summary_plot.png')
  summary_pdf <- base::paste0(sub_dir, '/summary_plot.pdf')
  return(base::list(abun=abun, updated=updated, group_1=group_1,
                    group_2=group_2, ref=ref, diff_trans=diff_trans,
                    diff_genes=diff_genes, summary_txt=summary_txt,
                    summary_png=summary_png, summary_pdf=summary_pdf))
}

find_gene_id <- function(file_paths, gene_names, gene_ids) {
  gene_name <- gene_names[[1]]
  gene_id <- gene_ids[[1]]
  if (gene_name == '') {
    gene_name <- "''"
  }
  if (gene_id == '') {
    gene_id <- "''"
  }
  script <- base::paste0(SHINY_SCRIPT_PATH, '/find_gene_id_name.py')
  args <- base::c(script, gene_name, gene_id, file_paths$updated, file_paths$ref)
  out_lines <- base::system2(PYTHON_PATH, args, stdout=TRUE)
  if (out_lines[1] == '') {
    return(NULL)
  }
  return(out_lines[1])
}

read_diff_transcripts <- function(path, gene_id) {
  parsed <- utils::read.table(path, header=TRUE, sep='\t')
  for_gene <- parsed[parsed$gene_id == gene_id, ]
  return(for_gene)
}

get_default_differential_transcripts_table <- function() {
  return(base::list(gene_id=base::c('...'),
                    feature_id=base::c('...')))
}

get_default_isoform_differences_table <- function() {
  return(base::list(transcript1=base::c('...'),
                    transcript2=base::c('...')))
}

get_default_image <- function() {
  return (base::paste0(SHINY_SCRIPT_PATH, '/clear_1px.png'))
}

get_gene_gtf <- function(file_paths, gene_id, output_dir) {
  updated_out_path <- base::paste0(output_dir, '/updated.gtf')
  ref_out_path <- base::paste0(output_dir, '/ref.gtf')
  combined_out_path <- base::paste0(output_dir, '/combined.gtf')
  base::file.create(updated_out_path)
  base::file.create(ref_out_path)
  base::file.create(combined_out_path)
  if (FORCE_775) {
    base::Sys.chmod(updated_out_path, mode='0775', use_umask=FALSE)
    base::Sys.chmod(ref_out_path, mode='0775', use_umask=FALSE)
    base::Sys.chmod(combined_out_path, mode='0775', use_umask=FALSE)
  }

  args <- base::c(gene_id, file_paths$updated)
  return_code <- base::system2('grep', args, stdout=updated_out_path)
  if (return_code != 0) {
    return(NULL)
  }

  args <- base::c(gene_id, file_paths$ref)
  return_code <- base::system2('grep', args, stdout=ref_out_path)
  if (return_code != 0) {
    return(NULL)
  }

  base::file.copy(updated_out_path, combined_out_path, overwrite=TRUE)
  ref_lines <- base::readLines(ref_out_path, warn=FALSE)
  con <- base::file(combined_out_path, open='a')
  base::writeLines(ref_lines, con)
  base::close(con)

  return(base::list(updated=updated_out_path, ref=ref_out_path,
                    combined=combined_out_path))
}

classify_isoform_differences <- function(gene_gtfs, transcript_id, output_dir) {
  out_tsv <- base::paste0(output_dir, '/isoform_differences.tsv')
  script <- base::paste0(RMATS_LONG_SCRIPT_PATH, '/classify_isoform_differences.py')
  args <- base::c(script, '--main-transcript-id', transcript_id,
                  '--updated-gtf', gene_gtfs$updated, '--gencode-gtf',
                  gene_gtfs$ref, '--out-tsv', out_tsv)
  return_code <- base::system2(PYTHON_PATH, args)
  if (return_code != 0) {
    return(NULL)
  }

  if (FORCE_775) {
    base::Sys.chmod(out_tsv, mode='0775', use_umask=FALSE)
  }
  isoform_differences <- utils::read.table(out_tsv, header=TRUE, sep='\t')
  return(isoform_differences)
}

get_command_env <- function() {
  old_path <- base::Sys.getenv('PATH')
  old_ld_path <- base::Sys.getenv('LD_LIBRARY_PATH')
  if (ENV_PATH != '') {
    new_path <- base::paste0(ENV_PATH, ':', old_path)
  } else {
    new_path <- old_path
  }
  if (ENV_LD_LIBRARY_PATH != '') {
    new_ld_path <- base::paste0(ENV_LD_LIBRARY_PATH, ':', old_ld_path)
  } else {
    new_ld_path <- old_ld_path
  }

  command_env <- base::c(base::paste0('PATH=', new_path),
                         base::paste0('LD_LIBRARY_PATH=', new_ld_path))
  if (R_LIBS_PATH != '') {
      command_env <- base::append(command_env,
                                  base::paste0('R_LIBS_USER=', R_LIBS_PATH))
      command_env <- base::append(command_env,
                                  base::paste0('R_LIBS_SITE=', R_LIBS_PATH))
  }

  return(command_env)
}

generate_plots <- function(file_paths, gene_id, output_dir, intron_scaling,
                           group_1_samples, group_2_samples, group_1_name,
                           group_2_name) {
  group_1_path <- base::paste0(output_dir, '/group_1.txt')
  group_2_path <- base::paste0(output_dir, '/group_2.txt')
  write_group_samples(group_1_samples, group_1_path)
  write_group_samples(group_2_samples, group_2_path)
  if (FORCE_775) {
    base::Sys.chmod(group_1_path, mode='0775', use_umask=FALSE)
    base::Sys.chmod(group_2_path, mode='0775', use_umask=FALSE)
  }
  quoted_group_1 <- base::paste0("'", group_1_name, "'")
  quoted_group_2 <- base::paste0("'", group_2_name, "'")

  script <- base::paste0(RMATS_LONG_SCRIPT_PATH, '/visualize_isoforms.py')
  command_env <- get_command_env()
  plot_file_type <- 'all'
  args <- base::c(script, '--gene-id', gene_id, '--abundance',
    file_paths$abun, '--updated-gtf', file_paths$updated, '--gencode-gtf',
    file_paths$ref, '--diff-transcripts', file_paths$diff_trans, '--out-dir',
    output_dir, '--plot-file-type', plot_file_type, '--intron-scaling',
    intron_scaling, '--group-1', group_1_path, '--group-2', group_2_path,
    '--group-1-name', quoted_group_1, '--group-2-name', quoted_group_2)
  return_code <- base::system2(PYTHON_PATH, args, env=command_env)
  return(return_code != 0)
}

generate_summary_plot <- function(file_paths, output_dir) {
  if (!base::file.exists(file_paths$summary_txt)) {
    return(base::paste0('Need ', file_paths$summary_txt,
                        ' to generate summary plot'))
  }

  script <- base::paste0(RMATS_LONG_SCRIPT_PATH, '/visualize_summary.R')
  command_env <- get_command_env()
  args <- base::c(script, file_paths$summary_txt, output_dir)
  return_code <- base::system2('Rscript', args, env=command_env)
  if (return_code != 0) {
    return('command to generate summary failed')
  }

  return(FALSE)
}

server <- function(input, output, session) {
  datasets <- get_dataset_names()
  shiny::updateSelectInput(session, 'dataset_select', choices=datasets)

  output$sample_names <- shiny::renderText(
    format_sample_names(parse_sample_names(
      get_file_paths(input$dataset_select)$abun)))

  shiny::observe({
    file_paths <- get_file_paths(input$dataset_select)
    parsed_sample_names <- parse_sample_names(file_paths$abun)
    group_1_selected <- parse_group_sample_names(file_paths$group_1)
    group_2_selected <- parse_group_sample_names(file_paths$group_2)
    shiny::updateSelectInput(session, 'group_1_samples',
                             choices=parsed_sample_names,
                             selected=group_1_selected)
    shiny::updateSelectInput(session, 'group_2_samples',
                             choices=parsed_sample_names,
                             selected=group_2_selected)
    output$download_diff_transcripts <- shiny::downloadHandler(
      filename='differential_transcripts.tsv',
      content=function(file) {base::file.copy(file_paths$diff_trans, file)})
    output$download_diff_genes <- shiny::downloadHandler(
      filename='differential_genes.tsv',
      content=function(file) {base::file.copy(file_paths$diff_genes, file)})

    default_image_path <- get_default_image()
    output$abundance_image <- shiny::renderImage(
      list(src=default_image_path, width='800', height='600',
           alt='default image'),
      deleteFile=FALSE)
    output$structure_image <- shiny::renderImage(
        list(src=default_image_path, width='800', height='300',
             alt='structure plot'),
        deleteFile=FALSE)
    if (base::file.exists(file_paths$summary_png)
        && base::file.exists(file_paths$summary_pdf)) {
      output$summary_image <- shiny::renderImage(
        list(src=file_paths$summary_png, width='1000', height='800',
             alt='summary plot'),
        deleteFile=FALSE)
      output$summary_png <- shiny::downloadHandler(
        filename='summary_plot.png',
        content=function(file) {base::file.copy(file_paths$summary_png, file)})
      output$summary_pdf <- shiny::downloadHandler(
        filename='summary_plot.pdf',
        content=function(file) {base::file.copy(file_paths$summary_pdf, file)})
    } else {
      output$summary_image <- shiny::renderImage(
              list(src=default_image_path, width='1000', height='800',
                   alt='default image'),
              deleteFile=FALSE)
    }

    output$differential_transcripts <- shiny::renderTable(
      get_default_differential_transcripts_table())
    output$isoform_differences <- shiny::renderTable(
      get_default_isoform_differences_table())

    bslib::nav_select('top_navset', selected='summary_tab',
                      session=session)
    bslib::nav_select('gene_navset', selected='gene_parameters_tab',
                      session=session)
  })

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('creating summary plot...', session=session,
                              duration=NULL, id='creating_summary')
      file_paths <- get_file_paths(input$dataset_select)
      output_dir <- get_output_dir(input$dataset_select)
      error <- generate_summary_plot(file_paths, output_dir)
      if (error) {
        message <- base::paste0('error generating summary plot: ', error)
        shiny::showNotification(message, session=session,
                                duration=20, type='error')
        shiny::removeNotification(id='creating_summary', session=session)
        return()
      }

      summary_png_path <- base::paste0(output_dir, '/summary_plot.png')
      summary_pdf_path <- base::paste0(output_dir, '/summary_plot.pdf')
      if (FORCE_775) {
        base::Sys.chmod(summary_png_path, mode='0775', use_umask=FALSE)
        base::Sys.chmod(summary_pdf_path, mode='0775', use_umask=FALSE)
      }
      base::file.copy(summary_png_path, file_paths$summary_png, overwrite=TRUE)
      base::file.copy(summary_pdf_path, file_paths$summary_pdf, overwrite=TRUE)
      if (FORCE_775) {
        base::Sys.chmod(file_paths$summary_png, mode='0775', use_umask=FALSE)
        base::Sys.chmod(file_paths$summary_pdf, mode='0775', use_umask=FALSE)
      }

      output$summary_image <- shiny::renderImage(
        list(src=file_paths$summary_png, width='1000', height='800',
             alt='summary plot'),
        deleteFile=FALSE)
      output$summary_png <- shiny::downloadHandler(
        filename='summary_plot.png',
        content=function(file) {base::file.copy(file_paths$summary_png, file)})
      output$summary_pdf <- shiny::downloadHandler(
        filename='summary_plot.pdf',
        content=function(file) {base::file.copy(file_paths$summary_pdf, file)})

      shiny::removeNotification(id='creating_summary', session=session)
    }),
    input$plot_summary_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('creating plots...', session=session, duration=NULL,
                              id='creating_plots')
      file_paths <- get_file_paths(input$dataset_select)
      gene_id <- find_gene_id(file_paths, input$gene_name, input$gene_id)
      if (base::is.null(gene_id)) {
        shiny::showNotification('no valid gene found from name and ID fields',
                                session=session, duration=20, type='error')
        shiny::removeNotification(id='creating_plots', session=session)
        return()
      } else {
        shiny::showNotification(base::paste0('using gene ID: ', gene_id),
                                session=session, duration=20)
      }
      diff_transcripts <- read_diff_transcripts(file_paths$diff_trans, gene_id)
      output$differential_transcripts <- shiny::renderTable(diff_transcripts)

      output_dir <- get_output_dir(input$dataset_select)
      error <- generate_plots(file_paths, gene_id, output_dir, input$intron_scaling,
        input$group_1_samples, input$group_2_samples, input$group_1_name,
        input$group_2_name)
      if (error) {
        shiny::showNotification('error generating plots', session=session,
                                duration=20, type='error')
        shiny::removeNotification(id='creating_plots', session=session)
        return()
      }
      image_path_base <- base::paste0(output_dir, '/', gene_id, '_')
      abundance_path_base <- base::paste0(image_path_base, 'abundance')
      abundance_path_png <- base::paste0(abundance_path_base, '.png')
      abundance_path_pdf <- base::paste0(abundance_path_base, '.pdf')
      structure_path_base <- base::paste0(image_path_base, 'structure')
      structure_path_png <- base::paste0(structure_path_base, '.png')
      structure_path_pdf <- base::paste0(structure_path_base, '.pdf')
      if (FORCE_775) {
        base::Sys.chmod(abundance_path_png, mode='0775', use_umask=FALSE)
        base::Sys.chmod(abundance_path_pdf, mode='0775', use_umask=FALSE)
        base::Sys.chmod(structure_path_png, mode='0775', use_umask=FALSE)
        base::Sys.chmod(structure_path_pdf, mode='0775', use_umask=FALSE)
      }

      output$abundance_image <- shiny::renderImage(list(src=abundance_path_png,
                                                        width='800', height='600',
                                                        alt='abundance plot'),
                                                   deleteFile=FALSE)
      output$structure_image <- shiny::renderImage(list(src=structure_path_png,
                                                        width='800', height='300',
                                                        alt='structure plot'),
                                                   deleteFile=FALSE)

      output$abundance_png <- shiny::downloadHandler(
        filename='abundance.png',
        content=function(file) {base::file.copy(abundance_path_png, file)})
      output$abundance_pdf <- shiny::downloadHandler(
        filename='abundance.pdf',
        content=function(file) {base::file.copy(abundance_path_pdf, file)})
      output$structure_png <- shiny::downloadHandler(
        filename='structure.png',
        content=function(file) {base::file.copy(structure_path_png, file)})
      output$structure_pdf <- shiny::downloadHandler(
        filename='structure.pdf',
        content=function(file) {base::file.copy(structure_path_pdf, file)})

      bslib::nav_select('gene_navset', selected='gene_abundance_tab',
                        session=session)
      shiny::removeNotification(id='creating_plots', session=session)
    }),
    input$plot_button
  )

  shiny::bindEvent(shiny::observe({
      shiny::showNotification('classifying isoforms...', session=session,
                              duration=NULL, id='classifying_isoforms')
      file_paths <- get_file_paths(input$dataset_select)
      gene_id <- find_gene_id(file_paths, input$gene_name, input$gene_id)
      if (base::is.null(gene_id)) {
        shiny::showNotification('no valid gene found from name and ID fields',
                                session=session, duration=20, type='error')
        shiny::removeNotification(id='classifying_isoforms', session=session)
        return()
      } else {
        shiny::showNotification(base::paste0('using gene ID: ', gene_id),
                                session=session, duration=20)
      }

      output_dir <- get_output_dir(input$dataset_select)
      gene_gtfs <- get_gene_gtf(file_paths, gene_id, output_dir)
      if (base::is.null(gene_gtfs)) {
        shiny::showNotification('error reading gtf files', session=session,
                                duration=20, type='error')
        shiny::removeNotification(id='classifying_isoforms', session=session)
        return()
      }
      output$updated_gtf <- shiny::downloadHandler(
        filename='updated.gtf',
        content=function(file) {base::file.copy(gene_gtfs$combined, file)})

      isoform_differences <- classify_isoform_differences(
        gene_gtfs, input$transcript_id, output_dir)
      if (base::is.null(isoform_differences)) {
        shiny::showNotification('error classifying isoforms', session=session,
                                duration=20, type='error')
        shiny::removeNotification(id='classifying_isoforms', session=session)
        return()
      }
      output$isoform_differences <- shiny::renderTable(isoform_differences)

      shiny::removeNotification(id='classifying_isoforms', session=session)
    }),
    input$classify_button
  )
}
