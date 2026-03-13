library('bslib')
library('shiny')

ui <- fluidPage(
  shiny::tags$head(
    shiny::tags$style(shiny::HTML("
      .shiny-image-output {
        border: 1px solid;
        background-color: grey;
      }"
    ))
  ),

  shiny::tags$h1('rMATS-long'),
  bslib::navset_tab(
    bslib::nav_panel(
      shiny::tags$h2('Data'),
      shiny::tags$br(),
      shiny::tags$p('Provide the path to a directory used as --out-dir by rmats_long.py'),
      shiny::textInput('main_directory', 'main-directory', value=''),
      shiny::actionButton('set_main_directory', 'set-main-directory'),
      shiny::tags$br(),
      shiny::tags$br(),
      shiny::tags$p('Most of the remaining values can be automatically filled in once the main directory is set.'),
      shiny::tags$p('The values can also be set manually.'),
      shiny::tags$p('The values will be used to run commands from other tabs.'),
      shiny::textInput('group_1_name', 'group-1-name', value='group-1'),
      shiny::textInput('group_2_name', 'group-2-name', value='group-2'),
      shiny::textInput('group_1', 'group_1.txt', value=''),
      shiny::textInput('group_2', 'group_2.txt', value=''),
      shiny::textInput('diff_asms', 'differential_asms.tsv', value=''),
      shiny::textInput('diff_isos', 'differential_isoforms.tsv', value=''),
      shiny::textInput('summary_plot', 'summary_plot.png', value=''),
      shiny::textInput('summary_txt', 'summary.txt', value=''),
      shiny::textInput('align_dir', 'alignments-dir', value=''),
      shiny::textInput('event_dir', 'events-dir', value=''),
      shiny::textInput('count_tsv', 'count.tsv', value=''),
      shiny::textInput('combined_gtf', 'combined.gtf', value=''),
      shiny::textInput('sample_total_tsv', 'sample_read_totals.tsv',
                       value=''),
      shiny::textInput('gencode_gtf', 'gencode.gtf', value=''),
      shiny::tags$br(),
      shiny::tags$p('When an abundance plot is created from the ASM tab, differential_isoforms_modified.tsv and asm.gtf will be created'),
      shiny::textInput('modified_diff_isos',
                       'differential_isoforms_modified.tsv', value=''),
      shiny::textInput('asm_gtf', 'asm.gtf', value=''),
      shiny::tags$br(),
      shiny::tags$p('Only if rmats_long.py was run with --abundance will abundance.esp and updated.gtf be needed'),
      shiny::textInput('abundance', 'abundance.esp', value=''),
      shiny::textInput('updated_gtf', 'updated.gtf', value=''),
      value='data_tab'
    ),
    bslib::nav_panel(
      shiny::tags$h2('Summary'),
      bslib::navset_tab(
       bslib::nav_panel(
         shiny::tags$h3('Significant ASMs by splicing event'),
         shiny::imageOutput('summary_pie_image', width=1002, height=802),
         shiny::actionButton('plot_summary_button', 'create plot'),
         value='summary_pie_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Workflow summary'),
         shiny::verbatimTextOutput('summary_text'),
         value='summary_text_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('ASM statistical test results'),
         shiny::numericInput('asm_max_rows', 'max-rows', 100,
                             min=1, step=1),
         shiny::selectInput('asm_sort_col_1', 'sort-column-1', base::c(''),
                            selectize=FALSE),
         shiny::checkboxInput('asm_sort_col_1_decreasing', 'sort-decreasing-1'),
         shiny::selectInput('asm_sort_col_2', 'sort-column-2', base::c(''),
                            selectize=FALSE),
         shiny::checkboxInput('asm_sort_col_2_decreasing', 'sort-decreasing-2'),
         shiny::actionButton('sort_asm_results', 'sort-table'),
         shiny::tableOutput('asm_stat_table'),
         value='summary_asms_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Isoform statistical test results'),
         shiny::numericInput('iso_max_rows', 'max-rows', 100,
                             min=1, step=1),
         shiny::selectInput('iso_sort_col_1', 'sort-column-1', base::c(''),
                            selectize=FALSE),
         shiny::checkboxInput('iso_sort_col_1_decreasing', 'sort-decreasing-1'),
         shiny::selectInput('iso_sort_col_2', 'sort-column-2', base::c(''),
                            selectize=FALSE),
         shiny::checkboxInput('iso_sort_col_2_decreasing', 'sort-decreasing-2'),
         shiny::actionButton('sort_iso_results', 'sort-table'),
         shiny::tableOutput('isoform_stat_table'),
         value='summary_isoforms_tab'
       ),
       id='summary_navset',
       selected='summary_pie_tab'
      ),
      value='summary_tab'
    ),
    bslib::nav_panel(
      shiny::tags$h2('ASM'),
      shiny::tags$br(),
      shiny::tags$p('If rmats_long.py was run with --abundance then asm-id does not need to be set.'),
      shiny::tags$p('In that case some plots will not be created: isoform structure with splice junction counts, region within gene, full splice graph.'),
      shiny::tags$br(),
      shiny::tags$p('IDs can be found in the statstical result sub-tabs on the Summary tab'),
      shiny::textInput('asm_id', 'asm-id', value=''),
      shiny::textInput('gene_id', 'gene-id', value=''),
      shiny::actionButton('find_gene_id_for_asm', 'find-gene-id-for-asm'),
      shiny::tags$br(),
      shiny::tags$br(),
      shiny::tags$p('When an abundance plot is created, colors.tsv will also be created'),
      shiny::textInput('isoform_colors', 'colors.tsv', value=''),
      bslib::navset_tab(
       bslib::nav_panel(
         shiny::tags$h3('Isoform abundance'),
         shiny::imageOutput('abundance_image', width=1002, height=702),
         ## TODO allow creating the pdf and provide a download link
         ## shiny::selectInput('abundance_plot_file_type', 'plot-file-type',
         ##                    base::c('.pdf', '.png', 'all'), selected='.png'),
         shiny::numericInput('abundance_intron_scaling', 'intron-scaling', 1,
                             min=1, step=1),
         shiny::numericInput('abundance_max_transcripts', 'max-transcripts', 7,
                             min=1, step=1),
         shiny::actionButton('plot_abundance_button', 'create plot'),
         value='asm_abundance_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Isoform structure'),
         shiny::imageOutput('structure_image', width=1002, height=335),
         shiny::imageOutput('structure_sj_image', width=1002, height=335),
         value='asm_structure_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Region within gene'),
         shiny::imageOutput('in_gene_image', width=1002, height=335),
         shiny::imageOutput('in_gene_sj_image', width=1002, height=335),
         ## TODO allow creating the pdf and provide a download link
         ## shiny::selectInput('in_gene_plot_file_type', 'plot-file-type',
         ##                    base::c('.pdf', '.png', 'all'), selected='.png'),
         shiny::numericInput('in_gene_intron_scaling', 'intron-scaling', 1,
                             min=1, step=1),
         shiny::numericInput('in_gene_max_transcripts', 'max-transcripts', 7,
                             min=1, step=1),
         shiny::actionButton('plot_in_gene_button', 'create plot'),
         value='asm_in_gene_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Splicing differences'),
         shiny::tableOutput('isoform_differences_table'),
         shiny::actionButton('create_diff_table_button', 'create table'),
         value='asm_differences_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Simple splice graph'),
         shiny::imageOutput('simple_graph_image', width=1002, height=502),
         shiny::numericInput('simple_intron_scaling', 'intron-scaling', 1,
                             min=1, step=1),
         shiny::checkboxInput('equal_spacing', 'equal-spacing'),
         shiny::checkboxInput('color_first', 'color-first-transcript-only'),
         shiny::actionButton('plot_simple_graph_button', 'create plot'),
         value='asm_simple_graph_tab'
       ),
       bslib::nav_panel(
         shiny::tags$h3('Splice graph'),
         shiny::imageOutput('graph_image', width=1002, height=335),
         shiny::tags$h3('This may take a long time for large or complex ASMs'),
         shiny::actionButton('plot_graph_button', 'create plot'),
         value='asm_graph_tab'
       ),
       id='asm_navset',
       selected='asm_abundance_tab'
      ),
      value='asm_tab'
    ),
    id='top_navset',
    selected='data_tab'
  )
)
