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
  shiny::selectInput('dataset_select', 'select dataset', choices=base::c('')),
  bslib::navset_tab(
    bslib::nav_panel(
      'Summary',
      shiny::tags$h2('Differential transcript usage result files'),
      shiny::downloadButton('download_diff_transcripts', 'download differential_transcripts.tsv'),
      shiny::downloadButton('download_diff_genes', 'download differential_genes.tsv'),
      shiny::tags$h2('Isoform Switch Splicing Events'),
      shiny::imageOutput('summary_image', width=1002, height=802),
      shiny::downloadButton('summary_png', 'download .png'),
      shiny::downloadButton('summary_pdf', 'download .pdf'),
      shiny::actionButton('plot_summary_button', 'create summary plot'),
      value='summary_tab'
    ),
    bslib::nav_panel(
      'Gene',
      bslib::navset_tab(
        bslib::nav_panel(
          'Parameters',
          shiny::textOutput('sample_names'),
          shiny::selectInput('group_1_samples', 'select group 1 samples',
                             choices=base::c(''), multiple=TRUE),
          shiny::textInput('group_1_name', 'group 1 name', value='group 1'),
          shiny::selectInput('group_2_samples', 'select group 2 samples',
                             choices=base::c(''), multiple=TRUE),
          shiny::textInput('group_2_name', 'group 2 name', value='group 2'),
          shiny::textInput('gene_name', 'gene name', value=''),
          shiny::textInput('gene_id', 'gene ID', value=''),
          shiny::numericInput('intron_scaling', 'intron scaling', 1, min=1, step=1),
          shiny::actionButton('plot_button', 'create plots'),
          value='gene_parameters_tab'
        ),
        bslib::nav_panel(
          'Abundance',
          shiny::imageOutput('abundance_image', width=802, height=602),
          shiny::downloadButton('abundance_png', 'download .png'),
          shiny::downloadButton('abundance_pdf', 'download .pdf'),
          value='gene_abundance_tab'
        ),
        bslib::nav_panel(
          'Structure',
          shiny::imageOutput('structure_image', width=802, height=302),
          shiny::downloadButton('structure_png', 'download .png'),
          shiny::downloadButton('structure_pdf', 'download .pdf'),
          value='gene_structure_tab'
        ),
        bslib::nav_panel(
          'Differential transcripts',
          shiny::tableOutput('differential_transcripts'),
          value='gene_diff_transcripts_tab'
        ),
        id='gene_navset',
        selected='gene_parameters_tab'
      ),
      value='gene_tab'
    ),
    bslib::nav_panel(
      'Transcript',
      shiny::tags$h3('Isoform differences'),
      shiny::textInput('transcript_id', 'transcript ID', value=''),
      shiny::actionButton('classify_button', 'classify isoform differences'),
      shiny::tableOutput('isoform_differences'),
      shiny::downloadButton('updated_gtf', 'download detected isoform .gtf'),
      value='transcript_tab'
    ),
    id='top_navset',
    selected='summary_tab'
  )
)
