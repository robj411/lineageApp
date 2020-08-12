shinyUI(
  bootstrapPage(
    theme = shinythemes::shinytheme("spacelab"),
    tags$head(
      tags$link(rel = "stylesheet", type = "text/css", href = "styles.css"),
      tags$link(rel = "shortcut icon", type = "image/x-icon",
                href = "favicon.ico")
    ),
    # include the js code
    shinyjs::useShinyjs(),
    includeScript("scripts.js"),
    column(12,
           HTML(
             "
             <meta name='keywords' content='infectious,disease,epidemiology,genome,coronavirus,SARS-CoV-2,phylogenetics,lineage,mutation,Spike,UK'>
             <meta name='author' content='Rob Johnson'>
             <h1>COG-UK Phylodynamics dashboard: UK Lineages and Spike 614</h1>
             "
           ),
      tabsetPanel(
        tabPanel("About",
          HTML(
             "
             <p>Samples from lineages are summarised over time and geography in the UK according to the 
             details of the constituent sequences. Scalable phylodynamic methods are used to estimate 
             time-scaled phylogenies and estimate effective population size through time.  
             </p>
             <p>Lineage definitions are based on the dynamic nomenclature by  
             <a href='https://cov-lineages.org/' target='_blank'>Rambaut et al. 2020</a>.</p>
             <p>All data were compiled by the 
             <a href='https://www.cogconsortium.uk/' target='_blank'>UK COVID-19 Sequencing consortium</a>.</p>
             "
             )
        ),
        tabPanel("Methods",
                 HTML(
                   "
             <p>Skygrowth curves show aspects of UK lineages constructed using <a href='https://github.com/emvolz-phylodynamics/sarscov2Rutils' target='_blank'>sarscov2Rutils</a>. 
             Graphs show the time-dependent effective reproduction number (the average number of secondary cases per primary case over time), the effective population size (the number 
             of individuals weighted by their contribution to producing the next generation) and the effective growth rate (the growth rate of the effective population size). </p>
             <p>Lineage skygrowth curves are annotated in terms of their S (spike) protein variant as D/G in position 614.</p>
             <p>Region skygrowth curves are constructed as weighted, smoothed sums of the lineages' curves, as described 
             <a href='https://github.com/robj411/lineageApp/' target='_blank'>here</a>.</p>
             <p>Phylogenetic trees are constructed by <a href='https://cran.r-project.org/web/packages/treedater/index.html' target='_blank'>treedater</a>.</p>
             <p>Tabulated sequences are annotated as S-variant D/G/X and by location of collection from
             <a href='www.cogconsortium.uk' target='_blank'>COG-UK</a>.</p>
             "
                 )
        ),
        tabPanel("References",
                 HTML(
                   "
                   <p>A Rambaut, EC Holmes, Á O’Toole et al. 
             <a href='https://doi.org/10.1038/s41564-020-0770-5' target='_blank'>A dynamic nomenclature proposal for SARS-CoV-2 lineages to assist genomic epidemiology.</a>
              Nat Microbiol (2020).</p>
             <p>EM Volz, SDW Frost 
             <a href='https://doi.org/10.1093/ve/vex025' target='_blank'>Scalable relaxed clock phylogenetic dating.</a>
              Virus Evolution, 3:2 (2017).</p>
              <p>EM Volz, X Didelot
             <a href='https://doi.org/10.1093/sysbio/syy007' target='_blank'> Modeling the growth and decline of pathogen effective population size provides insight into epidemic dynamics and drivers of antimicrobial resistance</a>
              Systematic Biology (2018).</p>
              <p>EM Volz, V Hill, JT McCrone et al. 
             <a href='link' target='_blank'>Evaluating the effects of SARS-CoV-2 Spike mutation D614G on transmissibility and pathogenicity</a>
               (2020).</p>
             "
                 )
        )
        ),
    tags$hr(style="background: #cccccc; height: 4px;")
    ),
    column(12, id = "plot",
           tabsetPanel(
             tabPanel("Regions and lineages",
                      tabsetPanel(
                        tabPanel("Regions by lineage",
                                 column(2, selectInput('ti_filename_region', label = 'Lineage',choices='UK5 (G)', selected='UK5 (G)'),align='center'),
                                 column(5, plotOutput( 'byregion', width = "90%", height = "auto"),align='center'),
                                 column(5, plotOutput( 'hist_by_location3', width = "100%", height = "400px")),
                                 column(12, plotOutput( 'heatmap_reg_date_lin', width = "1200px", height = "400px"))
                        ),
                        tabPanel("Lineages by region",
                                 column(5,fluidRow(
                                        selectInput('ti_region', label = 'Region', NULL),
                                        sliderInput("ti_window", label = "Time window for lineage samples",
                                                    min = as.Date("2020-02-03","%Y-%m-%d"),
                                                    max = as.Date(lubridate::today(),"%Y-%m-%d"),
                                                    value=c(as.Date("2020-02-03"),as.Date(lubridate::today())),
                                                    timeFormat="%Y-%m-%d"),align='center'
                                 )
                                 ,
                                 fluidRow(plotOutput( 'linbyregion', width = "90%", height = "auto"), align="left"),
                                 h4("Samples over time"),
                                 fluidRow(plotOutput( 'hist_by_location2', width = "100%", height = "400px"), align="left")
                                 )
                                 ,
                                 column(5,
                                        h4("Combined skygrowth curves"),
                                        fluidRow( plotOutput( 'GR_reg', width = "100%", height = "300px"))
                                   , fluidRow( plotOutput( 'R_reg', width = "100%", height = "300px" ) )
                                   , fluidRow( plotOutput( 'Ne_reg', width = "100%", height = "300px" ) )
                                 )
                                 
                        ),
                        tabPanel("Heatplot",
                                 fluidRow(selectInput('ti_heat_norm', label = 'Normalisation', c('None','By lineage','By region'),selected='None')),
                                 column(12, plotOutput( 'heatmap', width = "1000px", height = "400px"))
                                 #column(6, plotOutput( 'hist_by_location3', width = "100%", height = "400px"))
                        )
                      )
             ),
             tabPanel("Samples over time",
                      column(4, id = "menu",
                             selectInput('ti_hist', label = '', choices=c('Combined','Country','Region','County','Local authority'),selected='Combined'),
                             checkboxGroupInput('ti_hist_level2', label = '', NULL)
                             #checkboxGroupInput(inputId='ti_filename', label='Lineages', choices = NULL, selected = NULL)
                      ),
                      column(8, 
                             fluidRow( plotOutput( 'hist_by_location', width = "100%", height = "600px"), align="right")
                      )
             ),
             tabPanel("Lineages summarised",
                      downloadButton("save_lineage_table", "Save table"),
                      #tableOutput("estimated_r_output")
                      fluidRow(
                        #checkboxInput(inputId='show_all_sequences', label='Show all sequences', value=F),
                        DT::dataTableOutput(outputId = 'lineages')
                        , style = "font-size:15px", align="left")
             ),
             tabPanel("Skygrowth curves",
                      column(3, id = "menu",
                             div(id = "1.2",
                                 div(id = "incidence_data_type_error_box", class = "error_box",
                                     selectInput('ti_curve', label = '', c('By lineage','By region'),selected='By lineage'),
                                     checkboxInput('ti_ci', label = "Credible intervals", 1),
                                     checkboxInput(inputId='ti_log_size', label='Effective size on log scale', value=F)
                                     , checkboxGroupInput(inputId='ti_DGX', label='Genotypes', choices = NULL, selected = NULL)
                                     , checkboxGroupInput(inputId='ti_groups', label='', choices = NULL, selected = NULL)
                                     , checkboxGroupInput(inputId='ti_filename', label='Lineages', choices = NULL, selected = NULL)
                                     
                                 )
                             )
                             
                      ),
                      #         downloadButton("save_plot", "Save Image"),
                      column(6, fluidRow( plotOutput( 'GR', width = "100%", height = "300px"))
                             , fluidRow( plotOutput( 'R', width = "100%", height = "300px" ) )
                             , fluidRow( plotOutput( 'Ne', width = "100%", height = "300px" ) )),
                      column(3,plotOutput('legend'))
             ),
             tabPanel("Phylogenetic trees",
                      #downloadButton("save_incidence",
                      #               "Save Table"),
                      selectInput('ti_filename_tree', label = 'Lineage', NULL),
                      selectInput('ti_tree_colour', label = 'Colour by', choices=c('Genotype','Country','Location'),selected='Genotype'),
                      HTML(
                        "
             <p>Click on legend names to isolate tree tips.</p>"
                      ),
                      fluidRow( plotly::plotlyOutput( 'tree', width = "100%", height = "auto"), align="right")
                      #fluidRow( plotOutput( 'tree', width = "100%", height = "auto"), align="right")
             ),
             tabPanel("Sequence data",
                      downloadButton("save_table", "Save table"),
                      #tableOutput("estimated_r_output")
                      fluidRow(
                        checkboxInput(inputId='show_all_sequences', label='Show all sequences', value=F),
                        DT::dataTableOutput(outputId = 'sequences')
                        , style = "font-size:15px", align="left")
             )
           )
    )
  )
)
