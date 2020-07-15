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
             <meta name='keywords' content='SARS-CoV-2'>
             <meta name='author' content='RJ'>
             <h1>UK Lineages</h1>
             <p>Samples from lineages are summarised over time and geography in the UK according to the details of the constituent sequences.  
             </p>
             <p>Skygrowth curves show aspects of UK lineages constructed using <a href='https://github.com/emvolz-phylodynamics/sarscov2Rutils' target='_blank'>sarscov2Rutils</a>. 
             Graphs show the time-dependent effective reproduction number (the average number of secondary cases per primary case over time), the effective population size (the number of individuals weighted by their contribution to producing the next generation) and the effective growth rate. </p>
             <p>Lineages' curves are currently annotated in terms of their S (spike) protein variant as D/G in position 614 from
             <a href='https://microreact.org/project/cogconsortium-2020-06-19/da7c0270/' target='_blank'> 
https://microreact.org/project/cogconsortium-2020-06-19/da7c0270/</a>.</p>
             <p>Phylogenetic trees are contructed by <a href='https://cran.r-project.org/web/packages/treedater/index.html' target='_blank'>treedater</a>.</p>
             <p>Tabulated sequences are currently annotated as S-variant D/G/X and by location of collection from
             <a href='www.cogconsortium.uk' target='_blank'> 
CoG-UK</a>.</p>
             "
          )
    ),
    
    column(12, id = "plot",
           tabsetPanel(
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
             tabPanel("Regions and lineages",
                      tabsetPanel(
                        tabPanel("Heatplot",
                                 fluidRow(selectInput('ti_heat_norm', label = 'Normalisation', c('None','By lineage','By region'),selected='None')),
                                 column(12, plotOutput( 'heatmap', width = "1000px", height = "400px"))
                                 #column(6, plotOutput( 'hist_by_location3', width = "100%", height = "400px"))
                        ),
                        tabPanel("Lineages by region",
                                 column(2,
                                        selectInput('ti_region', label = 'Region', NULL),
                                        sliderInput("ti_window", label = "Time window",
                                                    min = as.Date("2020-02-03","%Y-%m-%d"),
                                                    max = as.Date(lubridate::today(),"%Y-%m-%d"),
                                                    value=c(as.Date("2020-02-03"),as.Date(lubridate::today())),
                                                    timeFormat="%Y-%m-%d"),align='center'
                                 )
                                 ,
                                 column(5, plotOutput( 'linbyregion', width = "90%", height = "auto"), align="center"),
                                 column(5, plotOutput( 'hist_by_location2', width = "100%", height = "400px"), align="left")
                        ),
                        tabPanel("Regions by lineage",
                                 column(2, selectInput('ti_filename_region', label = 'Lineage', NULL),align='center'),
                                 column(5, plotOutput( 'byregion', width = "90%", height = "auto"),align='center'),
                                 column(5, plotOutput( 'hist_by_location3', width = "100%", height = "400px")),
                                 column(12, plotOutput( 'heatmap_reg_date_lin', width = "1200px", height = "400px"))
                        )
                      )
             ),
             tabPanel("Skygrowth curves",
                      column(3, id = "menu",
                             div(id = "1.2",
                                 div(id = "incidence_data_type_error_box", class = "error_box",
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
