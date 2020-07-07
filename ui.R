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
             <p>This web application shows features of UK lineages constructed using <a href='https://github.com/emvolz-phylodynamics/sarscov2Rutils' target='_blank'>sarscov2Rutils</a>. 
             Graphs show the time-dependent effective reproduction number (the average number of secondary cases per primary case over time), the effective population size (the number of individuals weighted by their contribution to producing the next generation) and the effective growth rate. </p>
             <p>Lineages are currently annotated as D/G from
             <a href='https://microreact.org/project/cogconsortium-2020-06-19/da7c0270/' target='_blank'> 
https://microreact.org/project/cogconsortium-2020-06-19/da7c0270/</a>.</p>
             <p>Trees show phylogenies contructed by <a href='https://cran.r-project.org/web/packages/treedater/index.html' target='_blank'>treedater</a>. Sequences lists the IDs and dates of the tree tips.</p>"
           )
    ),
    
    column(12, id = "plot",
           tabsetPanel(
             tabPanel("Trajectories",
                      column(3, id = "menu",
                             #div(id = "status",
                             #     shiny::verbatimTextOutput("output"),
                             #    shiny::htmlOutput("progress")
                             #),
                             div(id = "titles",
                                 div(id = "incidence_title", h4("Trajectory options")),
                                 shinyjs::hidden(div(id = "si_title",
                                                     h1("Serial Interval (SI) Input")))
                             ),
                             div(id = "1.1",
                                 div(id = "incidence_data_type_error_box", class = "error_box",
                                     checkboxInput('ti_ci', label = "Credible intervals", 1),
                                     checkboxInput(inputId='ti_log_size', label='Effective size on log scale', value=F)
                                     , checkboxGroupInput(inputId='ti_DGX', label='Genotypes', choices = NULL, selected = NULL)
                                     , checkboxGroupInput(inputId='ti_groups', label='', choices = NULL, selected = NULL)
                                     , checkboxGroupInput(inputId='ti_filename', label='Lineages', choices = NULL, selected = NULL)
                                     
                                 )
                             )#,
                             #div(id = "control",
                             #     shinyjs::hidden(actionButton("stop", label = "Stop")),
                             #     shinyjs::disabled(actionButton("prev", label = "Previous")),
                             #     shinyjs::disabled(actionButton("nxt", label = "Next")),
                             #     shinyjs::hidden(actionButton("go", label = "Go")),
                             #     textOutput("error")
                             #)
                             
                      ),
                      #         downloadButton("save_plot", "Save Image"),
                      column(6, fluidRow( plotOutput( 'GR', width = "100%", height = "300px"))
                             , fluidRow( plotOutput( 'R', width = "100%", height = "300px" ) )
                             , fluidRow( plotOutput( 'Ne', width = "100%", height = "300px" ) )),
                      column(3,plotOutput('legend'))
             ),
             tabPanel("Tree",
                      #downloadButton("save_incidence",
                      #               "Save Table"),
                      selectInput('ti_filename_tree', label = 'Lineage', NULL),
                      selectInput('ti_tree_colour', label = 'Colour by', choices=c('Genotype','Country','Location'),selected='Genotype'),
                      #fluidRow( plotlyOutput( 'tree', width = "100%", height = "auto"), align="right")
                      fluidRow( plotOutput( 'tree', width = "100%", height = "auto"), align="right")
             ),
             tabPanel("Sequences",
                      #downloadButton("save_r", "Save Table"),
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
