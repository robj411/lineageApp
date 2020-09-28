library( lubridate )
library( shiny )
library(shinyjs)
library(grid)
#library(ape)
library(DT)
library(scales)
library(ggtree)
library( ggplot2 )
library(plotly)
library(dplyr)
library(RColorBrewer)
library(plotrix)

#scale_fill_discrete <- function(...) {
#  scale_fill_brewer(..., palette="Accent")
#}

load('datasets/lineageSetup.Rdata')


shiny::shinyServer(function(input, output, session) {
  ## initialise checkboxes
  combined_labs <- c('Combined "D"','Combined "G"','Combined "D" and "G"')
  updateCheckboxGroupInput(session, inputId='ti_DGX', label='Genotypes', choices = parms$geno_labs)
  updateCheckboxGroupInput(session, inputId='ti_groups', label='', choices = combined_labs)
  updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  updateSelectInput(session, inputId='ti_filename_tree', label = 'Lineage', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  updateSelectInput(session, inputId='ti_filename_region', label = 'Lineage', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
  updateSelectInput(session, inputId='ti_region', label = 'Region', choices = unique(parms$sequences_geog$region), selected = 'london')
  
  ## initialise reactive values
  re.filename <- reactiveVal( parms$filename )
  re.filename_tree <- reactiveVal( parms$filename_tree )
  re.hist <- reactiveVal( 'Combined' )
  
  ## unique with na.rm
  unique.narm <- function(x,column='lad'){
    nms <- sort(unique(x))
    unnms <- nms[!is.na(nms)]
    tallies <- sapply(unnms,function(x)sum(parms$sequences_geog[[column]]==x,na.rm=T))
    paste0(unnms,' (',tallies,')')
  }
  
  ## observe input events
  ## geography ########################################################
  ## which geography to hist
  observeEvent( input$ti_hist, {
    re.hist( input$ti_hist )
    geog <- input$ti_hist
    if(geog=='Combined')
      hide('ti_hist_level2')
      #updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = '', NULL)
    if(geog=='Country'){
      show('ti_hist_level2')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Country', choices = unique.narm(parms$sequences_geog$country,'country'), selected = 'england')
    }
    if(geog=='Region'){
      show('ti_hist_level2')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Region', choices = unique.narm(parms$sequences_geog$region,'region'), selected = 'london')
    }
    if(geog=='County'){
      show('ti_hist_level2')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'County', choices = unique.narm(parms$sequences_geog$county,'county'), selected = 'cambridgeshire')
    }
    if(geog=='Local authority'){
      show('ti_hist_level2')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Local authority', choices = unique.narm(parms$sequences_geog$lad), selected = 'city of london')
    }
  })
  
  ## plot hist for geography
  observe({
    geog <- input$ti_hist
    geog_levels <- NULL
    if(geog!='Combined')
      geog_levels <- sapply(input$ti_hist_level2,function(x)strsplit(as.character(x),' \\(')[[1]][1])
    output$hist_by_location <- renderPlot({
      hist_by_location(parms$sequences_geog,geog,geog_levels)
    })#, 
    #height = 3*hgt+500)
  })
  
  ## skygrowth ########################################################
  ## plot lineage or region?
  observeEvent( input$ti_curve, {
    if(input$ti_curve=='By region'){
      hide('ti_DGX')
      hide('ti_groups')
      reg.options <- names(parms$region_lineage)
      reg.colors <- parms$cols[1:length(parms$region_lineage)+length(as.character(parms$lineage_table$labels))]
      reg.fun <- lapply(reg.options,function(o)
        tags$div(
          HTML(paste(tags$span('―',style = paste0('font-weight: bold; color: ', reg.colors[which(reg.options == o)],';')),o, sep = " "))
          , encoding = 'UTF-8'
          )
      )
      updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Regions', 
                               #choices = names(parms$region_lineage),#reg.fun(),#
                               choiceNames = reg.fun,#names(parms$region_lineage),#
                               choiceValues = names(parms$region_lineage), #reg.colors,#
                               selected = 'london')
      output$legend <- renderPlot({
        par(mar=c(0,0,0,0)); plot.new()#plot(c(0,1),c(0,1)); 
      })
    }else{
      show('ti_DGX')
      show('ti_groups')
      updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = as.character(parms$lineage_table$labels), selected = as.character(parms$lineage_table$labels)[which(as.character(parms$lineage_table$ids)=='UK5')])
      output$legend <- renderPlot({
        par(mar=c(0,0,0,0)); plot.new()#plot(c(0,1),c(0,1)); 
        legend(x=0,y=1,col=unique(parms$cols),lwd=3,bty='n',
               legend=unique(as.character(parms$lineage_table$p614)),cex=1.25,title='S614 variant')
      })
    }
  })
  ## plot all D or all G?
  observeEvent( input$ti_DGX, {
    new_labs <- as.numeric(parms$geno_labs%in%input$ti_DGX)
    group_to_update <- which(new_labs!=parms$track_geno_groups)
    ## only if there was a change to selection
    if(length(group_to_update)>0){
      old_selected <- input$ti_filename
      labs_to_include_or_exclude <- as.character(parms$lineage_table$labels)[as.character(parms$lineage_table$p614)%in%unique(as.character(parms$lineage_table$p614))[group_to_update]]
      ## either add a set or remove a set
      if(new_labs[group_to_update]==1){
        new_selected <- unique(c(old_selected,labs_to_include_or_exclude))
      }else{
        new_selected <- old_selected[!old_selected%in%labs_to_include_or_exclude]
        if(length(new_selected)==0)
          new_selected <- input$ti_filename_tree
      }
      parms$track_geno_groups <<- new_labs
      updateCheckboxGroupInput(session,"ti_filename",selected=new_selected)
    }
  }, ignoreNULL = FALSE)
  
  ## which filenames are selected for trajectories?
  observeEvent( {
    input$ti_filename}, {
    re.filename( input$ti_filename )
    ## if just one trajectory, update tree
    if(length(input$ti_filename)==1&input$ti_curve=='By lineage')
      updateSelectInput(session,"ti_filename_tree",selected=input$ti_filename)
  })
  
  ## plot outputs
  update_lineage <- function(fname){
    if(input$ti_curve=='By region'){
      l_ind <- match(fname,names(parms$region_lineage)) 
      lineage <- lapply(l_ind,function(x)parms$region_lineage[[x]])
      l_ind <- l_ind + length(as.character(parms$lineage_table$labels))
    }else{
      l_ind <- match(fname,as.character(parms$lineage_table$labels))
      lineage <- lapply(l_ind,function(x)parms$lineages[[x]])
    }
    list(lineage,l_ind)
  }
  
  prep_plot <- function(metric='growth',ci=1,log_size=F,lin,region_only=F){
    log_size <- metric=='Ne' & log_size
    groups <- input$ti_groups
    if(is.null(groups)|input$ti_curve=='By region'|region_only){
      #print(names(lin[[1]]))
      sgs <- lin[[1]]
      ids <- lin[[2]]
      #if(!inherits(sgs[[1]], "sarscov2skygrowth")) print(lin[[2]])
      p0 <- gg_sarscov2skygrowth(sgs[[1]],metric=metric,log_size=log_size,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,parms$cols[ids[1]])
      if(length(ids)>1){
        for(i in 2:length(ids)){
          p0 <- add_line(p0,sgs[[i]],metric=metric,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,parms$cols[ids[i]])
        }
      }
      return(p0)
    }else{
      indices <- match(groups,combined_labs)
      p0 <- gg_sarscov2skygrowth(parms$all_combined[[indices[1]]],
                                 metric=metric,log_size=log_size,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,
                                 c('navyblue','turquoise','darkorange')[indices[1]])
      if(length(indices)>1){
        for(i in 2:length(indices)){
          p0 <- add_line(p0,parms$all_combined[[indices[i]]],
                         metric=metric,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,
                         c('navyblue','turquoise','darkorange')[indices[i]])
        }
      }
      return(p0)
    }
  }
  
  ## plot
  observeEvent({
    input$ti_filename
    input$ti_ci
    input$ti_log_size
    },{
    ci <- input$ti_ci
    fname <- re.filename() 
    lin <- update_lineage(fname)
    if(!is.na(lin[[2]][1])){
      output$GR <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='growth',ci=ci,lin=lin)))
      })
      output$R <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='R',ci=ci,lin=lin)))
      })
      log_size <- input$ti_log_size
      output$Ne <- renderPlot({
        suppressWarnings(plot( prep_plot(metric='Ne',ci=ci,log_size=log_size,lin=lin)))
      })
    }
  })
  
  ## lineage by region ########################################################
  ## bar plot lineage by region
  observeEvent({
    input$ti_region
    input$ti_window
    },{
    region <- input$ti_region
    timerange <- input$ti_window
    output$linbyregion <- renderPlot({
      barplot_by_region(parms$sequences_geog,region,'Lineage',timerange)
    },height=100+10*length(unique(parms$sequences_geog$Lineage[parms$sequences_geog$region==region])))
    ## plot hist for geography
    output$hist_by_location2 <- renderPlot({
      hist_by_location(parms$sequences_geog,geog='region',geog_levels=region)
    })
    ## plot skygrowth curves
    ci <- 1
    lin <- list(list(parms$region_lineage[[region]]) ,1)
    output$GR_reg <- renderPlot({
      suppressWarnings(plot( prep_plot(metric='growth',ci=ci,lin=lin,region_only=T)))
    })
    output$R_reg <- renderPlot({
      suppressWarnings(plot( prep_plot(metric='R',ci=ci,lin=lin,region_only=T)))
    })
    output$Ne_reg <- renderPlot({
      suppressWarnings(plot( prep_plot(metric='Ne',ci=ci,log_size=F,lin=lin,region_only=T)))
    })
  })
  
  ## regions by lineage ########################################################
  ## bar plot lineage by region
  observe({
    lin <- sapply(input$ti_filename_region,function(x)strsplit(as.character(x),' \\(')[[1]][1])
    output$byregion <- renderPlot({
      barplot_by_region(parms$sequences_geog,lin,'region')
    },height=100+10*length(unique(parms$sequences_geog$region[parms$sequences_geog$Lineage==lin])))
    ## plot hist for geography
    output$hist_by_location3 <- renderPlot({
      hist_by_location(parms$sequences_geog,geog='Lineage',geog_levels=lin)
    })
    ## plot heatmap
    output$heatmap_reg_date_lin <- renderPlot({
      par(mar=c(7,12,1,5.5))
      tab <- parms$region_date_lineage_table[[lin]]
      collabs <- as.Date(colnames(tab))
      blueheat(tab,collabs)
    })
  })
  
  ## regions and lineages ########################################
  observe({
    output$heatmap <- renderPlot({
      par(mar=c(5,12,1,5.5))
      normalisation <- input$ti_heat_norm
      tab <- parms$region_table
      collabs <- colnames(parms$region_table)
      if(normalisation=='By lineage') for(i in 1:ncol(tab)) tab[,i] <- tab[,i]/sum(tab[,i])
      if(normalisation=='By region') for(i in 1:nrow(tab)) tab[i,] <- tab[i,]/sum(tab[i,])
      blueheat(tab,collabs)
    })
  })
#  observe({
#    output$heatmap <- renderPlot({
#      par(mar=c(12,12,5,5.5))
#      blueheat(t(tab),collabs,collegend='top')
#    })
#  })
  
  ## tree ########################################################
  ## which filename is selected for tree?
  observeEvent( input$ti_filename_tree, {
    re.filename_tree( input$ti_filename_tree )
    updateCheckboxGroupInput(session,"ti_filename",selected=unique(c(input$ti_filename_tree,input$ti_filename)))
  })
  
  ## plot tree with height depending on number of sequences
  observe({
    anno <- 'd614g'
    if(input$ti_tree_colour =='Location') anno <- 'location'
    if(input$ti_tree_colour =='Country') anno <- 'country'
    fname <- re.filename_tree() 
    l_ind <- match(fname,as.character(parms$lineage_table$labels))
    #hgt <- length(parms$tree[[l_ind]]$Ti)
    #output$tree <- renderPlot({
    if(fname!='')
    output$tree <- renderPlotly({
      quick_annotated_treeplot(parms$tree[[l_ind]],annotation = anno) #%>% layout(height = 3*hgt+500)
    })#, 
    #height = 3*hgt+500)
  })
  
  ## lineages #####################################################
  output$lineages <- DT::renderDataTable({
    datatable(parms$lineage_table[,colnames(parms$lineage_table)%in%c("Lineage","p614","Number.of.samples", "First.date","Most.recent.date")],options = list("pageLength" = 500),rownames=F)
  })
  
  ## lineages: download button
  output$save_lineage_table <- downloadHandler(
    filename = function() { 
      paste("lineage-dataset-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(parms$lineage_table[,colnames(parms$lineage_table)%in%c("Lineage","p614","Number.of.samples", "First.date","Most.recent.date")], file, row.names = F)
    })
  
  ## sequences: show all sequences button ##############################
  observe({
    show_all_sequences <- input$show_all_sequences
    tab <- parms$sequences[,colnames(parms$sequences)%in%c("Sequence","Lineage","gisaid_epi_isl","adm1")]
    if(!show_all_sequences){
      fname <- re.filename_tree() 
      l_ind <- match(fname,as.character(parms$lineage_table$labels))
      tab <- subset(tab,Lineage%in%as.character(parms$lineage_table$ids)[l_ind])
    }
    output$sequences <- DT::renderDataTable({
      datatable(tab,options = list("pageLength" = 500),rownames=F)
    })
    ## sequences: download button
    output$save_table <- downloadHandler(
      filename = function() { 
        paste("sequences-dataset-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(tab, file, row.names = F)
      })
  })
  
  output 
  
}
)
