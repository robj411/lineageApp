library( deSolve ) 
library( lubridate )
library( ggplot2 )
library( shiny )
library(sarscov2)
library(grid)
library(ape)
library(ggtree)
library(DT)

if(file.exists('data/lineageSetup.Rdata')){
  load('data/lineageSetup.Rdata')
}else{
  add_line <- function(pl,x1,metric='Ne',  date_limits = c( as.Date( '2020-02-01'), NA ),ci ,col,... ){
    y = x1[[metric]]
    pldf <- data.frame( Date =  as.Date( y$time ) , reported=FALSE )
    pldf[[metric]] = y$pc50
    pl = pl  + geom_path( data = pldf, aes( x = Date , y = get(metric) ), lwd=0.8,col=col )
    if(ci==1) {
      pldf$`2.5%` = y$pc2.5
      pldf$`97.5%` = y$pc97.5
      pl <- pl + geom_ribbon( data = pldf, aes(x = Date, ymin=`2.5%`, ymax=`97.5%`),col=col , alpha = .1, lwd=0)
    }
    
    pl
  }
  
  
  gg_sarscov2skygrowth <- function(x, metric='growth', date_limits = c( as.Date( '2020-03-01'), NA ) ,ci,col,... )
  {
    require(ggplot2)
    require(lubridate)
    stopifnot( inherits( x, 'sarscov2skygrowth' ))
    y = x[[metric]]
    taxis = as.Date( y$time )
    
    if ( is.na( date_limits[2]) )
      date_limits[2] <- as.Date( date_decimal( max(taxis)  ) )
    #qs <- c( .5, .025, .975 )
    
    pldf <- data.frame( Date = taxis , reported=FALSE )
    pldf$out = y$pc50
    pldf$`2.5%` = y$pc2.5
    pldf$`97.5%` = y$pc97.5
    
    pldf <- pldf[ with( pldf, Date > date_limits[1] & Date <= date_limits[2] ) , ]
    pl = ggplot( pldf) + 
      geom_path( aes(x = Date, y = out ), lwd=0.8,col=col) 
    
    if(ci==1) {
      pl <- pl +  geom_ribbon( aes(x = Date, ymin=`2.5%`, ymax=`97.5%`) , alpha = .1 ,col=col,lwd=0) 
    }
    
    y_lab <- 'Effective growth rate'
    if(metric=='R') y_lab <- 'Effective reproduction number'
    if(metric=='Ne') y_lab <- 'Effective population size'
    
    if(metric%in%c('R','growth')) pl <- pl + geom_hline( aes(yintercept = 1 ), colour = 'red' )
    pl <- pl + xlab('')  +theme(axis.text=element_text(size=12),axis.title=element_text(size=14) ,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                                panel.background = element_blank())+ 
      ylab ( y_lab) 
    pl
  }
  
  files <- list.files('data/skygrowth1/')
  sgfiles <- files[sapply(files,function(x)grepl('\\.sg-',x))]
  gtdsfiles <- files[sapply(files,function(x)grepl('\\.gtds',x))]
  ids <- sapply(sgfiles,function(txt)gsub('.rds','',gsub('skygrowth1.sg-','',txt)))
  lineages <- readRDS('data/skygrowth1/skygrowth1.sgs.rds')
  #lengths <- min_time <- max_time <- c()
  trees <- list()
  sequences <- list()
  for(i in 1:length(sgfiles)){
    #lengths[i] <- lineages[[ids[i]]]$Ntip
    x <- readRDS(paste0('data/skygrowth1/',gtdsfiles[i]))
    #times <- range(sapply(x,function(y)y[[14]]))
    #min_time[i] <- times[1]
    #max_time[i] <- times[2]
    maxind <- which.max(sapply(x,function(y)y[[7]]))
    trees[[i]] <- x[[maxind]]
    sequences[[i]] <- as.data.frame(cbind(t(sapply(x[[1]]$tip.label,function(y)strsplit(y,'\\|')[[1]][1:2])),ids[i])) # sapply(x[[1]]$tip.label,function(y)gsub('England/','',y))
    colnames(sequences[[i]]) <- c('Sequence','Date','Lineage')
    sequences[[i]][,2] <- as.character(as.Date(date_decimal(as.numeric(as.character(sequences[[i]][,2])))))
    rownames(sequences[[i]]) <- NULL
  }
  cols <- read.csv('data/DGcols.csv',stringsAsFactors=F,header = F)[,1]
  cols <- readRDS('data/DGcols.Rds')
  parms <- list()
  parms$sequences <- do.call(rbind,sequences)
  parms$tree <- trees
  parms$filename <- parms$filename_tree <- ids[1]
  parms$ids <- unname(ids)
  parms$lineages <- lineages
  parms$cols <- cols
  p614 <- readRDS('data/DGp614.Rds')
  parms$p614 <- p614
  labels <- paste0(ids,' (',p614,')')
  parms$labels <- labels
  parms$geno_labs <- paste0('All "',unique(parms$p614),'"')
  parms$track_geno_groups <- rep(0,length(parms$geno_labs))
  rm(cols,files,gtdsfiles,i,ids,labels,lineages,maxind,p614,sequences,sgfiles,trees,x)
  save(list=ls(),file= 'data/lineageSetup.Rdata')
}


shiny::shinyServer(function(input, output, session) {
  
  updateCheckboxGroupInput(session, inputId='ti_DGX', label='Genotypes', choices = parms$geno_labs)
  updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = parms$labels, selected = parms$labels[1])
  updateSelectInput(session, inputId='ti_filename_tree', label = '', parms$labels)
  re.filename <- reactiveVal( parms$filename )
  re.filename_tree <- reactiveVal( parms$filename_tree )
  re.ci <- reactiveVal( 1 )
  re.DGX <- reactiveVal( NULL )
  
  observeEvent( input$ti_ci, {
    re.ci( input$ti_ci )
  })
  
  observeEvent( input$ti_DGX, {
    re.DGX( input$ti_DGX )
    new_labs <- as.numeric(parms$geno_labs%in%input$ti_DGX)
    group_to_update <- which(new_labs!=parms$track_geno_groups)
    if(length(group_to_update)>0){
      old_selected <- input$ti_filename
      labs_to_include_or_exclude <- parms$labels[parms$p614%in%unique(parms$p614)[group_to_update]]
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
  
  observeEvent( input$ti_filename, {
    re.filename( input$ti_filename )
  })
  
  observeEvent( input$ti_filename_tree, {
    re.filename_tree( input$ti_filename_tree )
    updateCheckboxGroupInput(session,"ti_filename",selected=unique(c(input$ti_filename_tree,input$ti_filename)))
  })
  
  update_lineage <- reactive({
    fname <- re.filename() 
    l_ind <- match(fname,parms$labels)
    lineage <- lapply(l_ind,function(x)parms$lineages[[x]])
    list(lineage,l_ind)
  })
  
  prep_plot <- function(metric='growth'){
    ci <- re.ci() 
    lin <- update_lineage()
    sgs <- lin[[1]]
    ids <- lin[[2]]
    req(inherits(sgs[[1]], "sarscov2skygrowth"))
    p0 <- gg_sarscov2skygrowth(sgs[[1]],metric=metric,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,parms$cols[lin[[2]][1]])
    if(length(ids)>1){
      for(i in 2:length(ids)){
        p0 <- add_line(p0,sgs[[i]],metric=metric,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,parms$cols[lin[[2]][i]])
      }
    }
    p0
  }
  
  output$GR <- renderPlot({
    plot( prep_plot(metric='growth'))
  })
  
  output$R <- renderPlot({
    plot( prep_plot(metric='R'))
  })
  
  output$Ne <- renderPlot({
    plot( prep_plot(metric='Ne'))
  })
  
  output$tree <- renderPlot({
    fname <- re.filename_tree() 
    l_ind <- match(fname,parms$labels)
    #plot.phylo(parms$tree[[l_ind]],show.tip.label = F)
    quick_region_treeplot(parms$tree[[l_ind]],'England')
  })
  
  output$sequences <- renderDataTable({
    #fname <- re.filename_tree() 
    #l_ind <- match(fname,parms$labels)
    datatable(parms$sequences,options = list("pageLength" = 400))
    
  })
  
  output$legend <- renderPlot({
    par(mar=c(0,0,0,0)); plot.new()
    legend(x=0,y=1,col=unique(parms$cols),lwd=3,bty='n',
           legend=unique(parms$p614),cex=1.25,title='aa 614')
  })
  
  output 
  
}
)