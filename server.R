library( lubridate )
library( ggplot2 )
library( shiny )
#library(sarscov2)
library(grid)
#library(ape)
library(ggtree)
library(DT)
library(scales)

if(file.exists('datasets/lineageSetup.Rdata')){
  load('datasets/lineageSetup.Rdata')
}else{
  
  quick_d614g_treeplot <- function( td , maxdate = NULL){ #date_decimal( max(tr$sts))
    cols <- c('navyblue','turquoise','darkorange')
    names(cols) <- c('D','G','X')
    colScale <- scale_colour_manual(name = "g614d",values = cols)
    tr = td 
    class( tr ) = 'phylo'
    btr = ggtree(tr, mrsd= maxdate, ladderize=TRUE)  + theme_tree2() 
    #tipdeme <-  grepl( tr$tip.label, pat = region ) 
    tipdata <- data.frame( 
      taxa = tr$tip.label, 
      d614g =  tr$d614g
    )
    tipdata$size <- .75
    #tipdata$size[ !tipdata$d614g ] <- 0
    #tipdata$d614g[ !tipdata$d614g ] <- NA
    btr <- btr %<+% tipdata 
    btr = btr + geom_tippoint( aes(color = d614g, size = size), na.rm=TRUE, show.legend=TRUE, size =1.25) 
    #+ theme_tree2( legend.position = "none" )
    
    btr + colScale + theme(legend.position='top', 
                           legend.justification='left',
                           legend.direction='vertical',
                           legend.title=element_text(size=14), 
                           legend.text=element_text(size=12))#, 
                           #legend.title = element_blank(), 
                           #legend.key = element_blank()) #+ ggplot2::ggtitle( region )
  }
  
  
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
  
  
  gg_sarscov2skygrowth <- function(x, metric='growth', log_size=F,date_limits = c( as.Date( '2020-03-01'), NA ) ,ci,col,... )
  {
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
    if(log_size) pl <- pl  +  scale_y_continuous(limits=c(1e-3,NA), trans='log',label=scientific_format(digits=2)) 
    
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
  
  ## get file names
  files <- list.files('datasets/skygrowth1/')
  gtdsfiles <- files[sapply(files,function(x)grepl('\\.gtds',x))]
  ids <- sapply(gtdsfiles,function(txt)gsub('.rds','',gsub('skygrowth1.gtds-','',txt)))
  lineages <- readRDS('datasets/skygrowth1/skygrowth1.sgs.rds')
                
  ## store information
  trees <- list()
  sequences <- list()
  metadata <- readRDS('datasets/lineage_metadata.Rds')
  for(i in 1:length(gtdsfiles)){
    x <- readRDS(paste0('datasets/skygrowth1/',gtdsfiles[i]))
    maxind <- which.max(sapply(x,function(y)y[[7]]))
    trees[[i]] <- x[[maxind]]
    sequences[[i]] <- as.data.frame(cbind(t(sapply(x[[1]]$tip.label,function(y)strsplit(y,'\\|')[[1]][1:2])),ids[i])) # sapply(x[[1]]$tip.label,function(y)gsub('England/','',y))
    colnames(sequences[[i]]) <- c('Sequence','Date','Lineage')
    sequences[[i]][,2] <- as.character(as.Date(date_decimal(as.numeric(as.character(sequences[[i]][,2])))))
    rownames(sequences[[i]]) <- NULL
    trees[[i]]$d614g <- metadata$d614g[match(sequences[[i]]$Sequence,metadata$Sequence)] 
  }
  
  ## get genotypes
  ##!! use old data to colour lineages
  metadata <- readRDS('datasets/lineage_metadata.Rds')
  metadata <- read.csv('https://microreact.org/api/viewer/data?project=cogconsortium-2020-06-19')
  metadata <- subset(metadata,uk_lineage%in%ids)[,colnames(metadata)%in%c('uk_lineage','d614g')]
  #GDs <- sapply(ids,function(x){subx <- subset(metadata,uk_lineage==x);(unique(subx$d614g))})
  Ds <- sapply(ids,function(x){subx <- subset(metadata,uk_lineage==x);sum(subx$d614g=='G')<0.1*nrow(subx)})
  Gs <- sapply(ids,function(x){subx <- subset(metadata,uk_lineage==x);sum(subx$d614g=='D')<0.1*nrow(subx)})
  cols <- rep('darkorange',length=length(ids))
  cols[Ds] <- 'navyblue'
  cols[Gs] <- 'turquoise'
  p614 <- rep('D/G/X',length=length(cols))
  p614[cols=='navyblue'] <- 'D'
  p614[cols=='turquoise'] <- 'G'
  
  ## save items
  parms <- list()
  parms$sequences <- do.call(rbind,sequences)
  metadata <- readRDS('datasets/lineage_metadata.Rds')
  parms$sequences$d614g <- metadata$d614g[match(parms$sequences$Sequence,metadata$Sequence)] 
  parms$tree <- trees
  parms$filename <- parms$filename_tree <- ids[1]
  parms$ids <- unname(ids)
  parms$lineages <- lineages
  parms$cols <- cols
  parms$p614 <- p614
  parms$labels <- paste0(ids,' (',p614,')')
  parms$geno_labs <- paste0('All "',unique(parms$p614),'"')
  parms$track_geno_groups <- rep(0,length(parms$geno_labs))
  ## clear environment
  rm(cols,files,gtdsfiles,i,ids,lineages,maxind,p614,sequences,trees,x,metadata,Ds,Gs)
  save(list=ls(),file= 'datasets/lineageSetup.Rdata')
}


shiny::shinyServer(function(input, output, session) {
  ## initialise checkboxes
  updateCheckboxGroupInput(session, inputId='ti_DGX', label='Genotypes', choices = parms$geno_labs)
  updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = parms$labels, selected = parms$labels[which(parms$ids=='UK5')])
  updateSelectInput(session, inputId='ti_filename_tree', label = '', choices = parms$labels, selected = parms$labels[which(parms$ids=='UK5')])
  re.filename <- reactiveVal( parms$filename )
  re.filename_tree <- reactiveVal( parms$filename_tree )
  re.ci <- reactiveVal( 1 )
  re.DGX <- reactiveVal( NULL )
  re.log_size <- reactiveVal( F )
  
  ## observe input events
  ## plot credible intervals?
  observeEvent( input$ti_ci, {
    re.ci( input$ti_ci )
  })
  ## log size?
  observeEvent( input$ti_log_size, {
    re.log_size( input$ti_log_size )
  })
  
  ## plot all D or all G?
  observeEvent( input$ti_DGX, {
    re.DGX( input$ti_DGX )
    new_labs <- as.numeric(parms$geno_labs%in%input$ti_DGX)
    group_to_update <- which(new_labs!=parms$track_geno_groups)
    ## only if there was a change to selection
    if(length(group_to_update)>0){
      old_selected <- input$ti_filename
      labs_to_include_or_exclude <- parms$labels[parms$p614%in%unique(parms$p614)[group_to_update]]
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
  observeEvent( input$ti_filename, {
    re.filename( input$ti_filename )
    ## if just one trajectory, update tree
    if(length(input$ti_filename)==1)
      updateSelectInput(session,"ti_filename_tree",selected=input$ti_filename)
  })
  
  ## which filename is selected for tree?
  observeEvent( input$ti_filename_tree, {
    re.filename_tree( input$ti_filename_tree )
    updateCheckboxGroupInput(session,"ti_filename",selected=unique(c(input$ti_filename_tree,input$ti_filename)))
  })
  
  ## plot outputs
  update_lineage <- reactive({
    fname <- re.filename() 
    l_ind <- match(fname,parms$labels)
    lineage <- lapply(l_ind,function(x)parms$lineages[[x]])
    list(lineage,l_ind)
  })
  
  prep_plot <- function(metric='growth'){
    ci <- re.ci() 
    log_size <- metric=='Ne' & re.log_size()
    lin <- update_lineage()
    sgs <- lin[[1]]
    ids <- lin[[2]]
    req(inherits(sgs[[1]], "sarscov2skygrowth"))
    p0 <- gg_sarscov2skygrowth(sgs[[1]],metric=metric,log_size=log_size,date_limits=c( as.Date( '2020-02-01'), NA ) ,ci,parms$cols[lin[[2]][1]])
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
    suppressWarnings(plot( prep_plot(metric='Ne')))
  })
  
  ## plot tree with height depending on number of sequences
  observe({
    fname <- re.filename_tree() 
    l_ind <- match(fname,parms$labels)
    hgt <- length(parms$tree[[l_ind]]$Ti)
    output$tree <- renderPlot({
      quick_d614g_treeplot(parms$tree[[l_ind]])
    }, 
    height = 3*hgt+500)
  })
  
  observe({
    show_all_sequences <- input$show_all_sequences
    output$sequences <- renderDataTable({
      tab <- parms$sequences
      if(!show_all_sequences){
        fname <- re.filename_tree() 
        l_ind <- match(fname,parms$labels)
        tab <- subset(tab,Lineage%in%parms$ids[l_ind])
      }
      datatable(tab,options = list("pageLength" = 500))
    })
  })
  
  output$legend <- renderPlot({
    par(mar=c(0,0,0,0)); plot.new()#plot(c(0,1),c(0,1)); 
    legend(x=0,y=1,col=unique(parms$cols),lwd=3,bty='n',
           legend=unique(parms$p614),cex=1.25,title='aa 614')
  })
  
  output 
  
}
)