library( lubridate )
library( shiny )
#library(sarscov2)
library(grid)
#library(ape)
library(DT)
library(scales)
library(ggtree)
library( ggplot2 )
library(plotly)
library(dplyr)

if(file.exists('datasets/lineageSetup.Rdata')){
  load('datasets/lineageSetup.Rdata')
}else{
  ## functions ##############################################################
  quick_annotated_treeplot <- function( td , annotation='d614g', maxdate = NULL){ #date_decimal( max(tr$sts))
    tr = td 
    len <- length(unique(tr$data[[annotation]]))
    if(annotation=='d614g') {
      cols <- c('navyblue','turquoise','darkorange')
      names(cols) <- c('D','G','X')
      leg.dir <- "vertical"
      shapes <- rep(19,3)
    }else{
      if(len<10){
        cols <- rainbow(len)
        names(cols) <- unique(tr$data[[annotation]])
        shapes <- rep(19,len)
      }else{
        pchs <- c(15:18,25,5,4)
        if(len<60)  pchs <- c(15:18,25,5)
        if(len<50)  pchs <- c(15:18,25)
        if(len<40)  pchs <- c(15:18)
        reps <- ceiling(len/length(pchs))
        cols <- rep(rainbow(reps),times=length(pchs))[1:len]
        names(cols) <- unique(tr$data[[annotation]])
        shapes <- rep(pchs,each=reps)[1:len]
        names(shapes) <- unique(tr$data[[annotation]])
      }
      leg.dir <- "horizontal"
    }
    colScale <- scale_colour_manual(name = annotation,values = cols)
    shapeScale <- scale_shape_manual(name = annotation,values = shapes) 
    class( tr ) = 'phylo'
    maxdate <- date_decimal( max(tr$sts))
    #tipdeme <-  grepl( tr$tip.label, pat = region ) 
    tipdata <- data.frame( 
      taxa = tr$tip.label, 
      anno =  tr$data[[annotation]],
      text = paste0(tr$data$d614g,' / ',tr$data$location,' / ',tr$data$Date)
    )
    tipdata$size <- .75
    #tipdata$size[ !tipdata$d614g ] <- 0
    #tipdata$d614g[ !tipdata$d614g ] <- NA
    btr = ggtree(tr, mrsd= maxdate, ladderize=TRUE, as.Date=T,aes(text=text))  + theme_tree2() 
    btr <- btr %<+% tipdata 
    metat <- btr$data 
    
    #btr = btr + geom_tippoint( aes(color = anno, pch=anno, size = size, text=text), na.rm=TRUE, show.legend=TRUE, size =1.5) 
    btr = btr + geom_point( data=na.omit(metat),aes(color = anno, pch=anno)) 
    
    btr = btr + colScale + shapeScale +
      theme(legend.position='top', 
                           legend.justification='left',
                           legend.title=element_text(size=14), 
                           legend.text=element_text(size=12),
                           #legend.direction=leg.dir,
                           axis.text=element_text(size=12)) +
      scale_x_date(date_labels = "%Y/%m/%d")

        ## add some space under legend otherwise it covers tree
    hgt <- length(tr$Ti)
    height = 3*hgt + 1000 + ceiling(len/7)*20
    yval <- 1.12+ceiling(len/7)*0.005 - hgt/13000 #ceiling(len/7)/10
    btrplotly <- ggplotly(btr,height=height, tooltip = "text") %>%
      layout(legend = list(orientation = "h",xanchor='center',
                           yanchor='top',
                           x=0.4,y=yval)) %>%
      config(displayModeBar = F)   %>%
         style(hoverinfo = "none", traces = 1:2)
    btrplotly
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
    #stopifnot( inherits( x, 'sarscov2skygrowth' ))
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
  
  
  hist_by_location <- function(geog,geog_levels){
    subseq <- parms$sequences
    label_column <- 1
    if(geog=='Local authority'){ 
      label_column <- 'lad'
      subseq <- subseq[subseq[['lad']]%in%geog_levels,]
    }else if(geog!='Combined'){
      label_column <- tolower(geog)
      subseq <- subseq[subseq[[tolower(geog)]]%in%geog_levels,]
    }
    plt <- ggplot(subseq, aes(x=as.Date(Date), fill=eval(parse(text=label_column)))) +
      geom_histogram(bins=30,alpha=0.5, position="identity", color="black")
    if(geog!='Combined') {
      plt <- plt + guides(fill=guide_legend(title=label_column)) + 
      theme(legend.text=element_text(size=14),legend.title=element_text(size=14))
    }else{
      plt <- plt + guides(fill=FALSE)
    }
    plt + xlab('Date') + ylab ('Count') +
      theme(axis.text=element_text(size=14),axis.title=element_text(size=14) ,panel.grid.major = element_blank(), panel.grid.minor = element_blank(),  panel.background = element_blank())
  }
  
  ## prepare data ##########################################################
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
    trees[[i]]$data <- metadata[match(sequences[[i]]$Sequence,metadata$Sequence),] 
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
  
  ## store items
  parms <- list()
  parms$sequences <- do.call(rbind,sequences)
  metadata <- readRDS('datasets/lineage_metadata.Rds')
  parms$sequences <- left_join(parms$sequences,metadata,by='Sequence')
  #parms$sequences$d614g <- metadata$d614g[match(parms$sequences$Sequence,metadata$Sequence)] 
  #parms$sequences$location <- metadata$location[match(parms$sequences$Sequence,metadata$Sequence)] 
  parms$tree <- trees
  parms$filename <- parms$filename_tree <- ids[1]
  parms$ids <- unname(ids)
  parms$lineages <- lineages
  parms$cols <- cols
  parms$p614 <- p614
  parms$labels <- paste0(ids,' (',p614,')')
  parms$geno_labs <- paste0('All "',unique(parms$p614),'"')
  parms$track_geno_groups <- rep(0,length(parms$geno_labs))
  
  ## create combined tree
  DGSlist <- list('D','G',c('D','G'))
  parms$all_combined <- lapply(1:length(DGSlist),function(dg){
    all_combined <- list()
    for(metric in c('growth','R','Ne')){
      
      blocks <- list()
      for(i in 1:length(lineages)){
        if(parms$p614[i]%in%DGSlist[[dg]])
          blocks[[length(blocks)+1]] <- lineages[[i]][[metric]]
      }
      combined <- do.call(rbind,blocks)
      daytimes <- combined$time
      combined$days <- sapply(daytimes,function(x)strsplit(as.character(x),' ')[[1]][1])
      uniquetimes <- unique(combined$days)
      combined$var = ((combined$pc97.5-combined$pc2.5)/2/qnorm(0.975))^2
      if(metric=='Ne'){
        combined$var <- (pmin(abs(combined$pc97.5-combined$pc50),abs(combined$pc2.5-combined$pc50),na.rm=T)/qnorm(0.975))^2
      }
      grouptraj <- as.data.frame(t(sapply(uniquetimes,function(j){
        temp <- subset(combined,days==j)
        mu <- mean(temp$pc50)
        sd <- sqrt(sum(temp$var)/nrow(temp)^2)
        if(metric=='Ne') {
          out <- c()
          temp <- subset(temp,!is.na(pc97.5))
          for(i in 1:1000) out[i] <- sum((rnorm(nrow(temp),mean=temp$pc50,sd=sqrt(temp$var))))
          c(j,quantile(out,c(0.025,0.5,0.975)))
        }else{
          c(j,qnorm(c(0.025,0.5,0.975),mean=mu,sd=sd))
        }
      })))
      colnames(grouptraj) <- colnames(blocks[[1]])
      for(i in 2:4) grouptraj[,i] <- as.numeric(as.character(grouptraj[,i]))
      grouptraj[,1] <- as.Date(grouptraj[,1])
      neworder <- sort.int(grouptraj[,1],decreasing = F,index.return = T)
      all_combined[[metric]] <- grouptraj[neworder$ix,]
    }
    all_combined
  })
  
  ## save ######################################################################################
  ## clear environment
  rm(cols,files,gtdsfiles,i,ids,lineages,maxind,p614,sequences,trees,x,metadata,Ds,Gs,DGSlist)
  save(list=ls(),file= 'datasets/lineageSetup.Rdata')
}


shiny::shinyServer(function(input, output, session) {
  ## initialise checkboxes
  combined_labs <- c('Combined "D"','Combined "G"','Combined "D" and "G"')
  updateCheckboxGroupInput(session, inputId='ti_DGX', label='Genotypes', choices = parms$geno_labs)
  updateCheckboxGroupInput(session, inputId='ti_groups', label='', choices = combined_labs)
  updateCheckboxGroupInput(session, inputId='ti_filename', label = 'Lineages', choices = parms$labels, selected = parms$labels[which(parms$ids=='UK5')])
  updateSelectInput(session, inputId='ti_filename_tree', label = 'Lineage', choices = parms$labels, selected = parms$labels[which(parms$ids=='UK5')])
  
  ## initialise reactive values
  re.filename <- reactiveVal( parms$filename )
  re.filename_tree <- reactiveVal( parms$filename_tree )
  re.tree_colour <- reactiveVal( 'Genotype' )
  re.hist <- reactiveVal( 'Combined' )
  re.ci <- reactiveVal( 1 )
  re.DGX <- reactiveVal( NULL )
  re.DGXgroup <- reactiveVal( NULL )
  re.log_size <- reactiveVal( F )
  
  ## unique with na.rm
  unique.narm <- function(x,column='lad'){
    nms <- sort(unique(x))
    unnms <- nms[!is.na(nms)]
    tallies <- sapply(unnms,function(x)sum(parms$sequences[[column]]==x,na.rm=T))
    paste0(unnms,' (',tallies,')')
  }
  
  ## observe input events
  ## which geography to hist
  observeEvent( input$ti_hist, {
    re.hist( input$ti_hist )
    geog <- input$ti_hist
    if(geog=='Combined')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = '', NULL)
    if(geog=='Country')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Country', choices = unique.narm(parms$sequences$country,'country'), selected = 'england')
    if(geog=='Region')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Region', choices = unique.narm(parms$sequences$region,'region'), selected = 'london')
    if(geog=='County')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'County', choices = unique.narm(parms$sequences$county,'county'), selected = 'cambridgeshire')
    if(geog=='Local authority')
      updateCheckboxGroupInput(session, inputId='ti_hist_level2', label = 'Local authority', choices = unique.narm(parms$sequences$lad), selected = 'city of london')
  })
  ## plot credible intervals?
  observeEvent( input$ti_ci, {
    re.ci( input$ti_ci )
  })
  ## log size?
  observeEvent( input$ti_log_size, {
    re.log_size( input$ti_log_size )
  })
  
  ## plot all D or all G?
  observeEvent( input$ti_groups, {
    re.DGXgroup( input$ti_groups )
  }, ignoreNULL = FALSE)
  
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
  
  ## how to colour the tree tips
  observeEvent( input$ti_tree_colour, {
    re.tree_colour( input$ti_tree_colour )
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
    groups <- input$ti_groups
    if(is.null(groups)){
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
  
  output$GR <- renderPlot({
    plot( prep_plot(metric='growth'))
  })
  
  output$R <- renderPlot({
    plot( prep_plot(metric='R'))
  })
  
  output$Ne <- renderPlot({
    suppressWarnings(plot( prep_plot(metric='Ne')))
  })
  
  output$legend <- renderPlot({
    par(mar=c(0,0,0,0)); plot.new()#plot(c(0,1),c(0,1)); 
    legend(x=0,y=1,col=unique(parms$cols),lwd=3,bty='n',
           legend=unique(parms$p614),cex=1.25,title='S614 variant')
  })
  
  ## plot hist for geography
  observe({
    geog <- input$ti_hist
    geog_levels <- NULL
    if(geog!='Combined')
      geog_levels <- sapply(input$ti_hist_level2,function(x)strsplit(as.character(x),' \\(')[[1]][1])
    output$hist_by_location <- renderPlot({
      hist_by_location(geog,geog_levels)
    })#, 
    #height = 3*hgt+500)
  })
  
  ## plot tree with height depending on number of sequences
  observe({
    anno <- 'd614g'
    if(re.tree_colour()=='Location') anno <- 'location'
    if(re.tree_colour()=='Country') anno <- 'country'
    fname <- re.filename_tree() 
    l_ind <- match(fname,parms$labels)
    #hgt <- length(parms$tree[[l_ind]]$Ti)
    #output$tree <- renderPlot({
    output$tree <- renderPlotly({
      quick_annotated_treeplot(parms$tree[[l_ind]],annotation = anno) #%>% layout(height = 3*hgt+500)
    })#, 
    #height = 3*hgt+500)
  })
  
  
  ## sequences: show all sequences button
  observe({
    show_all_sequences <- input$show_all_sequences
    tab <- parms$sequences
    if(!show_all_sequences){
      fname <- re.filename_tree() 
      l_ind <- match(fname,parms$labels)
      tab <- subset(tab,Lineage%in%parms$ids[l_ind])
    }
    output$sequences <- DT::renderDataTable({
      datatable(tab,options = list("pageLength" = 500))
    })
    ## sequences: download button
    output$save_table <- downloadHandler(
      filename = function() { 
        paste("dataset-", Sys.Date(), ".csv", sep="")
      },
      content = function(file) {
        write.csv(tab, file)
      })
  })
  
  output 
  
}
)