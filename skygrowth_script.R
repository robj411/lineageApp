'
v3:
- using june 19 r1 data
- deltrans clusters n > 40 
- fixing clock rate on rtt consensus
'

# library(devtools)
# install_github('https://github.com/emvolz/treedater')

path = 'datasets/uk_lineages/' # directory with tree files
mdfn = 'datasets/majora_trimmed.csv' # path to csv 

use_gtds = TRUE 
minsize = 50
maxsize = 4000 # will downsample to this value 

nthreads = 1
ncpu = 16

ntrees = 16
SGSTART <- as.Date('2020-01-15' )
MHSTEPS <- 5e5
RES = 60
CLOCKRATE <- 0.0008254
TAU0 <- 1e-5
library( ape )
library( lubridate )
library(skygrowth)
library( sarscov2 )
library( treedater )

trfns = list.files( patt='.newick', path = path, full=TRUE )
lineage_names <- regmatches( trfns, regexpr( trfns, pattern='del_trans_([0-9]+)' ))
lineage_names <- regmatches( trfns, regexpr( trfns, pattern='UK([0-9]+)' ))

tres = lapply( trfns, read.tree )

#fix duplicate sequence names, run ONCE 
md = read.csv(mdfn,  stringsAs=FALSE, header=TRUE)
md$sequence_name <- md$central_sample_id
tips <- unlist(sapply(tres,function(x)sapply(x$tip,function(y)strsplit(as.character(y),'/')[[1]][2])))
to_replace <- !md$sequence_name%in%tips
secondary <- sapply(md$secondary_identifier[to_replace],function(x)strsplit(as.character(x),'/')[[1]][3])
md$sequence_name[to_replace] <- secondary
tips[!tips%in%md$sequence_name]
md <- md [!is.na(md$sequence_name)& nchar( md$sequence_name) > 0 , ]
.md = do.call( rbind, lapply( split( md, md$sequence_name ), function(x) x[1,] ))
rownames(.md) <- .md$sequence_name 

md$sample_time <- decimal_date( as.Date( md$collection_date ))

tr2sts <- lapply( tres, function(tr){
  nm <- sapply(tr$tip,function(x)strsplit(as.character(x),'/')[[1]][2])
  setNames( md[ match(nm,md$sequence_name) , 'sample_time'] , tr$tip.label )
})

#fix seq id, remove missing sts
tres <- lapply( 1:length(tres), function(k){
  sts = tr2sts[[k]]
  tr = tres[[k]]
  todrop = names( sts ) [ is.na( sts ) ]
  if ( length( todrop ) > (Ntip(tr)-3))
    todrop <- todrop[1:(Ntip(tr)-3)]
  
  drop.tip(tr, todrop  )
})

names(tres ) = names(tr2sts) <- lineage_names 

ns = sapply( tres, Ntip )
print( sort( ns ))
keep <- (ns > minsize ) 
tres = tres[keep ]
last_sample_time <- max(unlist(tr2sts),na.rm=T)

tr2sts <- tr2sts [ names(tres ) ] 
lineage_names = names( tres ) 

#~ stop()

.cutie_treedater <- function (tr, ntres = 10, threads = 2, ncpu = 5, mrl = c(CLOCKRATE, CLOCKRATE+1E-5)) 
{
  sts <- sapply(strsplit(tr$tip.label, "\\|"), function(x) {
    as.numeric(tail(x, 2)[1])
  })
  names(sts) <- tr$tip.label
  trpl <- tr
  tr <- di2multi(tr, tol = 1e-05)
  tres <- lapply(1:ntres, function(i) {
    tr = unroot(multi2di(tr))
    tr$edge.length <- pmax(1/29000/5, tr$edge.length)
    tr
  })
  tds <- parallel::mclapply(tres, function(tr) {
    dater(unroot(tr)
          , sts[tr$tip.label]
          , s = 29000
          , meanRateLimits = mrl
          , ncpu = ncpu
          , searchRoot = ncpu  - 1
    )
  }, mc.cores = threads)
  tds
}


.lineage_skygrowth <- function(tr, sts)
{
  # labels 
  tr3 = tr
  tr3$tip.label <- paste(sep='|', tr3$tip.label, sts[tr3$tip.label], 'Il' )
  
  # time trees 
  st0 = Sys.time() 
  tds = .cutie_treedater(tr3, ntres = ntrees, threads = nthreads, ncpu = ncpu) 
  st1 <- Sys.time() 
  print( paste( 'treedater', st1 - st0 ))
  
  
  # smooth node times 
  if ( use_gtds ){
    a = capture.output({
      gtds = parallel::mclapply( tds, function(td) gibbs_jitter( td, returnTrees=2 )[[2]] , mc.cores = ncpu*nthreads )
    })
    sg0 = skygrowth1( gtds, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
  } else{
    sg0 = skygrowth1( tds, tau0 = TAU0, res = RES, ncpu = ncpu ,  tstart = decimal_date(SGSTART), mhsteps = MHSTEPS)
    gtds = tds 
  }
  
  # skygrowth 
  sg0$GRmat  <- sg0$GRmat[ , sample(1:ncol(sg0$GRmat), size = 50, replace=FALSE) ] # shrink it
  sg0$Rmat <- NULL 
  sg0$Ntip = Ntip( tr3 ) 
  
  list( sg = sg0, gtds = gtds, tds = tds  )
}

print( Sys.time() )
#~ st0 = system.time( { sg0 = .lineage_skygrowth(tres[[1]])  } ) #71 sec 

if(!(dir.exists( 'skygrowth3' ))) dir.create( 'skygrowth3' )

k=27
## test on 44 and 77 (k=27,46)
sgs <- lapply( 1:length(tres), function(k){
  tr = tres[[k]]
  sts <- tr2sts[[k]] 
  ln <- names(tres)[k]
  print(c(ln, Ntip( tr )))
  sg_filename <- paste0('skygrowth3/skygrowth3-sg-', ln, '.rds' )
  if(!file.exists(sg_filename)){
    # downsample as needed 
    if ( Ntip( tr ) > maxsize ){
      m = Ntip(tr) - maxsize 
      drop = sample( tr$tip.label, size = m, replace=FALSE)
      tr = drop.tip( tr,  drop )
      sts = sts[ tr$tip.label ]
    }
    
    st0 = Sys.time() 
    o = .lineage_skygrowth(tr, sts )
    sg = o$sg
    gtds = o$gtds
    tds = o$tds
    st1 = Sys.time() 
    print( paste( 'skygrowth', st1 - st0 ))
    
    saveRDS(gtds, file = paste0('skygrowth3/skygrowth3-gtds-', ln, '.rds' ) )
    saveRDS(tds, file = paste0('skygrowth3/skygrowth3-tds-', ln, '.rds' ) )
    saveRDS(sg, file = sg_filename )
    print( Sys.time() )
    print( paste( lineage_names[k], 'complete' ) )
  }else{
    sg <- readRDS(sg_filename)
  }
  sg
})

names(sgs) <- lineage_names
saveRDS(sgs, file=paste0('skygrowth3/skygrowth3-sgs', '.rds' ) )