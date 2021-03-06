library(raster)
require(flowCore)
require(maps)

FACSC20141216_files <- dir('.','.fcs',full.names = T)

findPeaks <- function(file,silent=F) {
  data <- read.FCS(file, transformation=FALSE,alter.names = T)
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])
  to_plot <- data_ <- log(to_plot) 
  
  coordinates(data_) <- ~ x + y
  
  template <- raster(res=0.1,ext=extent(data_),crs=NULL)
  test <- rasterize(data_,template,data_@coords[,1],fun='count')
  
  test_ <- na.exclude(as.data.frame(cbind(xyFromCell(test,cell = 1:ncell(test)),test[])))
  lines <- test_[test_[,3] > quantile(test_[,3],0.995,na.rm = T),1:2]
    
  simplify <- function(lines) {
    aaa <- lapply(1:nrow(lines), function(x) {
      test <- fields::rdist(lines[x,],lines[-x,]) <= 0.2
      sort(c(rownames(lines[x,]),rownames(lines[-x,])[test]))
    })
    
    if(nrow(lines[sapply(aaa,length) != 1,]) == 0) return(lines)
    
    bbb <- lines[sapply(aaa,length) == 1,]
    
    aaaa <- sapply(max(sapply(aaa,length)):2,function(x){
      step <- aaa[sapply(aaa,length) == x]
      if(length(step) == 0) return(NULL)
      apply(unique(do.call(rbind,step)),1,function(y){
        colMeans(lines[y,])
      })
    })
    if(class(aaaa) == 'list') {
      aaaaa <- t(do.call(cbind,aaaa))
    } else {
      aaaaa <- t(aaaa)
    }
    
    colnames(aaaaa) <- colnames(bbb) <- c('x','y')
    
    simplify(rbind(aaaaa,bbb))
  }
  
  cc <- simplify(lines)
  
  buffer <- extract(test,cc,buffer=0.25,fun=sum)
  values <- extract(test,cc,fun=mean)
  
  
      plot(test,xlab='log(PI-H)',ylab='log(DAPI-H)')
      sapply(cc[,1],function(x) abline(v=x,lty=2))
      sapply(cc[,2],function(x) abline(h=x,lty=2))
      #apply(cc,1,function(x) points(x=x[1],y=x[2],col='blue',cex=2))
      symbols(cc,circles=rep(0.25,nrow(cc)),inches=F,add=T)
  #cat(c("Peaks' coordinates:",apply(as.matrix(cc),1,function(x) paste0(round(x,2),collapse = ', '))),sep = '\n')
  response <- cbind(cc,peak.value=values,buffer)
  colnames(response) <- c('log(PI-H)','log(DAPI-H)','peak','buffer_0.25u')
  if(silent) response else {
    cat(paste0('***** ',basename(file),' *****\n'))
    print(response) 
  }
}


findPeaks(FACSC20141216_files[[1]])

