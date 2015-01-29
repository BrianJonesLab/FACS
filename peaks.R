library(raster)
require(flowCore)
require(maps)
require(ggplot2)

FACSC20141216_files <- dir('.','.fcs',full.names = T)

findPeaks <- function(file,silent=F,pimp_my_plot=FALSE,peak_max=6) {
  data <- read.FCS(file, transformation=FALSE,alter.names = T)
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])
  to_plot <- data_ <- to_plot #in order to make colMeans in simplify function below
  
  coordinates(data_) <- ~ x + y
  
  template <- raster(res=0.1,ext=extent(data_),crs=NULL)
  test <- rasterize(data_,template,data_@coords[,1],fun='count')
  
  test_ <- na.exclude(as.data.frame(cbind(xyFromCell(test,cell = 1:ncell(test)),test[])))
  lines <- test_[test_[,3] > quantile(test_[,3],0.995,na.rm = T),1:3]
    
  simplify <- function(lines) {
    aaa <- lapply(1:nrow(lines), function(x) {
      test <- fields::rdist(lines[x,],lines[-x,]) <= 0.2
      sort(c(rownames(lines[x,]),rownames(lines[-x,])[test]))
    })
    
    if(nrow(lines[sapply(aaa,length) != 1,]) == 0) return(lines)
    
    bbb <- lines[sapply(aaa,length) == 1,]
    
    #the problem of this part is that the scale is in log scale so colMeans dos not applies properly  
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
  
  cc<-cc[values,]
  if((nrow(cc)< peak_max)) warning('Maximum number of peaks reached, reduce the threshold value in order to find more peaks')
  
  
  #   ggplot(test_,aes(x,y,fill=V3)) + geom_raster()  + xlab('log(PI-H)') + ylab('log(DAPI-H)') + geom_vline(aes(xintercept=x),lines,linetype='longdash',colour='red') + geom_hline(aes(yintercept=y),lines,linetype='longdash',colour='red')
    if(pimp_my_plot){
    names(test_)[3] <- 'density'
    pimp<-ggplot(test_, aes(x=x, y=y)) + 
    geom_tile(aes(fill = density)) + 
    coord_equal()+
    labs(x='log(PI-H)',y='log(DAPI-H)')+
    geom_vline(xintercept=cc$x,linetype=2)+
    geom_hline(yintercept=cc$y,linetype=2)+
    geom_point(data = cc,aes(x,y),size=7.5)+
    theme_bw()
    print(pimp)
    } else {
      plot(test,xlab='log(PI-H)',ylab='log(DAPI-H)')
      sapply(cc[,1],function(x) abline(v=x,lty=2))
      sapply(cc[,2],function(x) abline(h=x,lty=2))
      #apply(cc,1,function(x) points(x=x[1],y=x[2],col='blue',cex=2))
      symbols(cc,circles=rep(0.25,nrow(cc)),inches=F,add=T)
    }
  #cat(c("Peaks' coordinates:",apply(as.matrix(cc),1,function(x) paste0(round(x,2),collapse = ', '))),sep = '\n')
  response <- cbind(cc,peak.value=values,buffer)
  colnames(response) <- c('log(PI-H)','log(DAPI-H)','peak','buffer_0.25u')
  if(silent) response else {
    cat(paste0('***** ',basename(file),' *****\n'))
    print(response) 
  }
}

findPeaks(FACSC20141216_files[[6]])

findPeaks(FACSC20141216_files[[2]],pimp_my_plot = T)
