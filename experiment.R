require(flowCore)
require(MASS)
require(ggplot2)
require(reshape2)
require(Rwave)
FACSC20141216_files <- dir('.','.fcs',full.names = T)

# findPeaks <- function(file,
#                       peaks=6,
#                       resol=500) {

  resol=500
  tolerance<-.995
  peaks=5
  data <- read.FCS(FACSC20141216_files[[1]], transformation=FALSE,alter.names = T)
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])
  plot(log(to_plot$x),log(to_plot$y))
  kern<-kde2d(log(to_plot$x),log(to_plot$y),n = resol)
  mask<-kern$z*as.numeric(kern$z>quantile(kern$z,tolerance))
  to_plot<-data.frame(x=rep(kern$x,resol),y=rep(kern$y,each=resol),z=melt(kern$z)[,3])
  
#this doesn't make sense to me, it should be colMeans...but is not ... i need some sleep#
  local_maximas_x<-kern$x[which(diff(sign(diff(rowMeans(kern$z*mask))))==-2)+1][1:peaks]
   
  ggplot(to_plot,aes(x,y))+
    geom_tile(aes(fill=z))+
    coord_equal()+
    labs(x='log(PI-H)',y='log(DAPI-H)')+
    geom_vline(xintercept=local_maximas_x,linetype=2,aes(colour=local_maximas_x))
  
  
