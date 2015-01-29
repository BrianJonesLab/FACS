require(flowCore)
require(MASS)
require(ggplot2)
require(reshape2)
FACSC20141216_files <- dir('.','.fcs',full.names = T)

# findPeaks <- function(file,
#                       peaks=6,
#                       resol=500) {

  resol=500
  data <- read.FCS(FACSC20141216_files[[1]], transformation=FALSE,alter.names = T)
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])
  plot(log(to_plot$x),log(to_plot$y))
  kern<-kde2d(log(to_plot$x),log(to_plot$y),n = resol)
  tol_x <- log(abs((exp(kern$x[1])-exp(kern$x[2])))*c(resol/4))
  tol_y <- log(abs((exp(kern$y[1])-exp(kern$y[2])))*c(resol/4))
  
  to_plot<-data.frame(x=rep(kern$x,resol),y=rep(kern$y,each=resol),z=melt(kern$z)[,3])
  first<-to_plot[which.max(to_plot$z),]
  
  lower<-c(log(exp(first$x)-exp(tol_x)))
  upper<-c(log(exp(first$x)+exp(tol_x)))
  
  second <- to_plot[!to_plot$x>lower&to_plot$x<upper,][which.max(to_plot[!to_plot$x>lower&to_plot$x<upper,'z']),]
  
  a<-ggplot(to_plot,aes(x,y))+
    geom_tile(aes(fill=z))+
    coord_equal()+
    labs(x='log(PI-H)',y='log(DAPI-H)')+
    geom_vline(xintercept=first$x,linetype=2)+
    geom_hline(yintercept=first$y,linetype=2)
  a
  
  a+geom_vline(xintercept=second$x,linetype=2,colour='red')+
    geom_hline(yintercept=second$y,linetype=2,colour='red')
   
    
  

