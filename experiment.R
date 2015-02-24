require(flowCore)
require(MASS)
require(ggplot2)
require(reshape2)
require(Rwave)
require(EBImage)

FACSC20141216_files <- dir('FACS_02FEB2015/.','.fcs',full.names = T)

findPeaks1 <- function(file,tolerance=.995,peaks=5, resol=250,plot=F) {

  data <- read.FCS(file, transformation=FALSE,alter.names = T)
  mask<-data@exprs[,'PI.H']/data@exprs[,'DAPI.H']>=0.8&data@exprs[,'PI.H']/data@exprs[,'DAPI.H']<=1.2
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])[mask,]
#   plot(log(to_plot$x),log(to_plot$y))
  kern<-kde2d(log(to_plot$x),log(to_plot$y),n = resol)

  mask<-kern$z*as.numeric(kern$z>quantile(kern$z,tolerance))
  
  to_plot_kern<-data.frame(x=rep(kern$x,resol),y=rep(kern$y,each=resol),z=melt(kern$z)[,3])
  
#this doesn't make sense to me, it should be colMeans...but is not ... i need some sleep#
#REMINDER: always read the help file!... 
#kern2d creates a matrix were x are rows....now peace reigns again#
  loc_max<-function(DATA) which(diff(sign(diff(DATA)))==-2)+1
  loc_max_posit_x<-loc_max(rowSums(kern$z*mask))
  local_maximas_x<-kern$x[loc_max_posit_x]
  if(length(na.exclude(local_maximas_x))<peaks) {
    warning(immediate. = T,call.=FALSE,paste('only',length(na.exclude(local_maximas_x)),' peaks were higher that the tolerance value'))
    local_maximas_x<-local_maximas_x[complete.cases(local_maximas_x)]
    peaks<-length(local_maximas_x)
  }
  loc_max_posit_y<-sapply(loc_max_posit_x[order(local_maximas_x)[1:peaks]],function(x) tmp<-which.max(kern$z[x,]))
  
  local_maximas_y<-sapply(loc_max_posit_x[order(local_maximas_x)[1:peaks]],function(x) {
    tmp<-which.max(kern$z[x,])
    kern$y[tmp]
  })
  local_maximas_x<-kern$x[loc_max_posit_x][1:peaks]

# ggplot(to_plot,aes(log(x),log(y)))+
#   geom_point()+
#   coord_equal()+
#   labs(x='log(PI-H)',y='log(DAPI-H)')
if(plot){
    print(ggplot(to_plot_kern,aes(x,y))+
    geom_tile(aes(fill=z))+
    coord_equal()+
    labs(x='log(PI-H)',y='log(DAPI-H)')+
    geom_vline(xintercept=local_maximas_x,linetype=2,colour=local_maximas_x)+
    geom_hline(yintercept=local_maximas_y,linetype=2,colour=local_maximas_x))
}
# pix_size_x <- abs(kern$x[1]- kern$x[2])
# pix_size_y <- abs(kern$y[1]- kern$y[2])
# 
# 
# buffer <- mapply(FUN = function(x,y) kern$z[x-x*(buffer_size/pix_size_x):x-buffer_size/pix_size_x,y],loc_max_posit_x,loc_max_posit_y)
results <- round(data.frame(x_value=local_maximas_x,y_value=local_maximas_y),3)
colnames(results) <- c('log(PI-H)','log(DAPI-H)')

cat(paste0('Peaks for sample: \n##', basename(file), '##'),'\n')
print(results)
cat('#################################')
}  

findPeaks1(file = FACSC20141216_files[[4]],tolerance = .98,resol = 100,peaks = 5,plot = T)

