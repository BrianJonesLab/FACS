library(raster)

FACSC20141216_files <- dir('.','.fcs',full.names = T)
FACSC20141216 <- lapply(FACSC20141216_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))

findPeaks <- function(data) {
  to_plot <- data.frame(x=data@exprs[,'PI.H'],y=data@exprs[,'DAPI.H'])
  to_plot <- data_ <- log(to_plot)
  
  coordinates(data_) <- ~ x + y
  
  template <- raster(res=0.25,ext=extent(data_))
  test <- rasterize(data_,template,data_@coords[,1],fun='count')
  
  test_ <- na.exclude(as.data.frame(cbind(xyFromCell(test,cell = 1:ncell(test)),test[])))
  lines <- test_[test_[,3] > quantile(test_[,3],0.99,na.rm = T),1:2]
  plot(test,xlab='log(PI-H)',ylab='log(DAPI-H)')
  sapply(lines[,'x'],function(x) abline(v=x))
  sapply(lines[,'y'],function(x) abline(h=x))
  #ggplot(test_,aes(x,y,fill=V3)) + geom_raster()  + xlab('log(PI-H)') + ylab('log(DAPI-H)') + geom_vline(aes(xintercept=x),lines,linetype='longdash',colour='red') + geom_hline(aes(yintercept=y),lines,linetype='longdash',colour='red')
  cat(c("Peaks' coordinates:",apply(as.matrix(lines),1,function(x) paste0(round(x,2),collapse = ', '))),sep = '\n')
}

findPeaks(FACSC20141216[[2]])
