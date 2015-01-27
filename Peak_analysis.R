require(flowCore)

FACSC20141216_files <- dir(pattern = 'fcs',full.names = T)

FACSC20141216 <- lapply(FACSC20141216_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))

localMax<-function(tt) which(diff(sign(diff(tt)))==-2)+1
require(plyr)

results_arf<-lapply(1:length(FACSC20141216_files),function(j){
  
data <-as.data.frame(FACSC20141216[[j]]@exprs[,c('PI.H','DAPI.H')])

up<-seq(0.6,1.05,0.01)
down<-seq(0.5,1,0.01)
combinations<-expand.grid(up,down)
combinations<-combinations[combinations[,1]>combinations[,2],]
combinations<-combinations[combinations[,1]-combinations[,2]>=0.1,]


results<-adply(combinations,1,function(i){
  selection<-data[data$PI.H/data$DAPI.H>=as.numeric(i[2]) & round(data$PI.H/data$DAPI.H,3)<=as.numeric(i[1]),]
  localm <- localMax(density(selection$PI.H,kernel='gaussian')$y)
  local_maximas<-density(selection$PI.H)$x[localm][1:5]
})

#So here we test all the combination of the threshold and then we use the peaks that repeat more,
# but if you have other idea, use it.

final<-sapply(3:7,function(k) table(round(results[,k],0))[which.max(table(round(results[,k],0)))])

})

#MArio was tire and wrote this long chorizo

FACSC20141216_RESULTS<-data.frame(file=dir(pattern = 'fcs'),do.call(rbind,lapply(results_arf,function(x) as.numeric(names(x)))))

colnames(FACSC20141216_RESULTS) <-c('File',1:5)
FACSC20141216_RESULTS

#This is a good approach, but we force to find 5 peaks and probably will be better to find the real ono.
#Also I need the number of how many points we have on each peak, so I can know the % of the total.
