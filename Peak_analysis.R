require(flowCore)

arf_files<-dir('C:/Users/lab/Documents/JUAN PABLO/FACS/FACSVerse 16 Dec 2014/External JuanM20141216/',pattern = 'arf',full.names = T)
Col0_files<-dir('C:/Users/lab/Documents/JUAN PABLO/FACS/FACSVerse 16 Dec 2014/External JuanM20141216/',pattern = 'Col0',full.names = T)

arf<- lapply(arf_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))
Col0<- lapply(arf_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))

localMax<-function(tt) which(diff(sign(diff(tt)))==-2)+1
require(plyr)

results_arf<-lapply(1:length(arf_files),function(j){
  
data <-as.data.frame(arf[[j]]@exprs[,c('PI.H','DAPI.H')])

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

final<-sapply(3:7,function(k) table(round(results[,k],0))[which.max(table(round(results[,k],0)))])

})


ARF_RESULTS<-data.frame(file=dir('C:/Users/lab/Documents/JUAN PABLO/FACS/FACSVerse 16 Dec 2014/External JuanM20141216/',pattern = 'arf'),do.call(rbind,lapply(results_arf,function(x) as.numeric(names(x)))))

colnames(ARF_RESULTS) <-c('File',1:5)
ARF_RESULTS
