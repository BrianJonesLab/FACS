require(flowCore)

arf_files<-dir('C:/Users/lab/Documents/JUAN PABLO/FACS/FACSVerse 16 Dec 2014/External JuanM20141216/',pattern = 'arf',full.names = T)
Col0_files<-dir('C:/Users/lab/Documents/JUAN PABLO/FACS/FACSVerse 16 Dec 2014/External JuanM20141216/',pattern = 'Col0',full.names = T)

arf<- lapply(arf_files,function(x) read.FCS(x, transformation=FALSE))
Col0<- lapply(arf_files,function(x) read.FCS(x, transformation=FALSE))


#check some properties of each object#
summary(arf[[1]])
arf[[1]]@exprs[1:10,1:10]

#check if this is higher than 0 #
keyword(object = arf[[1]],paste0('$P1E'))


#alter names for mor "R friendly names"#
#note ...the function alter the names first and after that looks for a column pattern#
arf<- lapply(arf_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))
arf[[1]]@exprs[1:5,1:5]





up<-seq(0.6,1.05,0.01)
down<-seq(0.5,1,0.01)
combinations<-expand.grid(up,down)

combinations<-combinations[combinations[,1]>combinations[,2],]

require(plyr)

results<-adply(combinations,1,function(x){
  selection<-data[data$PI.H/data$DAPI.H>=as.numeric(x[2]) & round(data$PI.H/data$DAPI.H,3)<=as.numeric(x[1]),]
  localm <- localMax(density(selection$PI.H,kernel='gaussian')$y)
  local_maximas<-density(selection$PI.H)$x[localm][1:3]
})

final<-sapply(3:7,function(x) table(round(results[,x],0))[which.max(table(round(results[,x],0)))])
