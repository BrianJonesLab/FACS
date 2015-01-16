####Install Bioconductor libraries and required packages####
source("http://bioconductor.org/biocLite.R")
biocLite()

biocLite(c("GenomicFeatures", "AnnotationDbi","flowCore",'flowViz'))


#####Read some flow data####

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

#Visualization#

require('flowViz')
#DON'T DO THIS...TO MANY GRAPHICS#
#plot(arf[[1]])

par()
plot(arf[[1]],c("PI.H", "DAPI.H"))

require(ggplot2)

#example of plot#
ggplot(as.data.frame(arf[[1]]@exprs),aes(PI.H,DAPI.H))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

#example of histogram#
ggplot(as.data.frame(arf[[1]]@exprs),aes(PI.H))+
  geom_density()+
  scale_x_log10()
data <-as.data.frame(arf[[1]]@exprs[,c('PI.H','DAPI.H')])

#did not work :(#
##
require(manipulate)
manipulate(
   ggplot(data[data$PI.H/data$DAPI.H>=down & round(data$PI.H/data$DAPI.H,3)<=up,],aes(PI.H))+
    geom_density(fill='red',alpha=.4)+
    scale_x_log10(limits=c(100, 80000)),  
  down=slider(min=.7,max=0.89,step=.01),
  up=slider(min = 0.9,max = 1.5,step = .01))

#end...by now#

manipulate(
  ggplot(data[data$PI.H/data$DAPI.H>=down & round(data$PI.H/data$DAPI.H,3)<=up,],aes(PI.H,DAPI.H))+
    geom_point(fill='red',alpha=.4)+
    scale_y_log10(limits=c(100, 80000))+
    scale_x_log10(limits=c(100, 80000)),  
  down=slider(min=.7,max=0.89,step=.01),
  up=slider(min = 0.9,max = 1.5,step = .01))




test_threshold <- function(up,down){
threshold <- c(up,down)

selection<-data[data$PI.H/data$DAPI.H>=threshold[1] & round(data$PI.H/data$DAPI.H,3)<=threshold[2],]

ggplot(selection,aes(PI.H,DAPI.H))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

peaks

# histograma <-ggplot(selection,aes(PI.H))+
#   geom_density(fill='red',alpha=.5)+
#   scale_x_log10()+
#   labs(title='Example of peaks')
}

test_threshold(up = .95,down = 1)


a<-density(selection$PI.H)$y

#Visualize local maximas#
require(zoo)

plot_function <- function(up, down){
  
#   localMax<- function(x) {
#     # Use -Inf instead if x is numeric (non-integer)
#     y <- diff(c(-.Machine$integer.max, x)) > 0L
#     rle(y)$lengths
#     y <- cumsum(rle(y)$lengths)
#     y <- y[seq.int(1L, length(y), 2L)]
#     if (x[[1]] == x[[2]]) {
#       y <- y[-1]
#     }
#     y
#   }
  
  localMax<-function(tt) which(diff(sign(diff(tt)))==-2)+1
  
  selection<-data[data$PI.H/data$DAPI.H>=down & round(data$PI.H/data$DAPI.H,3)<=up,]
  localm <- localMax(density(selection$PI.H,kernel='gaussian')$y)
  local_maximas<<-density(selection$PI.H)$x[localm][1:5]
  kernel_estimator<-data.frame(x=density(selection$PI.H,kernel='gaussian')$x,y=density(selection$PI.H,kernel='gaussian')$y)
  
  ggplot(selection,aes(PI.H))+
    geom_line(data = kernel_estimator,aes(x=x,y=y))+
    scale_x_log10(limits=c(10, 80000))+
    geom_vline(xintercept=local_maximas)+
    labs(title=paste0(1:3,paste0(' peak:',local_maximas),collapse = ','))
}
manipulate(plot_function(up,down),  
  down=slider(min =.7,max = 0.89,step = .01),up=slider(min = 0.9, max = 1.5,step = .01)
)





#####Optimize Local maximas based in abundance of peaks for different combinations####
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
