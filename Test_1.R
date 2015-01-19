####Install Bioconductor libraries and required packages####
source("http://bioconductor.org/biocLite.R")
## Use biocLite()

biocLite(c("GenomicFeatures", "AnnotationDbi","flowCore",'flowViz'))

#####Read FACS data####

require(flowCore)

# in some cases we use pattern = 'arf', to select some of the files
# Col0_files<-dir('C:/Users/lab/Documents/JUAN PABLO/FACS_Data/External JuanM20141216/',pattern = 'Col0',full.names = T)

FACSC20141216_files <- dir('C:/Users/lab/Documents/JUAN PABLO/FACS_Data/External JuanM20141216/',full.names = T)



#read the files
FACSC20141216 <- lapply(FACSC20141216_files,function(x) read.FCS(x, transformation=FALSE))


#check properties of each object#
summary(FACSC20141216[[1]])
FACSC20141216[[1]]@exprs[1:10,1:10]

# in the library of FACS says that we need to check if there is a parameter > 0, but we didn't find any
#check if this is higher than 0 #
#this is a test of one of them
keyword(object = FACSC20141216[[1]],paste0('$P1E'))


#alter names for mor "R friendly names"#
#note ...the function alter the names first and after that looks for a column pattern#
FACSC20141216<- lapply(FACSC20141216_files,function(x) read.FCS(x, transformation=FALSE,alter.names = T))
FACSC20141216[[1]]@exprs[1:5,1:5]

#Visualization#

require('flowViz')
#DON'T DO THIS...TO MANY GRAPHICS#
#plot(FACSC20141216[[1]])

FACSC20141216[[1]]

par()
plot(FACSC20141216[[1]],c("PI.H", "DAPI.H"))
#Would be better to use log scale in both axes


require(ggplot2)

#example of plot# the points in the middle are the nuclei
ggplot(as.data.frame(FACSC20141216[[1]]@exprs),aes(PI.H,DAPI.H))+
  geom_point()+
  scale_x_log10()+
  scale_y_log10()

#example of histogram# but if we plot all, we got the peaks and noise
ggplot(as.data.frame(FACSC20141216[[1]]@exprs),aes(PI.H))+
  geom_density()+
  scale_x_log10()
data <-as.data.frame(FACSC20141216[[1]]@exprs[,c('PI.H','DAPI.H')])

# Using manipulate we treshole the data
require(manipulate)

manipulate(
  ggplot(data[data$PI.H/data$DAPI.H>=down & round(data$PI.H/data$DAPI.H,3)<=up,],aes(PI.H,DAPI.H))+
    geom_point(fill='red',alpha=.4)+
    scale_y_log10(limits=c(100, 80000))+
    scale_x_log10(limits=c(100, 80000)),  
  down=slider(min=.7,max=0.89,step=.01),
  up=slider(min = 0.9,max = 1.5,step = .01))

# And then the same with the histogram, the peaks are much more clear

manipulate(
   ggplot(data[data$PI.H/data$DAPI.H>=down & round(data$PI.H/data$DAPI.H,3)<=up,],aes(PI.H))+
    geom_density(fill='red',alpha=.4)+
    scale_x_log10(limits=c(100, 80000)),  
  down=slider(min=.7,max=0.89,step=.01),
  up=slider(min = 0.9,max = 1.5,step = .01))

# Then the next step is to find a way to optimize this threshold checking all the possible combination
# and get the local max, but using the 2D histogram information.

#Next file is Pea_analysis