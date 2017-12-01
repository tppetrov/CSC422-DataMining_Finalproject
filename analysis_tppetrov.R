library(ggplot2)
num.id <- 6

# change this to your dir
determinedPath = 'C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/'

dt.id <-paste("DT", as.character(num.id), sep="")
# file shen gave us, this file contained summary information
filename <- paste(determinedPath,dt.id, "/same_small_feature/", dt.id,"_whole_measure.csv",sep="")
data <- read.csv(filename)

# remove the unstart ppl
data <- data[-which(is.na(data$L1.BadSteps)),]
###################
## create for loop to load all the data#####
###################

# files to be read
# I actually have the different cluster types separated so thats why the .linkType
#specify the file you are looking for for statistical analysis
filenames.list <- c('gaussian.cluster4.hclust.average',
                    'gaussian.cluster4.hclust.median', 'gaussian.cluster4.hclust.ward.D2')
filename <- paste(determinedPath,"/DT6/same_small_feature/kernel/cluster/", 'gaussian.cluster4.hclust.complete', ".csv",sep="")
cluster.result <- read.csv(filename)

cluster.result <- cluster.result[-1]
r <- cluster.result

##########################
## put all different links into one variable/matrix
## r
#########################
for (file in filenames.list){
  filename <- paste( determinedPath,"/DT6/same_small_feature/kernel/cluster/", file, ".csv",sep="")
  #read
  cluster.result <- read.csv(filename)
  #remove unecessary data
  cluster.result <- cluster.result[-(1:2)]
  ##########################
  ## put all the data into one variable/matrix
  r <- cbind(r, cluster.result)
}
colns = colnames(r)
colnames = colns[-1]
data1 <- merge(r, data, by=c("student"))

# displays and saves to a .txt file the:
# aov() - statistical analysis
#       Pr(>F) is the p-value
#
# table() - how many students per cluster
#

displayData <- function(){
  fileconn <- file('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/DT6/same_small_feature/kernel/cluster/aovOutput/gaussiansas10_cluster4.txt')
  for(myCol in colnames){
    sink(fileconn)
    cat("Kernel: ", myCol, "\n")
    results <- aov(data1$postTest ~ data1[,myCol])
    mytable <- table(data1[myCol])
    #view avona
    print(summary(results))
    # check clusters size
    #df <- data.frame(group = c(rep("1"), rep("2"), rep("3"), rep("4")))
    print(mytable)
  }
  sink(file = NULL)
}
