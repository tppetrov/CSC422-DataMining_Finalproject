#Import data
load.distanceK <- function(granularity, select.prefix){
  prefix <- 'C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/'
  filename <- paste(prefix, dt.id, '/same_small_feature/',select.prefix,'/', granularity, '.distanceG10.csv',sep="") #indicate which file is distance matrix here
  dis.matrix <- read.csv(filename)
  row.names(dis.matrix) <- names(dis.matrix)
  return(as.dist(dis.matrix))
}

#load a distance matrix from desired location
load.distance <- function(granularity, select.prefix){
  prefix <- 'C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/'
  filename <- paste(prefix, dt.id, '/same_small_feature/',select.prefix,'/', granularity, '.distance.csv',sep="")
  dis.matrix <- read.csv(filename)
  row.names(dis.matrix) <- names(dis.matrix)
  return(as.dist(dis.matrix))
}

#clusters distance matrix into specified number of clusters based off given method
clustering <- function(d, k, method, link="ward.D2"){
  # hierarchical clustering
  if (method =='hclust'){
    fit <- hclust(d, method=link) 
    plot(fit)
    #         rect.hclust(fit, k = 3, border = "red")
    #         source("http://addictedtor.free.fr/packages/A2R/lastVersion/R/code.R")
    #         op = par(bg = "#EFEFEF")
    #         A2Rplot(fit, k = 3, boxes = FALSE, col.up = "gray50", col.down = c("#FF6B6B", "#4ECDC4", "#556270","#337DFF"))
    groups <- cutree(fit, k)
    cluster.result <- data.frame('student'= names(groups), 'ward'=as.numeric(groups))
  }
  else if(method=='k.mediods'){
    library('cluster')
    # k-mediods clustering
    #         asw <- numeric(5)
    #         for (k in 2:5){
    #           asw[k] <- pam(d, k) $ silinfo $ avg.width
    #         }
    #         k.best <- which.max(asw)
    #         plot(1:5, asw, type= "h", main = "pam() clustering assessment",
    #              xlab= "k (# clusters)", ylab = "average silhouette width")
    #         axis(1, k.best, paste("best",k.best,sep="\n"), col = "red", col.axis = "red")
    #     groups <- pam(d, k.best)
    groups <- pam(d, k)
    cluster.result <- data.frame('student'= names(groups$clustering), 'k.mediods'= as.numeric(groups$clustering))
  }
  else if(method=='spectral'){
    library(kernlab)
    library(stats)
    
    rbf <- rbfdot(sigma = 0.05)
    S <-kernelMatrix(rbf, as.matrix(d))
    A <- make.affinity(S, 3)
    D <- diag(apply(A, 1, sum))
    
    # unormalized 
    U <- D - A
    
    # simple Laplacian
    L <- diag(nrow(as.matrix(d))) - solve(D) %*% A  
    
    # normalized Laplacian
    L <- (D %^% (-1/2)) %*% A %*% (D %^% (-1/2))  
    
    evL <- eigen(U, symmetric=TRUE)
    evL <- eigen(L, symmetric=TRUE)
    
    
    signif(evL$values,2)
    plot(1:10, rev(evL$values)[1:10], log="y")
    abline(v=5.5, col="red", lty=2) # there are 5 clusters as expected
    
    k   <- 5
    Z   <- evL$vectors[,(ncol(evL$vectors)-k+1):ncol(evL$vectors)]
    
    km <- kmeans(Z, centers=k, nstart=5)
    
    cluster.result <- data.frame('student' = colnames(as.matrix(d)), 'spectral'=as.numeric(km$cluster))
  }
  else if (methpd =='dbscan'){
    obj <- dbscan(d, eps=0.5)
  }
  return(cluster.result)
}

#I only worked at the problem level, this should not impact you. This was how I loaded my data and specified which file I needed.
load.student.data <- function(dt.id, granularity){
  prefix <<- 'C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/'
  
  # load data set
  if (granularity=='last.level'){
    granularity <- unlist(strsplit(granularity, split='.', fixed=TRUE))[2]
  }
  if (granularity=='problem.level'){
    granularity <- unlist(strsplit(granularity, split='.', fixed=TRUE))[1]
  }
  
  filename <- paste(prefix, dt.id,'/same_small_feature/', dt.id, '_continuous_', granularity,'.csv' , sep="")
  data <- read.csv(filename)
  
  return(data)
}

#This is where I normalized my data, this was preprocessing
standardization <- function(data, feature.list){
  
  mu.features <- rep(0, length(feature.list))
  sigma.features <- rep(0, length(feature.list))
  
  for (fID in 1:length(feature.list)){
    mu.features[fID] <- mean(data[,feature.list[fID]])
    sigma.features[fID] <- sd(data[,feature.list[fID]])
    data[,feature.list[fID]] <- (data[,feature.list[fID]] - mu.features[fID]) / sigma.features[fID]
  }
  return(data)
  
}

#I defined multiple Kernel Functions and the eucledian distance.  These Equations were used to weight the DTW algorithm.
euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))

linearl.kernel.dist <- function(x1, x2) abs(t(x1)%*%x2 + 100) #make abs if this doesn't work
#had a large c bias of 1000 and got all jank

polynomial.kernel.dist <- function(x1, x2) (t(x1)%*%x2 + 20)^2 #gabe had one that worked with just +100 using inner product, will repeat using linear

gaussian.kernel.dist <- function(x1, x2) exp( -(sqrt(sum((x1 - x2) ^ 2) ^ 2)) / 10 ) #there is a sigma^2 on the bottom that is set to 1 therefor we were dividing by 2

laplace.kernel.dist <- function(x1, x2) exp( -(sqrt(sum((x1 - x2) ^ 2))) / 100 )

sigmoid.kernel.dist <- function(x1, x2) tanh(t(x1)*x2)

#This is the actual DTW Algorithm.  Areas marked with modify were points where I needed to change the weighted function that modifies DTW and helps pick values.
dtw.pair.dist <- function(obj1, obj2){
  w.h <- 1      # 1   Horizontal
  w.v <- 1      # 2   Vertival   
  w.d <- 2      # 3   Diagonalc
  
  dmatrix <- matrix(Inf, nrow = nrow(obj1), ncol = nrow(obj2))
  
  # first step
  #******* modify*******
  #dmatrix[1,1] <- euc.dist(as.numeric(obj1[1,]), as.numeric(obj2[1,]))
  dmatrix[1,1] <- laplace.kernel.dist(as.numeric(obj1[1,]), as.numeric(obj2[1,]))
  
  # second step
  for( i in 2:nrow(obj1)){
    #******* modify*******
    dmatrix[i,1] <- dmatrix[i-1,1] + w.v * laplace.kernel.dist(as.numeric(obj1[i,]), as.numeric(obj2[1,]))
  }
  
  for(j in 2:nrow(obj2)){
    #******* modify*******
    dmatrix[1, j] <- dmatrix[1, j-1] +  w.h * laplace.kernel.dist(as.numeric(obj1[1,]), as.numeric(obj2[j,]))
  }
  # synamic programming
  for( i in 2:nrow(obj1)){
    for(j in 2:nrow(obj2)){
      
      #******* modify*******
      dvalue <- laplace.kernel.dist(as.numeric(obj1[i,]), as.numeric(obj2[j,]))
      
      choice1 <- dmatrix[i, j-1] + w.v*dvalue
      choice2 <- dmatrix[i-1, j] + w.h*dvalue
      choice3 <- dmatrix[i-1, j-1] + w.d*dvalue
      dmatrix[i,j] = min(choice1, choice2, choice3)
    }
  }
  return(dmatrix[nrow(obj1), nrow(obj2)]/(nrow(obj1)+nrow(obj2))) 
}

#This calculates the similarity matrix
dtw.dist <- function(data, feature.list){
  student.list <- as.character(unique(data$student))
  nstu<- length(student.list)
  disc.matrix <- matrix(0, nrow = nstu, ncol = nstu)
  
  rownames(disc.matrix) <- student.list
  colnames(disc.matrix) <- student.list
  
  # calculate the similarity matrix
  for(i in 1:(nstu-1)){
    stu1 <- data[which(data$student==student.list[i]),]
    stu1 <- stu1[,feature.list]
    
    for(j in (i+1):nstu){
      stu2 <- data[which(data$student==student.list[j]),]
      stu2 <- stu2[,feature.list]
      disc.matrix[i,j] <- dtw.pair.dist(stu1, stu2)
      disc.matrix[j,i] <- disc.matrix[i,j]
    }
  }
  
  return(as.dist(disc.matrix))
}


#Essentially my main function.
kernel.clustering.process <-function(){
  
  dt.id <- 'DT6'
  granularity <- 'problem'
  
  data <- load.student.data(dt.id, granularity)
  
  # generate the feature list used for calcuating distance
  if(granularity=='problem'){
    start.feature <- 4
  }else{
    start.feature <- 3}
  
  
  feature.list <- as.character(names(data)[start.feature : ncol(data)])
  
  #Step 1, Standardize
  
  data <- standardization(data, feature.list)
  
  #Step 2, DTW distance
  
  d <- dtw.dist(data, feature.list)
  filename <- paste('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/', dt.id, '/same_small_feature/kernel/', granularity,'.distanceL100.csv',sep="")
  write.csv(as.matrix(d), file = filename, row.names = FALSE)
  
  
  #Step 3, Clustering
  #loading of data
  num.cluster <- 4
  select.prefix <- 'pca_data2'
  kernel.type <- 'original'
  d <- load.distanceK('problem', 'kernel')
  #The difference between these two was whether I was using a Kernel or Eucledian distance.  Overall it was me being lazy
  d <- load.distance('problem', select.prefix)
  
  #clustering based off different statistical key points.
  cluster.result1 <- clustering(d, num.cluster, 'hclust', link = 'ward.D2')
  names(cluster.result1) <- c('student', 'problem.ward')
  #Saves the plot
  filename <- paste('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/DT6/same_small_feature/kernel/cluster/', kernel.type,'.cluster',num.cluster,'.hclust.ward.D2.csv',sep="")
  write.csv(as.matrix(cluster.result1), file = filename, row.names = TRUE)
  
  cluster.result2 <- clustering(d, num.cluster, 'hclust', link = 'average')
  names(cluster.result2) <- c('student', 'problem.average')
  filename <- paste('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/DT6/same_small_feature/kernel/cluster/', kernel.type,'.cluster',num.cluster,'.hclust.average.csv',sep="")
  write.csv(as.matrix(cluster.result2), file = filename, row.names = TRUE)
  
  cluster.result3 <- clustering(d, num.cluster, 'hclust', link = 'median')
  names(cluster.result3) <- c('student', 'problem.median')
  filename <- paste('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/DT6/same_small_feature/kernel/cluster/', kernel.type,'.cluster',num.cluster,'.hclust.median.csv',sep="")
  write.csv(as.matrix(cluster.result3), file = filename, row.names = TRUE)
  
  cluster.result4 <- clustering(d, num.cluster, 'hclust', link = 'complete')
  names(cluster.result4) <- c('student', 'problem.complete')
  filename <- paste('C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/DT6/same_small_feature/kernel/cluster/', kernel.type,'.cluster',num.cluster,'.hclust.complete.csv',sep="")
  write.csv(as.matrix(cluster.result4), file = filename, row.names = TRUE)
  
  
  
  # store result
  
  
  #Step 4, Evaluation
  
  # load evaluation data set
  dt.id <- 'DT6'
  filename <- paste("C:/Users/Thomas/Documents/NCSU_Classes/Research/Chi/clustering/clustering/",dt.id, "/same_small_feature/", dt.id,"_whole_measure.csv",sep="")
  data <- read.csv(filename)
  # remove the unstart ppl
  data <- data[-which(is.na(data$L1.BadSteps)),]
  
  
  # load clustering result
  filename <- ''
  cluster.result <- read.csv(filename)
  data1 <- merge(data, cluster.result1, by=c("student"))
  
  
  # statistical analysis
  result <- aov(data1$postTest ~ ward,data1)
  print(summary(result))
  #Made a second program for analysis
  
}

