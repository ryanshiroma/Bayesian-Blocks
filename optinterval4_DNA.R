
# R version of Dr. Jackson's MATLAB code from week 1 lecture notes



# NEW OBJECTIVE FUNCTION
# A function that computes the objective function value of a block 
# N is the length of the block and M is the number of data points in the block 

newobjective = function(N,M,c){
  value = N*(log(N/M)) - c
  value[is.nan(value)]<- -Inf
  return(value)
}

newobjective = function(N,M,c){
  value = N*log(N/M) + (M-N)*log((M-N)/M)-c
  value[M==N]=-Inf
  return(value)
}

newobjective = function(N,M,c){
  value = StirlingAprox(M,N)+N*log((N)/M)+(M-N)*log(1-(N)/M)-c 
  value[M==N]=-Inf
  #value[is.na(value)]<- -Inf
  return(value)
}

StirlingAprox <- function(n,k){
  if (n==k){
    return(1)
  }
  return (n*log(n)-(n-k)*log(n-k)-k*log(k)+0.5*log(n/(2*pi*(n-k)*k)))
  
}

#S3 generic plot function 

"plot.BB" <- function(x,points=TRUE,blocks=TRUE,hist=FALSE,bins=100) {
  data       <- x$data
  n          <- length(data)
  intensity  <- round(x$N/x$A,4)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2  
  edges      <- cumsum(c(lowerbound,x$A))
  
  
  #initialize empty plot
  if (hist==TRUE && points==FALSE){
    plot(NA, ylim = c(0, max(x$intensity)*1.4), xlim = c(lowerbound,upperbound),
         ylab = "CG Ratio",xlab="Position",main="Bayesian Blocks")   
    legend(x=x$data[1],y=max(x$intensity)*1.4,legend= c("Binned Data Points","Blocks"), lty=c(1,1), 
           lwd=c(3,3), col=c("grey80","red")) 
    
  } else{
    plot(NA, ylim = c(0, quantile(intensity,0.90)), xlim = c(lowerbound,upperbound),
         ylab = "CG Ratio",xlab="Position",main="Bayesian Blocks")
    legend(x=x$data[1],y=quantile(intensity,0.90),legend= c("Data Points","Blocks"), lty=c(1,1), 
           lwd=c(3,3), col=c("cornflowerblue","red")) 
    
    
  }
  
  
  #plot histogram if requested
  if (hist==TRUE){
    end    <- x$data[length(x$data)]
    width  <- end/bins
    i      <- 0:bins
    height <- hist(x$data, breaks=c(i,bins+1)*width,plot=FALSE)$counts/width
    rect(xleft=i*width,ybottom=0,xright= (i+1)*width, ytop=height,col="grey80",border="grey40")
    

    
  }
  
  
  #plot individual intensities if requested
  if (points==TRUE){ 
  segments(x0=edges[1:(length(edges)-1)], y0=intensity[1:length(intensity)], 
           x1=edges[2:length(edges)], y1=intensity[1:length(intensity)],
           col = "cornflowerblue", lwd = 3)
  }
  

  #plot constant blocks
  if (blocks==TRUE){
    indices <- which(x$type=="constant")
    segments(x0=edges[x$left[indices]], y0=x$intensity[indices], 
             x1=edges[x$right[indices]+1], y1=x$intensity[indices],
             col = "red", lwd = 3)
    
  }
  

}




#S3 Generic summary function
summary.BB <- function(x) {
  cat(length(x$data), 'data points\n\n',
      length(x$left),'total blocks:\n\t',
      sum(x$type=='constant'),'constant blocks\n\t',
      sum(x$type=='linear'),'linear blocks\n\t',
      sum(x$type=='exponential'),'exponential blocks\n\t',
      sum(x$type=='power'),'power function blocks\n\n',
      floor((length(x$pruned)*100)/length(x$data)),'% of the points were pruned')
}



# OPTINTERVAL FUNCTION 
# An O(N^2) algorithm for finding the optimal partition of N data points on an interval
# INPUT:
# data is the vector ofcells
# N is the number of "observations" per point. This is the same "N" we've been using all along
# c is the block penalty term
# OUTPUT:
# a list of objects pertaining to the "BB" object

optinterval = function(data,N,c){
  n <- length(N)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2
  
  xx <- c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound) # voronoi cell vertices
  A  <- diff(xx) # length of each voronoi cell
  
  opt           <- rep(0,n+1)
  lastchange    <- rep(1,n)
  changeA       <- rep(0,n) 
  changeN       <- rep(0,n)
  optint        <- rep(0,n)
  lasttype      <- rep("None",n)
  lastintensity <- rep(0,n)

  unpruned = NULL
  endobj = rep(0,n)
  optint = rep(0,n)

  for (i in 1:n){
    unpruned = c(unpruned,i)

    changeA[unpruned] = changeA[unpruned]+A[i]
    changeN[unpruned] = changeN[unpruned] + N[i]
    endobj[unpruned] = newobjective(changeN[unpruned],changeA[unpruned],c)
    optint[unpruned] = opt[unpruned] + endobj[unpruned]

    opt[i+1]         <- max(optint[unpruned])
    last = which.max(optint[unpruned])
    lastchange[i]    <- unpruned[last]
    unpruned = unpruned[((optint[unpruned]+c-opt[i+1])>0)];
    
    lastintensity[i] <- sum(N[lastchange[i]:i])/sum(A[lastchange[i]:i])
    lasttype[i]      <- "constant"
    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastintensity"= lastintensity,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  return(getStep(BBdata,n))
}



#this function returns a BB object for the optimal blocks at any given cell location "n" where 1 < n < N
getStep <- function(BBobj,n){
  
  lasttype      <- BBobj$lasttype[1:n]
  lastchange    <- BBobj$lastchange[1:n]
  lastintensity <- BBobj$lastintensity[1:n]
  
  lastchange<-lastchange-1
  
  left          <- vector()
  intensity     <- vector()
  type          <- vector()
  type[1]       <- lasttype[n]
  left[1]       <- lastchange[n]
  intensity[1]  <- lastintensity[n]
  i=1
  
  while (left[i] > 0) {
    left[i+1]      <- lastchange[left[i]]
    type[i+1]     <- lasttype[left[i]]
    intensity[i+1] <- lastintensity[left[i]]
    
    i <- i+1
  }
  left      <- rev(left)+1
  type      <- rev(type)
  intensity <- rev(intensity)
  
  if (length(left) == 1){
    right <- n
  }
  else {
    right <- c(left[2:length(left)]-1,n)
  }
  
  BBdata   <- list("data"         = BBobj$data[1:n],
                   "N"            = BBobj$N[1:n],
                   "A"            = BBobj$A[1:n],
                   "left"         = left,
                   "right"        = right,
                   "type"         = type,
                   "intensity"    = intensity,
                   "opt"          = BBobj$opt[1:n],
                   "lastchange"   = lastchange+1,
                   "lastintensity"= lastintensity,
                   "lasttype"     = lasttype,
                   "pruned"       = BBobj$pruned)
  
  BBobject <- structure(BBdata, class = "BB")
  
  return(BBobject)
}


#Non-homogeneous Poisson Process Generator!!
NHPois <- function(time_length, equation){
  
  expr<-function(x){
    x<-x
    eval(equation)
  }
  current <- 0
  i<-1
  data <- vector()
  maxrate = optim(par=time_length/2,fn=expr,method='L-BFGS-B',lower=0,upper=time_length,control = list(fnscale = -1))$value
  while(current < time_length){
    current<- current+ rexp(n=1,rate=maxrate)
    u <- runif(1)
    
    if (u < expr(current) / ( maxrate)){
      data[i] <- current
      i<-i+1
    }
    
  }
  return(data[c(-length(data))])
}

############################## SCRIPT STARTS HERE ############################################



setwd("~/R")
load("homosapiens.Rdata")
subset=0.002#get 0.2 % of the data from the first half of the chromosome
x=homo.sapiens[seq(1,length(homo.sapiens)/2,length.out=length(homo.sapiens)*subset)]*subset*2
x=x[2:length(x)]
N <- rep(1,length(x))#let all cells be of size 1

#run BB algorithm
start.time <- Sys.time()
model<-optinterval(x,N,3)#00001)
print(Sys.time() - start.time)


#check summary information and plot
summary(model)
plot(model,points=FALSE,hist=TRUE,blocks=TRUE,bins=600)

plot(model$opt)

#plot lengths of every other block
l=model$right-model$left
l2=l[seq(1,length(l),2)]
plot(l2[10:length(l2)])
