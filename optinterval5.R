

constantobjective = function(N,A,c){
  value = N*(log(N/A)) - c
  value[is.nan(value)]<-0
  return(value)
}



#######linear intensity fitter

linfit <-function(t,c){
  
  if (length(t)==0){ # if no data is recieved, set cost to 0
    print("nothin")
    return(list("a"=0,"b"=0,"cost"=0))
  }  
  
  if (length(t)==1){ #if only one point is supplied, treat as constant intensity
    #print("just one")
    return(list("a"=0,"b"=t[1],"cost"=-c))
  }
  
  
  epsilon  <- 1
  i        <- 1
  M        <- t[length(t)]-t[1]
  S        <- (t[length(t)]^2-t[1]^2)/2
  
  
  
  if (length(t)==2){ #if two points are provided, create a perfect fit
    cost=sum(log(diff(1/diff(c(0,t)))*t+1/t[1]))-diff(1/diff(c(0,t)))*S-1/t[1]*M-c
    return(list("a"=diff(1/diff(c(0,t))),"b"=1/t[1],"cost"=cost))
  }
  
  #start with some initial values for a and b
  coef     <- matrix(c(0,M/length(t)),nrow=2)
  new.cost <- -Inf
  
  while (abs(epsilon) >10e-4 && i<100){
    old.cost <- new.cost
    
    #first derivatives
    f=matrix(c(sum(t/(coef[1]*t+coef[2]))-S,sum(1/(coef[1]*t+coef[2]))-M),nrow=2)
    
    #hessian values
    fa <- -sum((t/(coef[1]*t+coef[2]))^2)
    fb <- -sum(t/(coef[1]*t+coef[2])^2)
    ga <- fb
    gb <- -sum(1/(coef[1]*t+coef[2])^2)
    
    #create the hessian
    hess <- matrix(c(fa,fb,ga,gb),nrow=2,ncol=2,byrow=TRUE) 
    
    #run newtons method
    coef=coef-solve(hess)%*%f 
    
    #calculate the new cost
    new.cost <- sum(log(coef[1]*t+coef[2]))-coef[1]*S-coef[2]*M
    epsilon  <- new.cost-old.cost
    
    i <- i+1
    
  }
  
  list("a"=coef[1],"b"=coef[2],"cost"=(new.cost-c))
  
}

#S3 generic plot function 

"plot.BB" <- function(x,plot,bins=100,ylim) {
  data       <- x$data
  n          <- length(data)
  intensity  <- round(x$N/x$A,4)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2  
  edges      <- cumsum(c(lowerbound,x$A))
  
  
  #initialize empty plot
  if ("hist" %in% plot && !("points" %in% plot)){
    plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),
         ylab = "Intensity",xlab="Time",main="Bayesian Blocks")   
    grid(NA,NULL)#add grid for y-axis
    legend(x=x$data[1],y=ylim[2]*0.9,legend= c("Binned Data Points","Blocks"), lty=c(1,1), 
           lwd=c(3,3), col=c("grey70","red")) 
    
  } else{
    plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),
         ylab = "Intensity",xlab="Time",main="Bayesian Blocks")
    grid(NA,NULL)#add grid for y-axis
    legend(x=x$data[1],y=ylim[2]*0.9,legend= c("Data Points","Blocks"), lty=c(1,1), 
           lwd=c(3,3), col=c("cornflowerblue","red")) 
  }
  
  
  #plot histogram if requested
  if ("hist" %in% plot){
    end    <- x$data[length(x$data)]
    width  <- end/bins
    i      <- 0:bins
    height <- hist(x$data, breaks=c(i,bins+1)*width,plot=FALSE)$counts/width
    rect(xleft=i*width,ybottom=0,xright= (i+1)*width, ytop=height,col="grey70",border="grey40")
  }
  
  
  #plot individual intensities if requested
  if ("points" %in% plot){ 
  segments(x0=edges[1:(length(edges)-1)], y0=intensity[1:length(intensity)], 
           x1=edges[2:length(edges)], y1=intensity[1:length(intensity)],
           col = "cornflowerblue", lwd = 3)
  }
  

  #plot constant blocks
  if ("blocks" %in% plot){
    indices <- which(x$type=="constant")
    segments(x0=edges[x$left[indices]], y0=sapply( x$params , "[[" , 1 )[indices], 
             x1=edges[x$right[indices]+1], y1=sapply( x$params , "[[" , 1 )[indices],
             col = "red", lwd = 3)
    
    indices <- which(x$type=="linear")
    a=sapply( x$params , "[[" , "a" )[indices]
    b=sapply( x$params , "[[" , "b" )[indices]
    lefts=edges[x$left[indices]]
    rights=edges[x$right[indices]+1]
    segments(x0=lefts, y0=b, 
             x1=rights, y1=a*(rights-lefts)+b,
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
      floor((length(x$pruned)*100)/length(x$data)),'% of the points were pruned\n\n')
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
  lastparams    <- list()
  
  unpruned = NULL
  endobj = rep(0,n)
  optint = rep(0,n)
  
  
  #begin looping through each point
  for (i in 1:n){
    
    
    unpruned          <- c(unpruned,i)
    
    ##### constant blocks ##########
    
    #changeA[unpruned] <- changeA[unpruned] + A[i]
    #changeN[unpruned] <- changeN[unpruned] + N[i]
    #endobj[unpruned]  <- constantobjective(changeN[unpruned],changeA[unpruned],c)
    #optint[unpruned]  <- opt[unpruned] + endobj[unpruned]
    #opt[i+1]          <- max(optint[unpruned])
    #last              <- which.max(optint[unpruned])
    #lastchange[i]     <- unpruned[last]
    #lasttype[i]       <- "constant"
    #lastparams[i]     <- list("b"= sum(N[lastchange[i]:i])/sum(A[lastchange[i]:i]))

    ################################
    

    ##### linear blocks ############
    linblocks=list()
    for(j in unpruned){
      linblocks[[j]]  <- linfit(data[j:i]-max(data[j-1],0),c)
      optint[j]       <- opt[j] + linblocks[[j]][["cost"]]
    }

    ################################
    opt[i+1]         <- max(optint[unpruned])
    last            <- which.max(optint[unpruned])
    lastchange[i]     <- unpruned[last]
    lastparams[[i]]   <- linblocks[[unpruned[last]]]
    lasttype[i]       <- "linear"
    unpruned          <- unpruned[((optint[unpruned]+c-opt[i+1])>0)];

    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  return(getStep(BBdata,n))
}



#this function returns a BB object for the optimal blocks at any given cell location "n" where 1 < n < N
getStep <- function(BBobj,n){
  
  lasttype      <- BBobj$lasttype[1:n]
  lastchange    <- BBobj$lastchange[1:n]
  lastparams    <- BBobj$lastparams[1:n]
  
  lastchange<-lastchange-1
  
  left          <- vector()
  params        <- list()
  type          <- vector()
  type[1]       <- lasttype[n]
  left[1]       <- lastchange[n]
  params[[1]]   <- lastparams[[n]]
  i=1
  print(left)
  while (left[i] > 0) {
    left[i+1]      <- lastchange[left[i]]
    type[i+1]      <- lasttype[left[i]]
    params[[i+1]]  <- lastparams[[ left[i] ]]
    
    i <- i+1
  }
  left      <- rev(left)+1
  type      <- rev(type)
  params <- rev(params)
  
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
                   "params"       = params,
                   "opt"          = BBobj$opt[1:n],
                   "lastchange"   = lastchange+1,
                   "lastparams"   = lastparams,
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




#create random data
set.seed(1)

equation <- expression(-20*x+150) #generate some data with intensity function y=2x^2+5
equation2 <- expression(20*x+30) #generate some data with intensity function y=2x^2+5
equation3 <- expression(110) #generate some data with intensity function y=2x^2+5

x <-NHPois(time_length=6,equation=equation)
x2 <-NHPois(time_length=4,equation=equation2)+x[length(x)]
x3 <-NHPois(time_length=5,equation=equation3)+x2[length(x2)]
x=c(x,x2,x3)

N <- rep(1,length(x))#let all cells be of size 1


#run BB algorithm
start.time <- Sys.time()
model <- optinterval(x, N, 4)
print(Sys.time() - start.time)


#check summary information and plot
summary(model)
plot(model,plot=c("hist","blocks"),bins=100,ylim=c(0,250))
lines(x[x<6],(eval(equation))[x<6],lw=2,col="grey2") #plot the real background noise
lines(x[x>6 &x<10],(eval(equation2))[x>6& x<10]-120,lw=2,col="grey2") #plot the signal
lines(x[x>10],rep(eval(equation3),sum(x>10)),lw=2,col="grey2") #plot the signal
