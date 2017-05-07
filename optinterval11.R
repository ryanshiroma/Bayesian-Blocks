library("cairoDevice")
#####constant fitter
# constant fitter
constantobjective = function(N,A){
  value = N*(log(N/A))-N
  value[is.nan(value)]<-0
  return(value)
}

##############################exponential fitter
powfunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  tn=t[length(t)]
  eb=exp(1)^b
  sum(log(a*t^eb+c))-((a*tn^(eb+1))/(eb+1)+c*tn)
}

powgunc= function(x,t){
  a=x[1]
  b=x[2]
  c=x[3]
  eb=exp(1)^b
  atbc=a*t^(eb)+c
  teb=t^eb
  
  tn=t[length(t)]
  logt=log(t)
  fa=sum(teb/atbc)-((tn^(eb+1))/(eb+1))
  fb=sum((a*eb*t^(eb+1)*(eb*logt+logt-1))/(eb+1)^2)-(a*eb*tn^(eb+1)*((eb+1)*log(tn)-1))/(eb+1)^2
  fc=sum(1/atbc)-tn
  return(c(fa,fb,fc))
}

powfit <- function(N,x){
  if(length(x)<10){
    return(list("a"=0,"b"=0,"c"=0,"cost"=-Inf))
  }
  x=rep(x,N)
  fit <- optim(par=c(1,0,length(x)/x[length(x)]),fn=powfunc,gr=powgunc,control = list(fnscale = -1),t=x)
  pars=fit$par
  a=pars[[1]]
  b=exp(1)^pars[[2]]
  c=pars[[3]]
  
  if(sum((a*c(0,x)^b+c)<0)>0){
    return(list("a"=0,"b"=0,"c"=0,"cost"=-Inf))
  }
  cost=fit$value
  return(list("a"=a,"b"=b,"c"=c,"cost"=cost))
}

##############################exponential fitter
func= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  f=a*exp(1)^(b*xk)+c
  total=sum(1/f)
  sum((f-Yk)^2/(f*total))
}

gunc= function(x,Yk,xk){
  a=x[1]
  b=x[2]
  c=x[3]
  fa=exp(1)^(b*xk)
  fb=a*xk*exp(1)^(b*xk)
  fc=1
  f=a*exp(1)^(b*xk)+c
  f2=f^2
  f3=f2*f
  y2=Yk^2
  total=sum(1/f)
  part1=total*(f2-y2)
  part2=(f3-2*y*f2+y2*f)
  denom=f2*total^2
  da=fa*part1-sum(fa/f2)*part2
  db=fb*part1-sum(fb/f2)*part2
  dc=fc*part1-sum(fc/f2)*part2
  return(c(da,db,dc)/denom)
}

expfit <- function(N,x){
  if(length(x)<20){
    return(list("a"=0,"b"=0,"c"=0,"cost"=-Inf))
  }
  x=rep(x,N)
  bins=as.integer(max(2,min(length(x)/10,20)))
  
  # number of breaks needs to depend on the time length to give accurate estimates
  histogram = hist(x, breaks =seq(0,x[length(x)],length.out=bins+1), plot = FALSE)
  L = x[length(x)]         # total interval length
  N = length(histogram$counts)      # number of blocks on this interval
  binwidth=L/N
  k = 1:N                           # k'th block index (1<=k<=N)
  xk = (k - 1/2)*binwidth            # midpoints of blocks (from page 367)
  Yk = histogram$counts/binwidth   
  fit <- optim(par=c(1,0,10),fn=func,gr=gunc,Yk=Yk,xk=xk)
  pars=fit$par
  a=pars[[1]]
  b=pars[[2]]
  c=pars[[3]]
  if(sum((a*exp(1)^(b*x)+c)<0)>0){
    return(list("a"=0,"b"=0,"c"=0,"cost"=-Inf))
  }
  cost=sum(log(a*exp(1)^(b*x)+c))-((a/b)*exp(1)^(b*x[length(x)])+c*x[length(x)]-a/b)
  return(list("a"=a,"b"=b,"c"=c,"cost"=cost))
}

#### iwls linear fitter ###

iwls = function(N,t) {
  x=rep(t,N)
  bins=as.integer(max(length(x)/7,2))
  # number of breaks needs to depend on the time length to give accurate estimates
  histogram = hist(x, breaks =seq(0,x[length(x)],length.out=bins+1), plot = FALSE)
  L = max(histogram$breaks)         # total interval length
  N = length(histogram$counts)      # number of blocks on this interval
  binwidth=L/N
  k = 1:N                           # k'th block index (1<=k<=N)
  xk = (k - 1/2)*binwidth            # midpoints of blocks (from page 367)
  Yk = histogram$counts/binwidth            # event counts in bins
  epsilon=1
  est.weights=rep(1,N)
  M=x[length(x)]
  S=M^2/2
  old.cost=-Inf
  while(epsilon>1e-4){
    m = lm(Yk ~ xk,weights=est.weights)
    est.weights = 1/m$fitted.values
    coef=m$coefficients
    new.cost=(sum(log(coef[[2]]*x+coef[[1]]))-coef[[2]]*S-coef[[1]]*M)
    if (is.nan(new.cost)){
      return(list("a"=0,"b"=0,"cost"=-Inf))
    }
    epsilon = new.cost-old.cost
    old.cost=new.cost
  }
  
  list("a"=coef[2],"b"=coef[1],"cost"=new.cost)
}

#######linear fitter

linfit <-function(N,t){
  
  if (length(t)==0){ # if no data is recieved, set cost to 0
    return(list("a"=0,"b"=0,"cost"=-Inf))
  }  
  if (length(t)==1){ #if only one point is supplied, treat as constant intensity
    return(list("a"=0,"b"=N/t[1],"cost"=log(N/t[1])))
  }
  
  t=rep(t,N)
  
  epsilon  <- 1
  n        <- length(t)
  i        <- 1
  M        <- t[n]
  S        <- t[n]^2/2
  
  #start with some initial values for a and b
  coef     <- matrix(c(0,length(t)/M),nrow=2)
  
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
    
    #create the inverse of the hessian
    invhess <- 1/(fa*gb-fb*ga)*matrix(c(gb,-fb,-ga,fa),nrow=2,ncol=2,byrow=TRUE) 
    
    #run newtons method
    coef=coef-invhess%*%f 
    
    if(coef[2]<0 || (M*coef[1]+coef[2])<0){#if best line has a negative intensity, dont use linear
      return(list("a"=0,"b"=0,"cost"=-Inf))
    }
    #calculate the new cost
    new.cost <- sum(log(coef[1]*t+coef[2]))-coef[1]*S-coef[2]*M
    epsilon  <- new.cost-old.cost
    
    i <- i+1
  }
  list("a"=coef[1],"b"=coef[2],"cost"=new.cost)
}

#S3 generic plot function 

"plot.BB" <- function(x,show=c("hist","blocks"),binwidth=NULL,bins=NULL,ylim=NULL,xlim=NULL,xact=NULL,main="Bayesian Blocks",xlab="Time",ylab="Intensity") {
  
  data       <- x$data
  n          <- length(data)
  intensity  <- round(x$N/x$A,4)
  legend     <- vector()
  legend.col <- vector()
  
  #check if xlim is supplied
  if (!is.null(xlim)){
    lowerbound=xlim[1]
    upperbound=xlim[2]
  }
  else{
    lowerbound <- data[1]-(data[2]-data[1])/2
    upperbound <- data[n]+(data[n]-data[n-1])/2  
  }
  
  #check if ylim is supplied
  if(is.null(ylim)){
    ylim=c(0,max(intensity)*0.75)
  }
  
  #check if the hist bins are supplied
  if (is.null(binwidth) & is.null(bins)){
    binwidth=data[n]/100
    bins=100
  }
  else if (is.null(bins)){
    bins   <- round(data[n]/binwidth,0)
  }
  else{
    binwidth=data[n]/bins
  }
  
  
  #initialize plot
  plot(NA, ylim = ylim, xlim = c(lowerbound,upperbound),xaxt=xact,bty="n",
       main=main,xlab=xlab,ylab=ylab)  
  grid(NA,NULL)#add grid for y-axis
  
  
  
  #plot histogram if requested
  if ("hist" %in% show){
    histdata <- rep(data,x$N)
    end      <- histdata[length(histdata)]
    i        <- 0:(end/binwidth)
    
    height   <- hist(histdata, breaks=c(i,bins+1)*binwidth,plot=FALSE)$counts/binwidth
    rect(xleft=i*binwidth,ybottom=0,xright= (i+1)*binwidth, ytop=height,col="grey70",border="grey40")
    legend=c(legend,"Binned Data")
    legend.col=c(legend.col,"grey70")
  }
  
  
  #plot individual intensities if requested
  if ("points" %in% show){ 
    segments(x0=data[1:(length(data)-1)], y0=intensity[1:length(intensity)], 
             x1=data[2:length(data)], y1=intensity[1:length(intensity)],
             col = "cornflowerblue", lwd = 3)
    legend=c(legend,"Data Points")
    legend.col=c(legend.col,"cornflowerblue")
  }
  
  
  #plot blocks if requested
  if ("blocks" %in% show){
    #plot constant blocks
    indices <- which(x$type=="constant")
    if (length(indices)>0){
      b=sapply( x$params [indices], "[[" , "b" )
      segments(x0=data[x$left[indices]], y0=b, 
               x1=data[x$right[indices]], y1=b,
               col = "red", lwd = 3)
    }
    #plot linear blocks
    indices <- which(x$type=="linear")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      segments(x0=lefts, y0=b, 
               x1=rights, y1=a*(rights-lefts)+b,
               col = "red", lwd = 3)
    }
    
    indices <- which(x$type=="exponential")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        lines(xv,a[i]*exp(1)^(b[i]*(xv-lefts[i]))+c[i],col="red",lwd=3)
      }
    }
    
    indices <- which(x$type=="power")
    if(length(indices)>0){
      a=sapply( x$params[indices], "[[" , "a" )
      b=sapply( x$params[indices], "[[" , "b" )
      c=sapply( x$params[indices], "[[" , "c" )
      lefts=data[x$left[indices]]
      rights=data[x$right[indices]]
      for (i in 1:length(a)){
        xv=seq(lefts[i],rights[i],length.out=100)
        lines(xv,a[i]*(xv-lefts[i])^b[i]+c[i],col="red",lwd=3)
      }
    }
    
    
    legend=c(legend,"Blocks")
    legend.col=c(legend.col,"red")
    
  }
  
  #plot legend
  legend(x=x$data[1],y=ylim[2]*0.98,legend= legend, lty=rep(1,length(legend)), 
         lwd=rep(3,length(legend)), col=legend.col) 
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

optinterval = function(data,N,c,type=c("constant"),pvalue=0.05,verbose=FALSE){
  start.time <- Sys.time()
  
  chi.lin    <- qchisq(1-pvalue,df=1)/2
  chi.pow    <- chi.lin*2
  if("power" %in% type){
    maxpen   <- chi.pow +c  
  }else if("linear" %in% type){
    maxpen   <- chi.lin+c
  }else{
    maxpen   <- c
  }
  
  n          <- length(N)
  percent    <- floor(n/100)
  lowerbound <- data[1]-(data[2]-data[1])/2
  upperbound <- data[n]+(data[n]-data[n-1])/2
  xx         <- c(lowerbound,(data[1:(length(data)-1)] + data[2:(length(data))])/2,upperbound) # voronoi cell vertices
  A          <- diff(xx) # length of each voronoi cell
  
  opt           <- rep(0,n+1)
  lastchange    <- rep(1,n)
  changeA       <- rep(0,n) 
  changeN       <- rep(0,n)
  optint1       <- rep(-Inf,n)
  optint2       <- rep(-Inf,n)
  optint3       <- rep(-Inf,n)
  lasttype      <- rep("None",n)
  lastparams    <- list()
  
  unpruned <- NULL
  endobj   <- rep(0,n)
  optint   <- rep(0,n)
  
  
  #begin looping through each point
  for (i in 1:n){
    
    unpruned          <- c(unpruned,i)
    
    ##### constant blocks ##########
    if ("constant" %in% type){
      changeA[unpruned]  <- changeA[unpruned] + A[i]
      changeN[unpruned]  <- changeN[unpruned] + N[i]
      optint1[unpruned]  <- opt[unpruned] + constantobjective(changeN[unpruned],changeA[unpruned])-c
      last1              <- which.max(optint1[unpruned])
    }
    ################################
    
    ##### linear blocks ############
    if ("linear" %in% type){
      linblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        linblocks[[j]]   <- linfit(N[j:i],x)
        optint2[j]       <- opt[j] + linblocks[[j]][["cost"]]-c-chi.lin
      }
      last2              <- which.max(optint2[unpruned])
    }
    
    # ##### exp blocks ############
    # if ("exponential" %in% type){
    #   expblocks=list()
    #   for(j in unpruned){
    #     x=data[j:i]-max(data[j-1],0)
    #     expblocks[[j]]   <- expfit(N[j:i],x)
    #     optint3[j]       <- opt[j] + expblocks[[j]][["cost"]]-c-chi.exp
    #   }
    #   last3              <- which.max(optint3[unpruned])
    # }
    #   ################################
    
    ##### exp blocks ############
    if ("power" %in% type){
      powblocks=list()
      for(j in unpruned){
        x=data[j:i]-max(data[j-1],0)
        powblocks[[j]]   <- powfit(N[j:i],x)
        optint3[j]       <- opt[j] + powblocks[[j]][["cost"]]-c-chi.pow
      }
      last3              <- which.max(optint3[unpruned])
    }
    ################################  
    
    ################################
    
    bestshape<-which.max(c(max(optint1),max(optint2),max(optint3)))
    
    if(bestshape==1){ #constant block is best
      lasttype[i]       <- "constant"
      lastchange[i]     <- unpruned[last1]
      lastparams[[i]]   <- list("b"= sum(N[lastchange[i]:i])/sum(A[lastchange[i]:i]))
    }
    
    else if (bestshape==2){#linear block is best
      lasttype[i]       <- "linear"
      lastchange[i]     <- unpruned[last2]
      lastparams[[i]]   <- linblocks[[unpruned[last2]]]
    }
    
    else if(bestshape==3){ #exp block is best
      lasttype[i]       <- "power"
      lastchange[i]     <- unpruned[last3]
      lastparams[[i]]   <- powblocks[[unpruned[last3]]]
    }
    
    
    
    
    if(lastchange[i]<max(lastchange[i-1],0)){
      cat("+")
    }else if(lastchange[i]!=max(lastchange[i-1],0)){
      cat("+")
    }
    
    optint=apply(matrix(c(optint1,optint2,optint3),ncol=3), 1, max)[1:i]
    opt[i+1]          <- max(optint[unpruned])
    
    
    unpruned          <- unpruned[((optint[unpruned]+maxpen-opt[i+1])>0)]
    
    if((verbose==TRUE) && (i %% percent==0)){#print out the progress of the algorithm
      cat("\n",round(100*i/n,0),"% of points completed",round(100*(i-length(unpruned))/i,0),"% pruned ")
    }
    
  }
  BBdata   <- list("data"         = data, 
                   "N"            = N,
                   "A"            = A,
                   "opt"          = opt,
                   "lastchange"   = lastchange,
                   "lastparams"   = lastparams,
                   "lasttype"     = lasttype,
                   "pruned"       = setdiff(1:i,unpruned))
  
  
  
  print(Sys.time()-start.time)
  model <- getStep(BBdata,n)
  summary(model)
  return(model)
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
  while (left[i] > 0) {
    left[i+1]      <- lastchange[left[i]]
    type[i+1]      <- lasttype[left[i]]
    params[[i+1]]  <- lastparams[[ left[i] ]]
    
    i <- i+1
  }
  left      <- rev(left)+1
  type      <- rev(type)
  params    <- rev(params)
  
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
set.seed(0)

equation1 <- expression(50*x^4+100)
equation2 <- expression(-200*x+900) 
equation3 <- expression(300) 

x1 <- NHPois(time_length =2,equation = equation1)
x2 <- NHPois(time_length = 3,equation = equation2)+x1[length(x1)]
x3 <- NHPois(time_length = 3,equation = equation3)+x2[length(x2)]

x=c(x1,x2,x3)
x=round(x,3)
x=unique(x)
N <- rep(1,length(x))#let all cells be of size 1

#run model
model <- optinterval(x, N, 60,type=c("constant","linear","power"),pvalue=0.05,verbose=TRUE)

plot(model,show=c("points","blocks"),bins=70,xlim=c(0,x[length(x)]),ylim=c(0,800),main="Bayesian Blocks")
lines(x,eval(equation1))
lines(x,eval(equation1),col="red",lwd=2)

save(file = "model.Rdata",model)

load(file="model.Rdata")
#check summary information and plot
summary(model)
plot(model,show=c("points","hist","blocks"),bins=70,xlim=c(0,x[length(x)]),ylim=c(0,2000),main="Bayesian Blocks")
expfunc(c(10,5,20),x)
expfunc(c(110,-266,428),x)

#plot the segmentation in steps
for (i in 1:20){
  #png(filename=paste0("~/ne3/",sprintf("%04d",i),".png"),800,640,type='cairo')
  plot(getStep(model,i*(length(x)/20)),show=c("hist","blocks"),binwidth=0.2,ylim=c(0,600),xlim=c(0,x[length(x)]))
  Sys.sleep(0.5)
  #dev.off()
}
plot(model,show=c("hist","blocks"),binwidth=0.2,ylim=c(0,400),xlim=c(0,x[length(x)]))

dev.off()

model$params[[2]]

