

cov.mariah <- function(h, theta){
  h[h>0] <- theta * log(1 + h[h>0]/theta)/ h[h>0]
  h[h==0] <- 1
  h
}

cov.exp <- function(h, theta)
  exp(- h/theta)

cov.spheric <- function(h, theta){
  h[h<=theta] <- 1 - 3/2 * h[h<=theta]/theta + h[h<=theta]^3 / theta^3 / 2
  h[h>theta] <- 0
  h}

cov.linear <- function(h, theta){
  h[h<=theta] <- 1 - h[h<=theta]/theta
  h[h>theta] <- 0
  h}

var.frac <- function(h,s,alpha){
  h^alpha/s}

Vario2EC <- function(x) {2 * pnorm( sqrt(x) / 2 )}

EC2Vario <- function(theta) { (2*qnorm(theta/2))^2 }



vario <- function(par,
                  riv.dist = RiverDis,
                  euc.coords = CatCntrWt,
                  riv.weights = Weight
)
{
  theta <- par[1]
  alpha <- par[2]
  s <- par[3]
  riv.lmb <- par[4]
  euc.lmb <- par[5]        
  beta <- par[6]
  c <- par[7]
  
  rot.mat <- cbind( c(cos(beta), c*sin(beta)),c(-sin(beta), c*cos(beta)))
  euc.dist <- as.matrix(dist(euc.coords %*% t(rot.mat)))        
  
  riv.lmb * 2 * (1 - cov.linear(riv.dist, 500*theta)*riv.weights) +
    euc.lmb * var.frac(euc.dist, s, alpha)  ##  * 500
}


vario.M1 <- function(par,
                     riv.dist = RiverDis,
                     euc.coords = CatCntrWt,
                     riv.weights = Weight
)
{
  theta <- par[1]
  alpha <- par[2]
  s <- par[3]
  riv.lmb <- par[4]
  euc.lmb <- par[5]        
  beta <- par[6]
  c <- par[7]
  
  rot.mat <- cbind( c(cos(beta), c*sin(beta)),c(-sin(beta), c*cos(beta)))
  euc.dist <- as.matrix(dist(euc.coords %*% t(rot.mat)))        
  
  euc.lmb * var.frac(euc.dist,s,alpha)
}

########################### Declustering Function ################################
#################################################################################
#################################################################################

ObtainMultVarEvents <- function(TSs, U, Lag, Year, StNames, Plotting, mfrow)
{
  # Obatain supposedly independent extreme events in a single year
  
  # TSs: common measurements of all stations and all years 
  # U: probability input for the quantile function
  # Lag: window half-width (total length is 2*lag + 1)
  # Year: year (1960-2010)
  # StNames: station indices (1-31)

  CommonDate <- TSs$Date
  YearAll <- as.numeric(substr(CommonDate, 1, 4))
  Data <- data.frame(TSs[, -1]) # drop the first column (date)
  NoSt <- length(StNames)
  Quant  <- apply(Data, 2, quantile, probs = U) # 10%-percentile of measurements at all times within each station 
  names(Quant) <- NULL
  
  # compute empirical cdfs (ranks) of measurements at all times within each station
  Ecdf <- Data
  for (i in 1:NoSt)
  {
    Temp1 <- ecdf(Data[, i])
    Ecdf[,i] <- Temp1(Data[, i])
  }
  
  # extract obs and ecdfs in a specific year
  Ind <- which(YearAll == Year) 
  DataYear <- data.frame(Data[Ind, ]) 
  EcdfYear <- data.frame(Ecdf[Ind, ]) 
  n <- nrow(DataYear)
  row.names(DataYear) <- c(1:n)
  row.names(EcdfYear) <- c(1:n)
  ObtainMultVarEvents <- function(Data, Ecdf, Quant, m)
  {
    # X <- DataYear
    n <- nrow(Data)
    Q  <- matrix(rep(Quant, nrow(Data)), ncol = NoSt, byrow = T)
    Data.copy <- Data
    Events <- numeric() # flood events (tilde X's)
    Location <- numeric() # day of the measurement for each station for each flood event
    CenterEventLoc <- numeric() # center of each flood event window
    while (sum(Data.copy > Q) > 0)
    {
      IndTemp <- which(Data.copy >= Q, arr.ind=TRUE)
      if (length(IndTemp) > 0)
      {
        Ind <- IndTemp[which(Ecdf[IndTemp] == max(Ecdf[IndTemp]))][1] # day of the largest measurement among all stations within a year
        lower <- max(1, Ind - m) # lower end of the window
        upper <- min(n, Ind + m) # upper end of the window
        CenterEventLoc <- c(CenterEventLoc, Ind)
        
        Data.window <- data.frame(Data.copy[lower:upper, ]) # measurements within the window
        Day.window <- lower:upper # days within the window
        Event <- apply(Data.window, 2, max) # each station's largest measurements within the window

        # extract the day corresonding to the largest measurement 
        # for each station within the window
        Loc <- numeric()
        for (j in 1:NoSt)
        {
          Temp2 <- which(Data.window[, j] == Event[j])[1]
          Loc[j] <- as.numeric(Day.window[Temp2])
        }
        
        # delete measurements in the window
        # IMPORTANT: NOT Data.copy[lower:upper, ] <- 0 due to NON-OVERLAPPING requirement
        Data.copy[max(lower - m, 0):min(upper + m, n), ] <- 0   
        
        Events <- rbind(Events, Event)
        Location <- rbind(Location, Loc)
      }
    }

    row.names(Events) <- NULL
    row.names(Location) <- NULL
    names(CenterEventLoc)<- NULL 
    EventYear <- list()
    EventYear$Events <- Events
    EventYear$Location <- Location 
    EventYear$CenterEventLoc <- CenterEventLoc
    return(EventYear)
  }
  
  AllEventsTemp <- numeric(NoSt)
  Events <- ObtainMultVarEvents(Data = DataYear, Ecdf = EcdfYear, Quant = Quant, m = Lag)
  EventsTemp <- Events$Events
  
  Result <- list()
  Result$AllEvMat <- t(Events$Events)
  Result$Location <- t(Events$Location)
  Result$CenterEventLoc <- Events$CenterEventLoc
  
  return(Result)
}

ExceedancesNormSum <- function(TSs,TSsOrig,U,Lag,Year,StNames,Plotting,mfrow)
{
  DateAll <-  TSs$Date
  TSNoDat <- TSs[,-1]
  TSOrigNoDat <- TSsOrig[,-1]
  
  row.names(TSNoDat) <- c(1:dim(TSNoDat)[1])
  row.names(TSOrigNoDat) <- c(1:dim(TSOrigNoDat)[1])
  
  TsSum <- apply(TSNoDat,1,sum)
  Quant <- quantile(TsSum ,U)
  
  Ind <- which((as.numeric(format(as.Date(DateAll , "%Y-%m-%d"),"%Y"))== Year))
  TSsYear <- TSNoDat [Ind,]
  TSsYearOrig <- TSOrigNoDat[Ind,]
  
  row.names(TSsYear) <- c(1:dim(TSsYear)[1])
  row.names(TSsYearOrig) <- c(1:dim(TSsYearOrig)[1])
  
  DateYear <- DateAll[Ind]
  NoSt <- dim(TSNoDat)[2]
  
  TsSumYear <- apply(TSsYear ,1,sum) 
  Exceed <- (TsSumYear > Quant)*TsSumYear 
  X <- Exceed
  NoRow <- length(X)
  k <- 1
  YearEvents <- list()
  while (sum(X) > 0)
  {
    EventAndLoc <- list()
    m <- which(X == max(X))
    m1 <- max(1,m[1]-Lag)	
    m2 <- min(NoRow,(m[1]+Lag))
    EventRun <- TSsYear[m1:m2,]
    Event <- apply(EventRun,2,max)
    Location <- numeric()
    for (NoLoc in 1:length(Event))
    {
      Location[NoLoc] <- m1 + which(EventRun[,NoLoc]==Event[NoLoc],arr.ind = TRUE)[1]-1
    }
    m11 <- max(1,m1-Lag)
    m22 <- min(NoRow,m2+Lag)
    X[m11:m22] <- 0
    EventAndLoc$Event <- Event
    EventAndLoc$Location <- Location
    YearEvents[[k]] <- EventAndLoc 
    k <- k+1
  }
  if (length(YearEvents)>0 & Plotting==1)
  {
    par(mfrow=mfrow, oma = c(1,1,1,1) + 0.1,
        mar = c(1,1,1,1) + 0.1)
    xlim <- c(1,NoRow)
    for(i in 1:NoSt)
    {
      ylim <- c(min(TSsYearOrig[,i])-1,max(TSsYearOrig[,i])+10)
      plot(TSsYearOrig[,i],xlim=xlim,ylim=ylim,xlab="",ylab="",axes=F,type="l",lty=2)
      mtext(StNames[i], side=4)
      #abline(h=Quant[i])
      axis(2)
      box() 
      for (j in 1:length(YearEvents))
      {
        par(new=T)
        plot(YearEvents[[j]]$Location[i],TSsYearOrig[YearEvents[[j]]$Location[i],i],col="red",xlim=xlim,ylim=ylim,xlab="",ylab="",axes = FALSE)
      }
    }
    axis(side=1, at = c(1,NoRow),
         labels = DateYear[c(1,NoRow)])
    
    PlotName <- paste("Events",as.character(Year),".pdf",sep="")
    SavingPlotAdress <- paste(getwd(),"Bavaria/Results/EventsPlots",PlotName,sep="/")
    dev.copy(pdf,SavingPlotAdress)
    dev.off()
  }
  n1 <- length(YearEvents)
  Result <- list()
  AllEvMat <- matrix(0,NoSt,n1)
  AllEvLoc <- matrix(0,NoSt,n1)
  
  if (n1 >0)
  {
    for (j in 1:n1)
    {
      AllEvMat[,j]<- YearEvents[[j]]$Event
      AllEvLoc[,j] <- YearEvents[[j]]$Location
    }
    Result$AllEvMat <- AllEvMat
    Result$AllEvLoc<- AllEvLoc
    #Result$Date <- DateYear
  }
  
  return(Result)
  
}

################## Obtaing Theta #################
library(evd)
library(plotrix)

GetTheta <- function(StsChos,Q)
{
  
  load("Data/AllPairEvents.RData")
  
  Names <- names(AllPairEvents)
  NoSts <- length(StsChos)
  Theta <- matrix(NA,NoSts,NoSts)
  for (i in 1:(NoSts-1))
  {
    for (j in (i+1):NoSts)
    {
      A1 <- paste(as.character(StsChos[i]),as.character(StsChos[j]),sep=",")
      A2 <- paste(as.character(StsChos[j]),as.character(StsChos[i]),sep=",")
      Ind <- which((Names==A1 | Names==A2  ))
      X <- AllPairEvents[[Ind]]
      Chi <- try(chiplot(X,which=1,nq = 2,qlim =c(0.6,Q),conf = 0),silent=TRUE)
      if (length(Chi)>2)
      {
        Theta[i,j] <- 2-Chi$chi[2,2]
        Theta[j,i] <- Theta[i,j]
      }
    }
  }
  for (i in 1:NoSts)
    Theta[i,i] <- 1
  
  return(Theta)
}


###################Transform the margisn to pareto ########################
###########################################################################

Margin2Pareto <- function(Data, quant = 0.95, clustering = FALSE, empirical = FALSE)
{
  require("evd")
  len <- length(Data)
  if (empirical)
  {
    trans.data <- double(len)
    trans.data <- 1/(1 - rank(Data) / (len+1))
  } else{            
    threshold <- quantile(Data, quant,na.rm=T)
    result <- fpot(x = Data, threshold = threshold,
                   model = "gpd", cmax = clustering, r = 2, std.err=FALSE)
    
    shape <- result$param[2]
    scale <- result$param[1]
    
    trans.data <- double(len)
    trans.data[Data <= threshold] <-
      1/(1 - rank(Data)[Data <= threshold] / (len+1))
    
    trans.data[Data > threshold] <- pmax(
      (1 + shape * (Data[Data > threshold]- threshold) / scale),
      0)^(1 / shape) / (1 - quant)
  }
  return(trans.data)
}


#------------------------------------------------------------------------------------------------------------------------------
# Fitting of Brown-Resnick process on river network using SPECTRAL DENSITY
#------------------------------------------------------------------------------------------------------------------------------
pareto_BR_River <- function(
    data,
    riv.dist,
    euc.coord,
    is.connected,
    riv.weights,
    u,
    init,
    model = 3,
    maxit = 100,
    method = method,
    verbose = FALSE)
{
  if (model == 1) {
    fixed = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)    # fix lambda_riv 
  } else if (model == 2) {
    fixed = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)    # fix lambda_riv 
  } else if (model == 3) {
    fixed = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
  } else {
    fixed = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)     # fix beta, c 
  }
  
  ti <- proc.time()
  d <- ncol(data)
  
  # log spectral density
  logdV <- function(x,par){
    if (is.vector(x)){d <- length(x)}
    if (is.matrix(x)){d <- ncol(x)}	
    i <- 1
    varioHalf <- 1/2 * vario(par,
                             riv.dist = riv.dist,
                             euc.coords = euc.coord,
                             riv.weights = riv.weights
    )
    vario0 <- matrix(rep(varioHalf[-i,i],d-1),ncol = d-1)
    S <- vario0 + t(vario0) - varioHalf[-i,-i]
    cholS <- chol(S)
    Sm1 <- chol2inv(cholS)
    logdetS <- 2*sum(log(diag(cholS))) 
    
    if (is.vector(x)){
      y <- (log(x/x[i])+ varioHalf[,i])[-i]
      logdv <- - sum(log(x)) - log(x[i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * t(y)%*%Sm1%*%y
    }
    
    if (is.matrix(x)){
      y <- (t(t(log(x/x[,i])) + varioHalf[,i]))[,-i]
      logdv <- - apply(log(x),1,sum) - log(x[,i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * diag(y%*%Sm1%*%t(y))
    }
    
    return(logdv)
  }
  
  # project on the 1-sphere above the threshold
  project <- function(x,u){
    y <- x[which(rowSums(x) > u),]
    y <- y / rowSums(y)            
  }
  
  data.u <- project(data,u)
  
  # negative log likelihood function
  nllik <- function(par, model = 3){
    theta <- par[1]
    alpha <- par[2]
    s <- par[3]
    riv.lmb <- par[4]
    euc.lmb <- par[5]
    beta <- par[6]
    c <- par[7]
    
    if (model == 1) {

      if (theta <= 0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        y <- sum(logdV(x=data.u, par=par))
        return(-y)
      }
      
    } else if (model == 2) {

      if (theta <= 0 |euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        y <- sum(logdV(x=data.u, par=par))
        return(-y)
      }
      
    } else if (model == 3) {
      
      if (theta <= 0 | riv.lmb <0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        y <- sum(logdV(x=data.u, par=par))
        return(-y)
      }
      
    } else {
      
      if (theta <= 0 | riv.lmb <0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        y <- sum(logdV(x=data.u, par=par))
        return(-y)
      }
    }
    
  }
  
  init2 <- init[!fixed]
  nllik2 <- function(x){y <- numeric(length(init)); y[!fixed] <- x;
  y[fixed] <- init[fixed];
  return(nllik(y, model))}
  
  # optimize likelihood
  opt <- optim(
    init2, 
    nllik2, 
    hessian = TRUE,
    control = list(maxit = maxit),
    method = method
  )
  
  z <- list()
  z$convergence <- opt$convergence
  z$par[fixed] <- init[fixed]
  z$par[!fixed] <- opt$par
  z$nllik <- opt$value
  z$hessian <- opt$hessian
  z$time <- (proc.time() - ti)[3]
  
  if (verbose) {
    print(z$nllik)
    print(z$par)
  }
 
  return(z)
}


#------------------------------------------------------------------------------------------------------------------------------
# Fitting of Brown-Resnick process on river network using CENSORED ESTIMATION
#------------------------------------------------------------------------------------------------------------------------------

pareto_BR_River_cen <- function(data,
                                riv.dist,
                                euc.coord,
                                riv.weights,
                                u,
                                init,
                                # fixed=c(FALSE,FALSE,FALSE,FALSE,FALSE),
                                model = 3,
                                maxit = 100,
                                method = "BFGS",
                                seed=123456,
                                verbose = FALSE){
  require("mvtnorm")
  ti <- proc.time()
	d <- ncol(data)
	if (length(u)==1){u <- rep(u,d)}
  
	if (model == 1) {
	  fixed = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)    # fix lambda_riv 
	} else if (model == 2) {
	  fixed = c(FALSE, FALSE, TRUE, TRUE, FALSE, FALSE, FALSE)    # fix lambda_riv 
	} else if (model == 3) {
	  fixed = c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE)
	} else {
	  fixed = c(FALSE, FALSE, TRUE, FALSE, FALSE, TRUE, TRUE)     # fix beta, c 
	}
	
	# exponent measure and its derivatives
	V <- function(x, par){
	          if (verbose) {
	            print(par)
	          }
            
            d <- length(x)
     
            varioHalf <- 1/2 * vario(par,
                                     riv.dist = riv.dist,
                                     euc.coords = euc.coord ,
                                     riv.weights = riv.weights
                                     )            

            f1 <- function(i,x){                
                vario0 <- matrix(rep(varioHalf[-i,i],d-1),ncol = d-1)
                S <- vario0 + t(vario0) - varioHalf[-i,-i]
                return(1/x[i]*pmvnorm(upper=(log(x/x[i])+varioHalf[,i])[-i],mean=rep(0,d-1),sigma= S)[1])
            }                      
            
            return(sum(apply(cbind(1:d),1,f1,x=x)))
	}


	logdV <- function(x,par){
		if (is.vector(x)){d <- length(x)}
		if (is.matrix(x)){d <- ncol(x)}	

                i <- 1
                
                varioHalf <- 1/2 * vario(par,
                                         riv.dist = riv.dist,
                                         euc.coords = euc.coord ,
                                         riv.weights = riv.weights
                                         )     
                
                vario0 <- matrix(rep(varioHalf[-i,i],d-1),ncol = d-1)
                S <- vario0 + t(vario0) - varioHalf[-i,-i]
                cholS <- chol(S)
		Sm1 <- chol2inv(cholS)
		logdetS <- 2*sum(log(diag(cholS)))
                if (is.vector(x)){
                    y <- (log(x/x[i])+ varioHalf[,i])[-i]
                    logdv <- - sum(log(x)) - log(x[i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * t(y)%*%Sm1%*%y
                }
                if (is.matrix(x)){
                    y <- (t(t(log(x/x[,i])) + varioHalf[,i]))[,-i]
                    logdv <- - apply(log(x),1,sum) - log(x[,i]) -((d-1)/2)*log(2*pi) -1/2*logdetS  - 1/2 * diag(y%*%Sm1%*%t(y))
                }
                
                return(logdv)                
          
	}

	logdVK <- function(x,K,par){
            
                d <- length(x)                
            
                k <- length(K)
                
                i <- min(K)
                idxK <- which(K == i)
                
                varioHalf <- 1/2 * vario(par,
                                         riv.dist = riv.dist,
                                         euc.coords = euc.coord ,
                                         riv.weights = riv.weights
                                         )   
                vario0 <- matrix(rep(varioHalf[,i],d),ncol = d)
                S <- vario0 + t(vario0) - varioHalf
                if(k>1){
                    SK <- S[K[-idxK],K[-idxK]]                
                    cholSK <- chol(SK)
                    SKm1 <- chol2inv(cholSK)                
                    logdetSK <- 2*sum(log(diag(cholSK)))
                    idxK <- which(K == i)                
                    yK <- (log(x[K]/x[i])+ varioHalf[K,i])[-idxK]
                    logdvK <- - sum(log(x[K])) - log(x[i]) -((k-1)/2)*log(2 *pi) - 1/2*logdetSK - 1/2 * t(yK)%*%SKm1%*%yK
                    SnK <- S[-K,-K]
                    SnKK <- S[-K,K[-idxK]]
                    SKnK <- t(SnKK)
                    muCondK <- -varioHalf[-K,i] + SnKK %*% SKm1 %*% yK
                    if(k < d-1)
                        SCondK <- SnK - SnKK %*% SKm1 %*% SKnK
                    if(k == d-1)
                        SCondK <- SnK - t(SnKK) %*% SKm1 %*% t(SKnK)
                    logdvnK <- log(pmvnorm(upper=c(log(x[-K]/x[i])-muCondK),sigma=SCondK)[1])
                    logdv <- logdvK + logdvnK
                }
                if(k==1){
                    logdvK <- - 2*log(x[i])
                    logdvnK <- log(pmvnorm(upper=c(log(x[-K]/x[i]) + varioHalf[-K,i]),sigma=S[-K,-K])[1])
                    logdv <- logdvK + logdvnK
                }                

                return(logdv)
                
         
	}
	
	# censor below the (multivariate) threshold
	censor <- function(x,u){
		f2 <- function(x,u){
			y <- c()
			y[x>u] <- x[x>u]
			y[x<=u] <- u[x<=u]
			return(y)
			}
		return(t(apply(x,1,f2,u)))	
	}

	data.u <- censor(data,u)

	r <- nrow(data.u)

	L <- apply(data.u>matrix(u,ncol=d,nrow=r,byrow=TRUE),1,which)
	I <- which(lapply(L,length)>0 & lapply(L,length)<d)
	J <- which(lapply(L,length)==d)

	# negative log likelihood function
	nllik <- function(par){
	  set.seed(seed)
    theta <- par[1]
    alpha <- par[2]
    s <- par[3]
    riv.lmb <- par[4]
    euc.lmb <- par[5]
    beta <- par[6]
    c <- par[7]
                
    if (model == 1) {

      if (theta <= 0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        
        if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
        else {y2 <- 0}
        
        y <- sum(y1)+sum(y2) + (r-(length(I)+length(J)))*log(1-V(u,par=par))
        if (verbose) {
          print(-y)
          print(par)
        }
        return(-y)
      }
    } else if (model == 2) {

      if (theta <= 0 |euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        
        if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
        else {y2 <- 0}
        
        y <- sum(y1)+sum(y2) + (r-(length(I)+length(J)))*log(1-V(u,par=par))
        if (verbose) {
          print(-y)
          print(par)
        }
        return(-y)
      }
      
    } else if (model == 3) {
      
      if (theta <= 0 | riv.lmb <0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        
        if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
        else {y2 <- 0}
        
        y <- sum(y1)+sum(y2) + (r-(length(I)+length(J)))*log(1-V(u,par=par))
        if (verbose) {
          print(-y)
          print(par)
        }
        return(-y)
      }
      
    } else {
      
      if (theta <= 0 | riv.lmb <0 | euc.lmb <0 | alpha > 2 |
          alpha <= 0 | s <= 0) # iterates are infeasible
      {
        return(10^50)
      }                
      else{	
        if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
        else {y1 <- 0}
        
        if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
        else {y2 <- 0}
        
        y <- sum(y1)+sum(y2) + (r-(length(I)+length(J)))*log(1-V(u,par=par))
        if (verbose) {
          print(-y)
          print(par)
        }
        
        return(-y)
      }
    }
                
# 		if (theta<=0 | riv.lmb <0 | euc.lmb <0 | alpha > 2 | 
# 		    alpha <= 0 | s <= 0 | beta > 3*pi/4 | beta < pi/4 | c < 0){
# 		  return(10^50) 
# 		
# 		} else {
# 	
# 		if (length(I)>0){y1 <- mapply(logdVK,x=as.list(data.frame(t(data.u)))[I],K=L[I],MoreArgs=list(par=par))}
# 		else {y1 <- 0}
# 
# 		if (length(J)>0){y2 <- logdV(x=data.u[J,],par=par)}
# 		else {y2 <- 0}
# 	
# 		y <- sum(y1)+sum(y2) + (r-(length(I)+length(J)))*log(1-V(u,par=par))
#                 print(-y)
#                 print(par)
#                 return(-y)
# 		}
	}

  init2 <- init[!fixed]
  nllik2 <- function(x){
    y <- numeric(length(init)); 
    y[!fixed] <- x;
    y[fixed] <- init[fixed];
    return(nllik(y))
  }
	
	# optimize likelihood
	opt <- optim(init2,
	             nllik2,
	             hessian = TRUE,
	             control = list(maxit = maxit), 
	             method = method)
	
	z <- list()
	z$convergence <- opt$convergence
  z$par[fixed] <- init[fixed]
	z$par[!fixed] <- opt$par
	z$nllik <- opt$value
	z$hessian <- opt$hessian
        z$time <- (proc.time()-ti)[3]
	
	return(z)
}



######################## Plotting Functions ##############################
##########################################################################


plotECF <- 
  function(ECemp,
           ECtheo,
           Dist = RivDistancesNew,
           is.con = IsConnectedNew,
           StsIdx = StsChos,
           which.plots = c(TRUE,TRUE,TRUE),
           which.labels = c(TRUE,TRUE,TRUE),
           PDF = FALSE,
           filename = "")
  {
    
    if(PDF) pdf(filename, width = sum(which.plots) * 7)
    par(cex = 1.3, cex.lab = 2.3, cex.axis = 2.3, cex.main = 1.5, pch = 19,
        mar = c(5,5,4,2) +.1)
    len <- length(StsChos)
    RivLabel <- matrix("NA", ncol = len, nrow = len)
    for (i in 1:(len-1))
      for (j in (i+1):len)
        RivLabel[i,j] <- paste(as.character(StsIdx[i]),as.character(StsIdx[j]),sep=",")
    
    upper.idx <- upper.tri(ECemp)
    con.idx <- as.logical(upper.idx * is.con)
    not.con.idx <- as.logical(upper.idx * !(is.con))
    par(mfrow=c(1,sum(which.plots)))
    ## Plot ECemp
    if(which.plots[1]){
      plot(Dist[not.con.idx],ECemp[not.con.idx], xlab = "Hydrological Distance (km)", ylab="Extremal Coefficient", ylim = c(1,2),xlim = c(0,max(Dist)), mgp=c(3.3, 1.3, 0))     
      points(Dist[con.idx],ECemp[con.idx],col="blue",pch=4,cex = 1.3)
      if(which.labels[1]) thigmophobe.labels(Dist[upper.idx], ECemp[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)      
    }
    ## Plot ECtheo
    if(which.plots[2]){
      plot(Dist[not.con.idx],ECtheo[not.con.idx], xlab = "Hydrological Distance (km)", ylab="Extremal Coefficient",ylim = c(1,2),xlim = c(0,max(Dist)), mgp=c(3.3, 1.3, 0)) 
      points(Dist[con.idx],ECtheo[con.idx],col="blue",pch=4,cex = 1.3)
      if(which.labels[2]) thigmophobe.labels(Dist[upper.idx], ECtheo[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)
    }
    ## Plot QQ
    if(which.plots[3]){
      plot(ECtheo[not.con.idx], ECemp[not.con.idx], xlim = c(1,2), ylim = c(1,2),ylab="Empirical", xlab="Fitted Model", mgp=c(3.3, 1.3, 0)) 
      abline(b = 1, a = 0)
      points(ECtheo[con.idx],ECemp[con.idx],col="blue",pch=4,cex = 1.3)
      if(which.labels[3]) thigmophobe.labels(ECtheo[upper.idx], ECemp[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)
    }
    if(PDF) dev.off()
  }

##################################################
plotEucVsHyd <- function(ECemp,
                         ECtheo,
                         DistEuc = EucDisChos,
                         DistCatch = CatchDistancesNew,
                         is.con = IsConnectedNew,
                         StsIdx = StsChos,
                         which.plots = c(TRUE,TRUE,FALSE),
                         which.labels = c(TRUE,TRUE,TRUE),
                         PDF = FALSE,
                         filename = "")
{
  
  if(PDF) pdf(filename, width = sum(which.plots) * 7)
  par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      mar = c(5,5,4,2) +.1)
  len <- length(StsChos)
  RivLabel <- matrix("NA", ncol = len, nrow = len)
  for (i in 1:(len-1))
    for (j in (i+1):len)
      RivLabel[i,j] <- paste(as.character(StsIdx[i]),as.character(StsIdx[j]),sep=",")
  
  upper.idx <- upper.tri(ECemp)
  con.idx <- as.logical(upper.idx * is.con)
  not.con.idx <- as.logical(upper.idx * !(is.con))
  par(mfrow=c(1,sum(which.plots)))
  
  ## Plot ECemp
  if(which.plots[1]){
    plot(DistEuc[not.con.idx],ECemp[not.con.idx], xlab = "Euclidean Distance (km)", ylab="Extremal Coefficient", ylim = c(1,2),xlim = c(0,max(DistEuc)))   
    points(DistEuc[con.idx],ECemp[con.idx],col="blue",pch=4,cex = 1.3)
    if(which.labels[1]) thigmophobe.labels(DistEuc[upper.idx], ECemp[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)
    #legend("topleft", col=c("blue","black"), pch=c(4,19), legend=c("flow-connected pairs","flow-unconnected pairs"), bty = "n", cex =2, pt.cex = 1.5)
    
  }
  
  
  ## Plot ECtheo
  if(which.plots[2]){
    plot(DistCatch[not.con.idx],ECtheo[not.con.idx], xlab = "Hydrological Distance (km)", ylab="Extremal Coefficient",ylim = c(1,2),,xlim = c(0,max(DistCatch))) 
    points(DistCatch[con.idx],ECtheo[con.idx],col="blue",pch=4,cex = 1.3)
    if(which.labels[2]) thigmophobe.labels(DistCatch[upper.idx], ECtheo[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)
    #legend("topleft", col=c("blue","black"), pch=c(4,19),  legend=c("flow-connected pairs","flow-unconnected pairs"), bty = "n", cex = 1.5, pt.cex = 1.5)
  }
  
  
  ## Plot QQ
  if(which.plots[3]){
    plot(ECemp[upper.idx], ECtheo[upper.idx], xlim = c(1,2), ylim = c(1,2),main="Model/Empirical",xlab="Extremal Coefficient - Emperical", ylab="Extremal Coefficient - Fitted Model") 
    abline(b = 1, a = 0)
    points(ECemp[con.idx],ECtheo[con.idx],col="blue",pch=4,cex = 1.3)
    if(which.labels[3]) thigmophobe.labels(ECemp[upper.idx], ECtheo[upper.idx], labels = RivLabel[upper.idx], cex=0.6, offset=0.5)
    legend("topleft", col=c("blue","black"), pch=c(4,19), legend=c("flow-connected pairs","flow-unconnected pairs"), bty = "n", cex =2, pt.cex = 1.5)
  }
  if(PDF) dev.off()
}
################################################################################

exp.measure <- function(w0,w1,lambda)
  1/w0 * pnorm(lambda + log(w1/w0) / (2*lambda)) +
  1/w1 * pnorm(lambda + log(w0/w1) / (2*lambda))

exp.measureYY <- function(w0,w1,lambda)
  1/w0 + 1/w1 - exp.measure(w0,w1,lambda)

exp.measureYN <- function(w0,w1,lambda)
  1/w0 - exp.measureYY(w0,w1,lambda)

exp.measureNY <- function(w0,w1,lambda)
  1/w1 - exp.measureYY(w0,w1,lambda)

spec.dens.HR <- function(w0,w1,lambda)
  1/(2*lambda*w1*w0^2) * dnorm(lambda + log(w1/w0) / (2*lambda))

censored.dens.HR <- function(w0,w1,lambda)
  1/(w0^2) * pnorm(lambda + log(w1/w0) / (2*lambda))


SpecEstimationHR <- function(Data, thresholds = seq(.8,.999, len = 30), do.plot =TRUE)
{
  X0 <- Margin2Pareto(Data[,1],empirical = TRUE)
  X1 <- Margin2Pareto(Data[,2], empirical = TRUE)
  Z <- X0 + X1
  lambda.vec <- numeric(length(thresholds))
  for (i in 1:length(thresholds))
  { idx <- which(Z > quantile(Z, thresholds[i]))
  logllh <- function(lambda) -sum(log(spec.dens.HR(X0[idx], X1[idx], lambda)))
  opt <- optimize(logllh, interval =c(0,50))
  lambda.vec[i] <- opt$minimum
  }
  EC <- Vario2EC(4*lambda.vec^2)
  if(do.plot)
    plot(thresholds, EC)
  return(EC)
}

CensoredEstimationHR <- function(Data, thresholds = seq(.8,.999, len = 30), normalize = TRUE, do.plot = FALSE)
{
  if(normalize){
    X0 <- Margin2Pareto(Data[,1], empirical = TRUE)
    X1 <- Margin2Pareto(Data[,2], empirical = TRUE)
  }
  else{
    X0 <- Data[,1]
    X1 <- Data[,2]
  }     
  
  n <- length(X0)
  lambda.vec <- numeric(length(thresholds))
  for (i in 1:length(thresholds)){
    thr <- thresholds[i]
    u <- quantile(X0, thr)
    logllh <- function(lambda){
      -CensoredLLH.HR(X0, X1, u, rep(lambda, times=n))}       
    opt <- optimize(logllh, interval = c(0, 4))
    lambda.vec[i] <- opt$minimum
  }
  EC <- Vario2EC(4 * lambda.vec ^ 2)
  if (do.plot)
    plot(thresholds, EC)
  return(EC)
}

CensoredLLH.HR <- function(X0,X1,u,lambda)
{
  idxYN <- which(X0 > u & X1 < u)
  idxNY <- which(X1 > u & X0 < u)
  idxYY <- which(X0 > u & X1 > u)
  
  y <- rep(1, times = length(X0))
  y[idxYY] <- spec.dens.HR(X0[idxYY],X1[idxYY],lambda[idxYY]) / exp.measure(u,u,lambda[idxYY])
  y[idxYN] <- censored.dens.HR(X0[idxYN],u,lambda[idxYN]) / exp.measure(u,u,lambda[idxYN])
  y[idxNY] <- censored.dens.HR(X1[idxNY],u,lambda[idxNY])/exp.measure(u,u,lambda[idxNY])
  
  if(sum(log(y)) == -Inf) print("LLH is INFINITE")
  return(sum(log(y)))
}


# plot a scatter plot
ScatterPlot <- function(ScSts, StsInfo, X, filename) {
  pdf(filename)
  par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      mar = c(5,5,4,2) +.1,pty="s")
  xlab <- paste(StsInfo[ScSts[1], 1], " (",as.character(ScSts[1]),")", sep="")
  ylab <- paste(StsInfo[ScSts[2], 1], " (",as.character(ScSts[2]),")", sep="")
  plot(X[,ScSts[1]], X[,ScSts[2]], xlab=xlab, ylab=ylab)
  dev.off()
}

# plot a box plot
BoxPlots <- function(Par.spec,Par.cens,filename) {
  theta1 <- 2.2
  alpha1 <- 1.7
  s1 <- 289.73
  riv.lmb1 <- 0.3
  euc.lmb1 <- .05
  beta1 <- 1.4
  c1 <- .5
  init1 <- c(theta1, alpha1, s1, riv.lmb1, euc.lmb1, beta1, c1)
  nb.rep <- nrow(Par.spec)
  init1Mat <- matrix(rep(init1,nb.rep), nrow = nb.rep, byrow =TRUE)
  Par.specNorm <- Par.spec
  Par.specNorm[,c(2,3,6)] <- log(Par.spec[,c(2,3,6)] / init1Mat[,c(2,3,6)])
  Par.specNorm[,-c(2,3,6)] <- log(Par.spec[,-c(2,3,6)] / init1Mat[,-c(2,3,6)])
  pdf(filename, width=2*7)
  par(mfrow=c(1,2))
  LWD <- 2
  boxplot(Par.specNorm[,c(4,5,1,2,6,7)],
        ylim = c(-1,1),
        boxlwd = LWD,
        staplelwd = LWD,
        whisklwd = LWD,
        main = "Spectral Estimation",
        names = c(expression(lambda[Riv]),
                  expression(lambda[Euc]),
                  expression(tau),
                  expression(alpha),
                  expression(beta),
                  expression(c)),
        pars=list(cex.axis = 1.5, las=2, cex.main = 1.5),
        ylab="",
        outline=FALSE
        )
  abline(h = 0, lwd = 1.6, lty = 2)
  Par.censNorm <- Par.cens
  Par.censNorm[,c(2,3,6)] <- log(Par.cens[,c(2,3,6)] / init1Mat[,c(2,3,6)])
  Par.censNorm[,-c(2,3,6)] <- log(Par.cens[,-c(2,3,6)] / init1Mat[,-c(2,3,6)])
  boxplot(Par.censNorm[,c(4,5,1,2,6,7)],
        ylim = c(-1,1),
        boxlwd = LWD,
        staplelwd = LWD,
        whisklwd = LWD,
        main = "Censored Estimation",
        names = c(expression(lambda[Riv]),
                  expression(lambda[Euc]),
                  expression(tau),
                  expression(alpha),
                  expression(beta),
                  expression(c)),
        pars=list(cex.axis = 1.5, las=2, cex.main = 1.5, ylim=c(0.55, 1.6)),
        ylab="",
        outline=FALSE
        )
  abline(h = 0, lwd = 2, lty = 2)
  dev.off()
}

#########################################################################
####################### Plotting declustered Events #####################

PlottingDeclusterdEvents <- function(Year, Stations, AllRes, ComTSs, StsInfo, filename=filename) {
  Q85 <- numeric()
  for (i in 1:31)
  {
    Q85[i] <- quantile(ComTSs[,(i+1)], 0.85)
  }
  StationsNameTemp <- StsInfo[Stations, 1]
  StationsName <- character()
  for (i in 1:length(Stations)){
    StationsName <- c(StationsName,paste(StsInfo[Stations[i],1],"(",Stations[i],")",sep=""))
  }
  Evs <- AllRes[[which(Years==Year)]]
  Dates <- ComTSs[,1]
  Ind <- which(as.numeric(substr(Dates,1,4))==Year)
  
  TSStations <- ComTSs[Ind,Stations+1]
 # par(mfrow=c(length(Stations),1),oma = c(1,1,1,1) + 0.1,mar = c(0,1,0,1))
  
  Locs <- numeric()
  for (j in 1:length(Stations))
  {
    Locs <- rbind(Locs,Evs$Location[Stations[j],])
  }
  
  LocsCenter <- Evs$CenterEventLoc
  AA <- matrix(0,length(LocsCenter),2)
  for (i in 1:length(LocsCenter))
  {
    AA[i,] <- c(max(1,LocsCenter[i]-4),min(92,LocsCenter[i]+4))
  }
  
  B <- sort(AA[,1],index.return=T) 
  Verticals <- t(AA[B$ix,])
  pdf(filename)
  par(mfrow=c(length(Stations),1),cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19, mar = c(0,1,0,1),oma = c(2,2,2,2))
  
  InGrQ85 <- numeric()
  for (j in 1:length(Stations))
  {
    Ind1 <- (sort(Evs$Location[Stations[j],],index.return=T))$ix
    X1 <- Evs$Location[Stations[j],][Ind1]
    Y1 <- Evs$AllEvMat[Stations[j],][Ind1]
    InGrQ85 <- c(InGrQ85 ,which(Y1 >=Q85[Stations[j]]))
  }
  InGrQ85 <- unique(InGrQ85)
  
  for (j in 1:length(Stations)) {
    X <- TSStations[,j]
    xlim <- c(1,92)
    ylim <- c(min(X),max(X))
    plot(X,type="l",xlim=xlim,ylim=ylim,xlab="",ylab="",axes=F,col="blue",lwd=2)
    box()
    abline(h=Q85[Stations[j]],lty=3,lwd=3)
    
    for (k in InGrQ85)
    {
      abline(v=Verticals[1,k],lty=3)
      abline(v=Verticals[2,k],lty=3)
      polygon(x = c(Verticals[1,k],Verticals[2,k],Verticals[2,k],Verticals[1,k]),
              y = c(0,0,max(X),max(X)),density=10, angle=45, border=NA,col = "gray")
    }
    mtext(StationsName[j], side=4,line = .5,cex = 1.5)
    Label <- round(seq(min(X)+.25*min(X),max(X)-.25*max(X),length.out=2),0)
    
    axis(2,at=Label, labels= Label)
    Ind1 <- (sort(Evs$Location[Stations[j],],index.return=T))$ix
    X1 <- Evs$Location[Stations[j],][Ind1]
    Y1 <- Evs$AllEvMat[Stations[j],][Ind1]
    par(new=T)
    plot(X1[InGrQ85],Y1[InGrQ85],xlim=xlim,ylim=ylim,col="black",xlab="",ylab="",axes=F,pch=19)
  }
  
  Label <- c(paste("Jun 1,",Year,sep=""),paste("July 15,",Year,sep=""),paste("August 31,",Year,sep=""))
  At <- c(4,48,88)
  mtext(Label ,at=At, side=1,line = .5,cex = 1.5)
  mtext(expression(Discharge (m^3/s)), at=min(X)+.3*min(X)+(max(X)-min(X))*length(Stations)/2,side=2,line = 3,cex = 1.5)
  dev.off()
}


#########################################################################################################
#########################################################################################################
################################### Univariate Functions ################################################


############################# Fitting a GEV distribution by Joint PPP likelihood#########################
PPFit <- function(Data, u, NoYears, Init, CovarMu, CovarSc, CovarXi, LogModel, method=method, control=control) {
  NoSt <- length(Data)
  NoParMu <- ncol(CovarMu)
  NoParSc <- ncol(CovarSc)
  NoParXi <- ncol(CovarXi)
  NoPar <- NoParMu + NoParSc + NoParXi
  
  PP.lik.linear <- function(Params) {
    Out <- 0
    for (i in 1:NoSt) {
      mu <- sum(Params[1:NoParMu]*(CovarMu[i, ]))
      sc <- sum(Params[(NoParMu+1):(NoParMu+NoParSc)]*CovarSc[i, ])
      xi <- sum(Params[(NoParMu+NoParSc+1):NoPar]*CovarXi[i, ])
      if (LogModel == TRUE) {
        mu <- exp(mu)
        sc <- exp(sc)
        if (sc==Inf)
          sc <- 0
      }
      y1 <- 1 + xi * ((u[i] - mu)/sc)
      y2 <- 1 + xi * ((Data[[i]] - mu)/sc)
      InterceptSc <- Params[NoParMu + 1]
      InterceptMu <- Params[1]
      
      if ((sc <= 0) | (min(y1) <= 0) | (min(y2) <= 0) | InterceptSc <= 0 | InterceptMu <0) { 
        l <- 10^6
      } else {
        l <- NoYears[i]*(y1 ^ (-1/xi)) + length(Data[[i]]) * log(sc) + (1 / xi + 1) * sum(log(y2))
      }
      Out <- Out+l
    }
    return(Out)
  }
  
  x <- optim(Init, PP.lik.linear, hessian = TRUE, control = control, method = method )
  output <- x
  return(output)
}


ECDF <- function (x) {
  x <- sort(x)
  n <- length(x)+1
  if (n < 1) 
    stop("'x' must have 1 or more non-missing values")
  vals <- unique(x)
  rval <- approxfun(vals, cumsum(tabulate(match(x, vals)))/n, 
                    method = "constant", yleft = 0, yright = 1, f = 0, ties = "ordered")
  class(rval) <- c("ecdf", "stepfun", class(rval))
  assign("nobs", n, envir = environment(rval))
  attr(rval, "call") <- sys.call()
  rval
}


QQPlot <- function (a, dat,main,Bounds,filename) {
  DatSort <- sort(dat)
  EmpCdfFun <- ECDF(DatSort)
  EmpCdf <- EmpCdfFun(DatSort)
  TheQuant <- gevq(a,1-EmpCdf)
  
  xlim1 <- c(min(c(DatSort,TheQuant)),max(c(DatSort,TheQuant)))
  xlim <- c(floor(xlim1[1]/100)*100,ceiling(xlim1[2]/100)*100)
  ylim <- xlim
  
  pdf(filename)
  par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      mar = c(5,5,4,2) +.1,pty="s")
 # par(cex = 1.3, cex.lab = 1.5, cex.axis = 1.5, cex.main = 1.5, pch = 19,
      #mar = c(5,5,4,2) +.1)
  Labels <- round(seq(xlim[1],xlim[2],length.out=4),0)
  plot(TheQuant,DatSort,xlim=xlim,ylim=ylim,xlab="Model Quantile",ylab="Data Quantile",main=main,axes=F)
  box()
  axis(1, at=Labels ,labels=Labels )
  axis(2, at=Labels ,labels=Labels )
  
  abline(0, 1, col = 4)
  par(new=T)
  plot(TheQuant ,Bounds[,1],xlim=xlim,ylim=ylim,main="",xlab="",ylab="",type="l",lty=2,col="blue",axes=F)
  par(new=T)
  plot(TheQuant ,Bounds[,2],xlim=xlim,ylim=ylim,main="",xlab="",ylab="",type="l",lty=2,col="blue",axes=F)
  dev.off()
}


ConfinSimul <- function(Par,R,Data) {
  Data1 <- Data  
  NoData <- length(Data1)
  Temp <- matrix(0,NoData,R)
  for (r in 1:R)
  {
    X <- rgev(NoData,loc=Par[1],scale=Par[2],shape=Par[3])
    Temp[,r] <- sort(X)
  }
  Bounds <- matrix(0,NoData,2)
  
  for (j in 1:NoData)
  {
    Bounds[j,1] <- quantile(Temp[j,],0.025)
    Bounds[j,2] <- quantile(Temp[j,],0.975)
  }
  return(Bounds)
}
