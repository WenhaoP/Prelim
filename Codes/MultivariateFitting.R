source("Codes/Functions.R")

### Loading Data ###
FileList <- list.files(path = "Data")
for (i in 1:length(FileList)){
  load(paste("Data", FileList[i], sep="/"))
}

StsChos <- c(1:31) # Station indices
NoSt <- length(StsChos) # Total number of stations

# CatchEucDis <- as.matrix(dist(CatCntr))
CatchEucDisWt <- as.matrix(dist((CatCntrWt))) # Hydrological distances

### Declustering and obtaining events ###
Years <- unique(as.numeric(substr(ComTSs[, 1], 1, 4)))
U <- 0.1
Lag <- 4
Plotting <- 0
StNames <- as.character(StsChos)
YearsWithEvent <- numeric()
AllEvMat <- matrix(0, length(StNames), 1)
AllRes <- list()
for (i in 1:length(Years)){
  cat(paste("Declustering year", Years[i], "\n"))
  Res <- ObtainMultVarEvents(
    TSs = ComTSs, U = U, Lag = Lag, Year = Years[i],
    StNames = StNames, Plotting = Plotting, mfrow = c(NoSt,1))  ## For Censored Likelihood
  Events <- Res$AllEvMat
  Locations <- Res$Location
  if (length(Events) > 0) {
    YearsWithEvent <- c(YearsWithEvent, rep(Years[i], dim(Events)[2]))
    AllEvMat <- cbind(AllEvMat, Events)
  }
  AllRes[[i]] <- Res 
}

DataEvents <- t(AllEvMat[, -1])
rownames(DataEvents) <- YearsWithEvent 

### Compute EC for HR model with censored estimation ###
Theta.cen <- matrix(1, ncol = NoSt, nrow = NoSt)

for(i in 1:(NoSt-1))
    for(j in (i+1):NoSt) {
    Data.biv <- DataEvents[,c(i,j)]
    Theta.cen[i,j] <- CensoredEstimationHR(Data.biv, thresholds = .90)
    Theta.cen[j,i] <- Theta.cen[i,j]
}

### Transform to Pareto ###
TSNew <- matrix(0,dim(DataEvents)[1],dim(DataEvents)[2])
for (i in 1:dim(DataEvents)[2]) {
  TSNew[,i] <- Margin2Pareto(DataEvents[,i], empirical = TRUE) 
}

### Fitting of the max-stable process on river network ###

theta0 <- 2 #1*max(RiverDisChos)
alpha0 <- 1.6
s0 <- 289.73 ##max(CatCntrWtKmChos) ## FIX
riv.lmb0 <- .3
euc.lmb0 <- .05
beta0 <- 1.1
c0 <- .5

### SPECTRAL ESTIMATION ###

fit2 <- pareto_BR_River(
  data = TSNew, 
  riv.dist = RiverDis,
  euc.coord = CatCntrWt, 
  is.connected = FlowCon,
  riv.weights = Weight,
  u=quantile(rowSums(TSNew), 0.87),
  init= c(theta0, alpha0, s0, riv.lmb0, euc.lmb0, beta0, c0),
  fixed=c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
  maxit = 1000,
  method= "BFGS") 

### CENSORED ESTIMATION: This is computationally expensive (to skip it, just go to "par  <- fitM4$par" for the result which is already loaded)###

fitM4 <- pareto_BR_River_cen(
    data = TSNew,
    riv.dist = RiverDis,
    euc.coord = CatCntrWt,
    riv.weights = Weight,
    u=quantile(TSNew, .90),   
    init= c(theta0, alpha0, s0, riv.lmb0, euc.lmb0, beta0, c0),
    fixed=c(FALSE, FALSE, TRUE, FALSE, FALSE, FALSE, FALSE),
    maxit = 200,
    method= "BFGS")

par  <- fitM4$par 

### Transform estimates to the scale used in the paper ###
par2paper <- function(par)
  return(c(2 * par[4], par[5] / par[3], par[1] * 500, par[2], par[6], par[7]))
Par.paper <- par2paper(par)

### Estimates obtained from bootstraping the data ###
Par.boot.paper <- apply(Par.boot.cens[1:100,], MARGIN = 1, FUN = par2paper)
SDs <- apply(Par.boot.paper, MARGIN = 1, FUN = sd)


### Codes to obtain the plots in the paper ###

### Fig. 2 ###

# Parameters of the transformation matrix
beta <- par[6]
c <- par[7]

rot.mat <- cbind(c(cos(beta), c*sin(beta)),c(-sin(beta), c*cos(beta))) # Transformation matrix
Catch.dist <- as.matrix(dist(CatCntrWt %*% t(rot.mat)))  # Hydrological distance with the transformation


plotEucVsHyd(ECemp = ThetaOfBlockMaxima,
             ECtheo = ThetaOfBlockMaxima,
             DistEuc = EucDis,
             DistCatch = Catch.dist,
             is.con = FlowCon,
             StsIdx = StsChos,
             which.plots = c(TRUE,TRUE,FALSE),
             which.labels = c(FALSE,FALSE,TRUE),
             PDF = TRUE,
             filename = paste("Plots","EucVSCatch.pdf",sep="/"))


### Fig. 5 ###

PlottingDeclusterdEvents(Year = 1989,
                         Stations = c(1,2,11,13),
                         AllRes = AllRes,
                         ComTSs = ComTSs,
                         StsInfo = StsInfo,
                         filename = paste("Plots", "MultiEv1.pdf", sep="/"))

### Fig. 7 ###

### Left ###
ScSts <- c(2,5)
                                        
filename <- paste(paste(as.character(ScSts[1]),as.character(ScSts[2]),sep="-"),".pdf",sep="")

ScatterPlot(ScSts=ScSts,
            StsInfo,
            X=TSNew,
            filename = paste("Plots",filename,sep="/"))

### Left ###
ScSts <- c(13,16)
                                        
filename <- paste(paste(as.character(ScSts[1]),as.character(ScSts[2]),sep="-"),".pdf",sep="")

ScatterPlot(ScSts=ScSts,
            StsInfo,
            X=TSNew,
            filename = paste("Plots",filename,sep="/"))

### Fig. 8 ###

euc.dist <- as.matrix(dist(CatCntrWt %*% t(rot.mat)))    

### Left ###

plotECF(ECemp = ThetaOfBlockMaxima,
        ECtheo = Theta.cen,
        Dist = CatchEucDisWt,
        is.con = FlowCon,
        StsIdx  = StsChos,
        which.plots = c(FALSE,FALSE,TRUE),
        which.labels = c(FALSE, FALSE,FALSE),
        PDF = TRUE,
        filename = paste("Plots","ECF_Emp_HR.pdf",sep="/"))

### Center ###

plotECF(ECemp = ThetaOfBlockMaxima,
        ECtheo = Vario2EC(
            vario(par,
                  riv.dist = RiverDis,
                  euc.coords = CatCntrWt,
                  riv.weights = Weight
                  )
            ),
        Dist = euc.dist,
        is.con = FlowCon,
        StsIdx  = StsChos,
        which.plots = c(FALSE,TRUE,FALSE),
        which.labels = c(FALSE, FALSE,FALSE),
        PDF = TRUE,
        filename = paste("Plots","EmpVSModel2.pdf",sep="/")) 


### Right ###

plotECF(ECemp = ThetaOfBlockMaxima,
        ECtheo = Vario2EC(
            vario(par,
                  riv.dist = RiverDis,
                  euc.coords = CatCntrWt,
                  riv.weights = Weight
                  )
            ),
        Dist = euc.dist,
        is.con = FlowCon,
        StsIdx  = StsChos,
        which.plots = c(FALSE,FALSE,TRUE),
        which.labels = c(FALSE, FALSE,FALSE),
        PDF = TRUE,
        filename = paste("Plots","EmpVSModel3.pdf",sep="/")) 

### Fig. 9 ###

require("SpatialExtremes")
require("mvtnorm")
library("ismev")
Maxima <- CommonMaxima

### N=NoSt=31 reproduce the bottom-right plot of Fig. 9. For other random groups, change N ###

N <- 31
StsForPlot <- sort(sample(StsChos,N))

### Transformation to Unit Frechet ###

TransMax <- matrix(0,nrow(Maxima),N) 
for (i in 1:N){
    Pars <- gev.fit(Maxima[,StsForPlot[i]])$mle
    TransMax[,i] <- gev2frech(Maxima[,StsForPlot[i]], Pars[1], Pars[2], Pars[3], emp = FALSE)
}

### Here the data quantiles are computed ###
block.size <- 1
nb.blocks <- nrow(Maxima)
TSmax <- double(nb.blocks)
TSmax <- apply(TransMax,1,max)
quants.TSmax <- sort(TSmax)

### Here the model quantiles are computed ###

d <- ncol(TransMax)
V <- function(x, par){
    d <- length(x)    
    varioHalf <- 1/2 * vario(par,
                             riv.dist = RiverDis[StsForPlot,StsForPlot],
                             euc.coords = CatCntrWt[StsForPlot,],
                             riv.weights = Weight[StsForPlot,StsForPlot]
                             )            
    
    f1 <- function(i,x){                
        vario0 <- matrix(rep(varioHalf[-i,i],d-1),ncol = d-1)
        S <- vario0 + t(vario0) - varioHalf[-i,-i]
        return(1/x[i]*pmvnorm(upper=(log(x/x[i])+varioHalf[,i])[-i],mean=rep(0,d-1),sigma= S)[1])
    }
    return(sum(apply(cbind(1:d),1,f1,x=x)))
}


quant.level <- ((1:nb.blocks) - 1/2) / nb.blocks

V11 <- V(x = rep(1, times = d),par= fitM4$par)

quants.modelMax <-  - V11 * 1 / log(quant.level)

quants.indep <- - length(StsForPlot) * 1 / log(quant.level) ## quantiles if all points were independent
quants.dep <- - 1 / log(quant.level) ## quantiles if all points complete dependent

pdf(paste("Plots/QQMaxima_",length(StsForPlot),".pdf", sep=""))
par(cex = 1.3, cex.lab = 2.3, cex.axis = 2.3, cex.main = 1.5, pch = 19,
    mar = c(5,5,4,2) +.1)
plot(log(quants.modelMax),log(quants.TSmax),xlab="Model Quantile",ylab="Data Quantile")
abline(a = 0, b =1,col="blue")
lines(log(quants.modelMax), log(quants.indep),col="blue",lty=2) 
lines(log(quants.modelMax), log(quants.dep),col="blue",lty=3)
dev.off()




### Fig. 1 of supplementary material ###

BoxPlots(Par.spec=Par.sim.spec,
         Par.cens=Par.sim.cens,
         filename = paste("Plots","Sim_Study.pdf",sep="/"))


