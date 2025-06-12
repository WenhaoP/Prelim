library(ismev)
library(evd)

source("Codes/Functions.R")

### Loading Data ###
FileList <- list.files(path="Data", pattern = "\\.RData$")
for (i in 1:length(FileList)){
  load(paste("Data",FileList[i],sep="/"))
}

StsChos <- c(1:31)
NoSt <- length(StsChos)

### GEV analysis  ###
method <- "BFGS" 
control <- list(maxit = 1000, reltol=10^(-30), abstol=0.0001, trace=2)

Maxima <- SummerMaxima

ResGev <- list()
GevPars <- matrix(0, NoSt, 3)
FreeGevNLLTemp <- numeric()
for (i in 1:NoSt){
    ResGev[[i]] <- gev.fit(Maxima[[i]])
    GevPars[i,] <- ResGev[[i]]$mle
    FreeGevNLLTemp[i] <- ResGev[[i]]$nll
}

### Extraction of univariate events  ###
U <- 0.9
Lag <- 5
StNames <- as.character(1:31)
mfrow <- c(1,1)
YearsWithEvents <- list()
Threshold <- numeric()
AllEvents <- list()
NoOfYears <- numeric()
for (i in 1:NoSt) {
	TSs <- TSMonths[[i]]
	X <- TSMonths[[i]]$Val
	Threshold[i] <- quantile(X,U)
	Years <- unique(as.numeric(substr(TSs[,1],1,4)))
	Events <- numeric()
	YearEvents <- numeric()
	for (j in 1:length(Years)) {
		X1 <- ObtainMultVarEvents(TSs = TSMonths[[i]],U=U,Lag=Lag,Year=Years[j],StNames=StNames[i],Plotting=0,mfrow)
		X2 <- as.numeric(X1$AllEvMat)
		Events  <- c(Events,X2)
		YearEvents <- c(YearEvents,rep(Years[j],length(X2)))
	}
	AllEvents[[i]] <- Events
	YearsWithEvents[[i]] <- YearEvents
	NoOfYears[i] <- length(Years)
}

### Fitting a GEVD to each station (by maximizing the joint Poisson process likelihood; cf. formula (21) in the paper) ###
PPRes <- list()
ParsPP <- matrix(0,NoSt,3)
PPLLFree <- numeric()
for (i in 1:NoSt) {
	Init <- GevPars[i,]
	PPRes[[i]] <- PPFit(Data = AllEvents[i], u = Threshold[i], NoYears = NoOfYears[i],
	                    Init = Init, CovarMu = matrix(1,1,1), CovarSc = matrix(1,1,1),
	                    CovarXi = matrix(1,1,1), LogModel = FALSE, method = method, control = control)
	ParsPP[i,] <- PPRes[[i]]$par
	PPLLFree[i] <- PPRes[[i]]$value
}

M1LL <- sum(PPLLFree)

### Regionalized model with covariates; cf. formula (32) ###
Grp1 <- c(11,12,16,17,18,19,21,22)
Grp2 <- c(13,28:31)
Grp3 <- c(1:10,14,15,20)
Grp4 <- c(23,24,25,26,27)
Grp <- list(Grp1,Grp2,Grp3,Grp4)

Lat <- StsCenter[,2]
CovarStsAll <- cbind(rep(1,length(StsArea)),StsArea,StsAlt,StsSlope,StsDensity,Lat)
PPCovariateModels <- list()
CovarLL <- numeric()
ParsCovModel <- matrix(0,31,3)
Normalized <- FALSE
CovarLog <- cbind(CovarStsAll[,1],log(CovarStsAll[,2:6]))

### The covariates which are siginficant, obtained by log-likelihood ratio test###
MuCovars <- list(c(2,3,4,6), c(2,4,6), c(2,3,4,6), c(2,4,6))
ScCovars <- list(c(2,3,4), c(2,4), c(2,3,4), c(2,4))
InitAll <- list(c(1.1, 0.9, 0.7, 1.2, 0, -0.1), 
                c(0.6, 0.5, 0.4, 1.7, 0.1, -0.3), 
                c(1.2, 0.7, 0.9, 0.7, -0.1, -0.5), 
                c(0.1, 1, -0.1, 1.3, 0.5, 0.1))

for (i in 1:length(Grp)) {
	GrSts <- Grp[[i]]
	CovarMu <- CovarLog[,MuCovars[[i]]]
	InitMu <- InitAll[[i]][MuCovars[[i]]]
	CovarSc <- CovarLog[,ScCovars[[i]]]
	InitSc <- InitAll[[i]][ScCovars[[i]]]
	InitXi <- mean(ParsPP[GrSts,3])
	CovarMuReg <- CovarMu[GrSts,]
	CovarScReg <-  CovarSc[GrSts,]
	CovarXi <- as.matrix(CovarLog[GrSts ,1])
	Init <- c(InitMu,InitSc,InitXi)
	PPCovariateModels[[i]] <- PPFit(Data = AllEvents[GrSts], u = Threshold[GrSts],
	                                NoYears = NoOfYears[GrSts], Init = Init, CovarMu = CovarMuReg, CovarSc = CovarScReg,
	                                CovarXi = CovarXi, LogModel = TRUE, method = method, control = control)
	CovarLL[i] <- PPCovariateModels[[i]]$value
	ParsModelTemp <- PPCovariateModels[[i]]$par
	MuCov <- exp((CovarMuReg %*% ParsModelTemp [(1):(ncol(CovarMu))]))
	ScCov <- exp((CovarScReg %*%ParsModelTemp[(ncol(CovarMu)+1):(ncol(CovarMu)+ncol(CovarSc))]))
	XiCov <- as.matrix(CovarXi) %*%  ParsModelTemp[(ncol(CovarMu)+ncol(CovarSc)+1):length(ParsModelTemp)]
	ParsCovModel[Grp[[i]], ] <- cbind(MuCov, ScCov, XiCov)
}

M2LL <- sum(CovarLL)


M1AIC <- 2*(NoSt*3+M1LL)
M2AIC <- 2*(length(Grp)+length(unlist(ScCovars))+length(unlist(MuCovars))+M2LL)

### Fig. 2 (Left) in supplementary material ###
ChosStForPlot<- c(7,13)

k <- 1
R <- 100000 ## put to 100000 to reproduce the figure in the paper
for(i in ChosStForPlot){	
    Bounds <- ConfinSimul(Par=ParsCovModel[i,],R=R,Data=ResGev[[i]]$data)
    Temp1 <- paste(StsInfo[i,1]," (",sep="")
    Temp2 <- paste("Station",as.character(i))
    main <- paste(Temp1,Temp2,")",sep="")
    par(pty="s")
    QQPlot(a=ParsCovModel[i,],dat=ResGev[[i]]$data,main=main,Bounds=Bounds,
           filename=paste("Plots",paste("PPCovar",as.character(i),".pdf",sep=""),sep="/"))
}

### Fig. 2 (Right) in supplementary material ###
### Compute 100y return levels on the whole network; cf. Fig. 6 ###
NoIn <- nrow(IntPolPointCenter)

GrpIn1 <- which(IntPolPointCenter[,1] < 11.6 & IntPolPointCenter[,2] < 47.72)
GrpIn2 <-  which(IntPolPointCenter[,1] >= 11.6 & IntPolPointCenter[,2] < 47.72 )
GrpIn4 <- which(IntPolPointCenter[,2] > 48.65 )
GrpIn3 <- c(1:NoIn)[-c(GrpIn1,GrpIn2,GrpIn4)]
GrpIn <- list(GrpIn1,GrpIn2,GrpIn3,GrpIn4)

CovarIn <- cbind(rep(1,length(IntPolPointArea)),IntPolPointArea,IntPolPointAlt,IntPolPointSlope,IntPolPointDensity,IntPolPointCenter[,2])
CovarInLog <- cbind(CovarIn[,1],log(CovarIn[,2:6]))

ParsIn <- matrix(0,nrow(IntPolPointCenter),3)

for (i in 1:length(GrpIn)) {
  GrReg <- GrpIn[[i]]
	CovarMuIn <- CovarInLog[GrReg,MuCovars[[i]]] 
	CovarScIn <- CovarInLog[GrReg,ScCovars[[i]]] 
	CovarXi <- as.matrix( CovarInLog[GrReg,1])
	ParsModelTemp <- PPCovariateModels[[i]]$par
	MuCovIn <- exp((CovarMuIn %*% ParsModelTemp [1:length(MuCovars[[i]])]))
	ScCovIn <- exp(CovarScIn %*% ParsModelTemp[(length(MuCovars[[i]])+1):(length(MuCovars[[i]])+length(ScCovars[[i]]))])
	XiCovIn <-CovarXi  %*%  ParsModelTemp[length(ParsModelTemp)]
	ParsIn[GrpIn[[i]],] <- cbind(MuCovIn,ScCovIn,XiCovIn)
}

HundRetLev <- numeric()
for (i in 1:NoIn) {
  HundRetLev[i] <- qgev(.99,loc=ParsIn[i,1],scale=ParsIn[i,2],shape=ParsIn[i,3])
}

HundRetLevelSts <- numeric()
for (i in 1:NoSt){
  HundRetLevelSts [i] <- qgev(.99,loc=ParsPP[i,1],scale=ParsPP[i,2],shape=ParsPP[i,3])
}
