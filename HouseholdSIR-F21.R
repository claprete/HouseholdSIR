# Household SIR Model - Fall 2021
# Carolyn LaPrete
# 12/28/2021

# Function file ####
source("/Users/JayLaPrete/Desktop/HouseholdSIR-funcs-F21.R")

# Define parameters ###########################################################
## total population size (approx desired) ####
tPop <- 1000

## household size ####
n <- 4

## number of households ####
m <- tPop/n

## total population vector split into groups ####
N <- rep(n,m)

## total population size (calculated) ####
tPop <- sum(N)

## total individuals initially infected (to be spread 1 per group) ####
infected_init <- round(0.01*tPop)

## Rates ####
### recovery rate ####
gamma <- 1/14

R0 <- 2.5

### transmission rate ####
alpha <- R0*gamma

### probability of infection (given) local/global ####
k <- 1

### local alpha ####
aL <- alpha/((n-1)+((1/k)*(tPop-n)))

### global alpha ####
aG <- (1/k)*aL

### contact matrix of alphaLs and alphaGs ####
A <- matrix(nrow=m,ncol=m)
A <- fill_A_det(A,m,N,alpha,aL,aG)

# Define times to solve #######################################################
tmax <- 270
times <- seq(from=0,to=tmax,by=1)

# Define initial condition ####################################################
Iinit <- rep(0,m)       # m households with no infectives
while(sum(Iinit) < infected_init){
  for(i in 1:m){
    if((Iinit[i] < N[i]) & (sum(Iinit) < infected_init))
      Iinit[i] = Iinit[i] + 1
  }
}

Sinit <- N - Iinit      # m households of n individuals all susceptible
Rinit <- rep(0,m)       # m households with no removed

init <- c(Sinit,Iinit,Rinit)

# Assign names to initial conditions
Snms <- paste0("S",1:m)
Inms <- paste0("I",1:m)
Rnms <- paste0("R",1:m)
names(init) <- c(Snms,Inms,Rnms)

base_init <- c(tPop-infected_init,infected_init,0)
names(base_init) <- c("baseS","baseI","baseR")

init # print check
base_init

# Solve the equation numerically ##############################################
## initialize ####
time <- 1
So <- Sinit
Io <- Iinit
Ro <- Rinit
names(So) <- c(Snms)
names(Io) <- c(Inms)
names(Ro) <- c(Rnms)
SIRout <- data.frame(time=1,S=sum(Sinit),I=sum(Iinit),R=sum(Rinit),t(init))

bSo <- tPop-infected_init
bIo <- infected_init
bRo <- 0
baseSIRout <- data.frame(time=1,S=bSo,I=bIo,R=bRo)

## run ####
SIRout <- stoch_disc(So,Io,Ro,A,gamma,SIRout)
baseSIRout <- base_stoch(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout,max(SIRout$time))
esc_hh <- 0

# for averages
count <- 1
for(i in 1:49){
  SIRout1 <- data.frame(time=1,S=sum(Sinit),I=sum(Iinit),R=sum(Rinit),t(init))
  baseSIRout1 <- data.frame(time=1,S=bSo,I=bIo,R=bRo)
  SIRout1 <- stoch_disc(So,Io,Ro,A,gamma,SIRout1)
  baseSIRout1 <- base_stoch(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout1,max(SIRout1$time))

  count <- count + 1
  SIRout <- SIRout + SIRout1
  baseSIRout <- baseSIRout + baseSIRout1
  
  esc_hh1 <- 0
  Ss <- SIRout1[tmax,c(5:(4+m))]
  for(s in Ss){
    if(s==n){
      esc_hh1 = esc_hh1 + 1
    }
  }
  esc_hh1
  esc_hh <- esc_hh + esc_hh1
}
SIRout <- SIRout/count
baseSIRout <- baseSIRout/count
esc_hh <- esc_hh/count

# Output data to .csv and parameters to .txt ##################################
write.csv(SIRout,"/Users/JayLaPrete/Desktop/SIRout-stoch.csv")
write.csv(baseSIRout,"/Users/JayLaPrete/Desktop/baseSIRout-stoch.csv")
sink("/Users/JayLaPrete/Desktop/param-stoch.txt")
cat(" tPop ========================================================= \n")
tPop
cat("\n house size =========================================== \n")
n
cat("\n infected_init ================================================ \n")
infected_init
cat("\n SIRout ======================================================= \n")
tail(SIRout[,c(1,2,3,4)],1)
cat("\n baseSIRout =================================================== \n")
tail(baseSIRout[,c(1,2,3,4)],1)
cat("\n max I ======================================================== \n")
max(SIRout$I)
cat(" at time ")
which.max(SIRout$I)
cat("\n max base I =================================================== \n")
max(baseSIRout$I)
cat(" at time ")
which.max(baseSIRout$I)
cat("\n alphaL =============================================== \n")
mean(diag(A))
cat("\n alphaG =============================================== \n")
mean(A[row(A)!=col(A)])
cat("\n N ============================================================ \n")
N
cat("\n base init ==================================================== \n")
base_init
cat("\n init ========================================================= \n")
init
cat("\n A ============================================================ \n")
A
cat("\n")
sink()

# Graph #######################################################################
fname = paste("/Users/JayLaPrete/Desktop/k",round(k),"n",n,"P",tPop,".pdf",sep="")
pdf(file = fname)

# Graph base
plot(baseSIRout$S ~ time,baseSIRout,type="l",lty=1,lwd=2, ylim=c(0,tPop),
     ylab="Individuals",xlab="Time (days)",
     cex.axis=1,cex.lab=1, col="greenyellow")
lines(baseSIRout$I ~ time,baseSIRout,lty=1,lwd=2,col="pink",
      cex.axis=1,cex.lab=1)
lines(baseSIRout$R ~ time,baseSIRout,lty=1,lwd=2,col="light blue",
      cex.axis=1,cex.lab=1)
# Total SIR
lines(SIRout$S ~ time,SIRout,lty=2,lwd=2,col="dark green",
      cex.axis=1,cex.lab=1)
lines(SIRout$I ~ time,SIRout,lty=2,lwd=2,col="dark red",
      cex.axis=1,cex.lab=1)
lines(SIRout$R ~ time,SIRout,lty=2,lwd=2,col="dark blue",
      cex.axis=1,cex.lab=1)
legend("right",c("S","I","R","Standard S","Standard I","Standard R"),cex=0.7,
       col=c("dark green","dark red","dark blue","greenyellow","pink","light blue"),
       lty=c(2,2,2,1,1,1),lwd=c(2,2,2,2,2,2),bg="white")
title(main=paste0("SIR vs. Standard SIR \n P = ",tPop,", n = ",n))

dev.off()

# print to console
k
aL
aG

ss=tail(baseSIRout$S,1)/head(baseSIRout$S,1)
-log(ss)/(1-ss)
sh=tail(SIRout$S,1)/head(SIRout$S,1)
-log(sh)/(1-sh) # calculated R0