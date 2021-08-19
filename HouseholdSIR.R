# Household SIR Model
# Carolyn LaPrete
# 8/17/2021

# Function file ####
source("/Users/JayLaPrete/Desktop/Household SIR Model/HouseholdSIR-funcs.R")

# Define parameters ###########################################################
## total population size (approx desired) ####
tPop <- 1000

## household size ####
n <- 2.42

## number of households ####
m <- tPop/n

## total population vector split into groups ####
# N <- rep(n,m)
N <- round(rpois(m,n-1)+1)  # (Poisson distributed)
N # print check

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
# also used as starting point for alphaL/G (normalize for per capita rates)

### distance between groups i and j are at d[i,j] ####
# d <- matrix(0,m,m)
# d <- fill_d(d,m)
# d # print check

### contact matrix of alphaLs and alphaGs ####
A <- matrix(nrow=m,ncol=m)
A <- fill_A_det(A,m,N,alpha)
# A <- fill_A_stoch_varG(A,m,N,alpha)
# A <- fill_A_stoch_funcd(A,m,N,alpha)
A # print check

# Define times to solve #######################################################
tmax <- 140
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
#SIRout <- det_disc(So,Io,Ro,A,gamma,SIRout)
#baseSIRout <- base_det(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout,max(SIRout$time))
SIRout <- stoch_disc(So,Io,Ro,A,gamma,SIRout)
baseSIRout <- base_stoch(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout,max(SIRout$time))

# for averages
for(i in 1:9){
  SIRout1 <- data.frame(time=1,S=sum(Sinit),I=sum(Iinit),R=sum(Rinit),t(init))
  baseSIRout1 <- data.frame(time=1,S=bSo,I=bIo,R=bRo)
  SIRout1 <- stoch_disc(So,Io,Ro,A,gamma,SIRout1)
  baseSIRout1 <- base_stoch(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout1,max(SIRout1$time))
  
  SIRout <- SIRout + SIRout1
  baseSIRout <- baseSIRout + baseSIRout1
}
SIRout <- SIRout/10
baseSIRout <- baseSIRout/10

# print test (tPop stays same)
rowSums(SIRout[,c(2,3,4)])
rowSums(baseSIRout[,c(2,3,4)])
tPop

# Output data to .csv and parameters to .txt ##################################
write.csv(SIRout,"/Users/JayLaPrete/Desktop/SIRout-stoch3.csv")
write.csv(baseSIRout,"/Users/JayLaPrete/Desktop/baseSIRout-stoch3.csv")
sink("/Users/JayLaPrete/Desktop/param-stoch3.txt")
cat(" tPop ========================================================= \n")
tPop
cat("\n average house size =========================================== \n")
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
cat("\n average alphaL =============================================== \n")
mean(diag(A))
cat("\n average alphaG =============================================== \n")
mean(A[row(A)!=col(A)])
cat("\n N ============================================================ \n")
N
cat("\n base init ==================================================== \n")
base_init
cat("\n init ========================================================= \n")
init
#cat("\n d ============================================================ \n")
#d
cat("\n A ============================================================ \n")
A
cat("\n")
sink()

# Graph #######################################################################
pdf(file = "/Users/JayLaPrete/Desktop/plot3.pdf")

# if(any(Iinit == 0)){ # 3 graphs if any group starts without infected
#   attach(mtcars)
#   layout(matrix(c(1,1,2,3), 2, 2, byrow = TRUE))
# } else { # 2 graphs if all groups start with at least one infected
#   par(mfrow=c(2,1))
# }

# Graph base
plot(baseSIRout$S ~ time,baseSIRout,type="l",lty=1,lwd=2, ylim=c(0,tPop),
     ylab="Individuals",xlab="Time (days)",
     cex.axis=1,cex.lab=1, col="greenyellow")
lines(baseSIRout$I ~ time,baseSIRout,lty=1,lwd=2,col="pink",
      cex.axis=1,cex.lab=1)
lines(baseSIRout$R ~ time,baseSIRout,lty=1,lwd=2,col="light blue",
      cex.axis=1,cex.lab=1)
# Total SIR
lines(SIRout$S ~ time,SIRout,lty=1,lwd=2,col="dark green",
      cex.axis=1,cex.lab=1)
lines(SIRout$I ~ time,SIRout,lty=1,lwd=2,col="dark red",
      cex.axis=1,cex.lab=1)
lines(SIRout$R ~ time,SIRout,lty=1,lwd=2,col="dark blue",
      cex.axis=1,cex.lab=1)
legend("right",c("S","I","R","Standard S","Standard I","Standard R"),cex=0.7,
       col=c("dark green","dark red","dark blue","greenyellow","pink","light blue"),
       lty=c(1,1,1,1,1,1),lwd=c(2,2,2,2,2,2),bg="white")
title(main=paste0("SIR vs. Standard SIR \n P = ",tPop,", n = ",n))

# # SIR 1 group with infected start
# plot(SIRout$S1 ~ time,SIRout,type="l",lty=1,lwd=2, ylim=c(0,max(N)),
#      ylab="Individuals",xlab="Time (days)",
#      cex.axis=1,cex.lab=1, col="green3")
# lines(SIRout$I1 ~ time,SIRout,lty=1,lwd=2,col="red",
#       cex.axis=1,cex.lab=1)
# lines(SIRout$R1 ~ time,SIRout,lty=1,lwd=2,col="blue",
#       cex.axis=1,cex.lab=1)
# legend("right",c("S","I","R"),cex=0.7,
#        col=c("green3","red","blue"),
#        lty=c(1,1,1),lwd=c(2,2,2),bg="white")
# title(main="SIR 1 group \n with infected start")
# 
# # SIR 1 group without infected start (only used if exist)
# if(any(Iinit == 0)){
#   plot(SIRout$S20 ~ time,SIRout,type="l",lty=1,lwd=2, ylim=c(0,max(N)),
#        ylab="Individuals",xlab="Time (days)",
#        cex.axis=1,cex.lab=1, col="green3")
#   lines(SIRout$I20 ~ time,SIRout,lty=1,lwd=2,col="red",
#         cex.axis=1,cex.lab=1)
#   lines(SIRout$R20 ~ time,SIRout,lty=1,lwd=2,col="blue",
#         cex.axis=1,cex.lab=1)
#   legend("right",c("S","I","R"),cex=0.7,
#          col=c("green3","red","blue"),
#          lty=c(1,1,1),lwd=c(2,2,2),bg="white")
#   title(main="SIR 1 group \n without infected start")
# }

dev.off()

