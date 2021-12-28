# Household SIR Model (functions) - Fall 2021
# Carolyn LaPrete
# 12/28/2021

# fill A with matrix of alphaL at diagonal and alphaG everywhere else
fill_A_det <- function(A,m,N,alpha,aL,aG){
  A <- matrix(aG,m,m)
  for(i in 1:m){
    A[i,i] <- aL
  }
  return(A)
}

# stochastic discrete
stoch_disc <- function(So,Io,Ro,A,gamma,SIRout){
  while(time<tmax){ # for averages
    #while(sum(Io)>0 & time<tmax){
    deltat <- 1
    time <- time + deltat
    
    p_SI <- abs(1 - exp(-(A%*%Io))) # probability of transition S to I
    p_IR <- 1 - exp(-gamma) # probability of transition I to R
    deltaSI <- rbinom(So,So,c(p_SI))
    deltaIR <- rbinom(Io,Io,p_IR)
    
    Sn <- So - deltaSI
    In <- Io + deltaSI - deltaIR
    Rn <- Ro + deltaIR
    SIRout <- rbind(SIRout,c(time,sum(Sn),sum(In),sum(Rn),t(c(Sn,In,Rn))))
    
    So <- Sn
    Io <- In
    Ro <- Rn
  }
  return(SIRout)
}

# base (no groups) for comparison - stochastic
base_stoch <- function(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout,maintime){
  while(time<tmax){ # for averages
    #while(sum(Io)>0 & time<tmax & time < maintime){
    deltat <- 1
    time <- time + deltat
    
    bp_SI <- abs(1 - exp(-((alpha/tPop)*bIo))) # probability of transition S to I
    bp_IR <- 1 - exp(-gamma) # probability of transition I to R
    deltabSI <- rbinom(1,bSo,c(bp_SI))
    deltabIR <- rbinom(1,bIo,bp_IR)
    
    bSn <- bSo - deltabSI
    bIn <- bIo + deltabSI - deltabIR
    bRn <- bRo + deltabIR
    baseSIRout <- rbind(baseSIRout,c(time,bSn,bIn,bRn))
    
    bSo <- bSn
    bIo <- bIn
    bRo <- bRn
  }
  return(baseSIRout)
}
