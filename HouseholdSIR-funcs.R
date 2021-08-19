# Household SIR Model (functions)
# Carolyn LaPrete
# 8/5/2021

d_func <- function(d){
  return(abs(exp(-0.02*d)))
}

# fill d
fill_d <- function(d,m){
  for(i in 1:m){
    for(j in 1:m){
      if(i != j){
        if(d[i,j] == 0){
          d[i,j] <- round(rnorm(1,abs(i-j),abs(i-j/3)))
          d[j,i] <- d[i,j]
        }
      }
    }
  }
  return(d)
}

# fill A with matrix of alphaL at diagonal and alphaG everywhere else
fill_A_det <- function(A,m,N,alpha){
  A <- matrix(alpha,m,m)
  A <- A/(sum(N)-N)
  for(i in 1:m){
    A[i,i] <- alpha/N[i]
  }
  return(A)
}

# fill A with matrix of variable alphaL at diagonal and variable alphaG
fill_A_stoch_varG <- function(A,m,N,alpha){
  for(i in 1:m){
    for(j in 1:m){
      if(i == j){
        A[i,j] <- rnorm(1,alpha,alpha/3)/N[i]
      } else{
        if(is.na(A[i,j])){
          A[i,j] <- rnorm(1,alpha,alpha/3)/(sum(N)-N[i])
          if(A[i,j] < 0)
            A[i,j] <- 0
          A[j,i] <- A[i,j]
        }
      }
    }
  }
  return(A)
}

# fill A with matrix of variable alphaL at diagonal and alphaG as func of distance
fill_A_stoch_funcd <- function(A,m,N,alpha){
  for(i in 1:m){
    for(j in 1:m){
      if(i == j){
        A[i,j] <- rnorm(1,alpha,alpha/3)/N[i]
      } else{
        if(is.na(A[i,j])){
          A[i,j] <- rnorm(1,alpha,alpha/3)*d_func(d[i,j])/(sum(N)-N[i])
          if(A[i,j] < 0)
            A[i,j] <- 0
          A[j,i] <- A[i,j]
        }
      }
    }
  }
  return(A)
}

# deterministic discrete
det_disc <- function(So,Io,Ro,A,gamma,SIRout){
  while(sum(Io)>0 & time<tmax){
    deltat <- 1
    time <- time + deltat
    
    StoI <- So*(A%*%Io)
    ItoR <- gamma*Io
    Sn <- So - ceiling(StoI)
    In <- Io + ceiling(StoI) - floor(ItoR)
    Rn <- Ro + floor(ItoR)

    if(any(Sn < 0)){ # control for negative S groups
      for(i in 1:m){
        if(Sn[i] < 0){
          extra <- Sn[i]
          Sn[i] <- 0 # S group stops at 0
          In[i] <- In[i] + extra # remove extra that were already added to I
        }
      }
    }

    SIRout <- rbind(SIRout,c(time,sum(Sn),sum(In),sum(Rn),t(c(Sn,In,Rn))))

    So <- Sn
    Io <- In
    Ro <- Rn
  }
  return(SIRout)
}

# base (no groups) for comparison - deterministic
base_det <- function(bSo,bIo,bRo,alpha,tPop,gamma,baseSIRout,maintime){
  while(sum(Io)>0 & time<tmax & time < maintime){
    deltat <- 1
    time <- time + deltat
    
    bStoI <- bSo*(alpha/tPop)*bIo
    bItoR <- gamma*bIo
    bSn <- bSo - ceiling(bStoI)
    bIn <- bIo + ceiling(bStoI) - floor(bItoR)
    bRn <- bRo + floor(bItoR)

    if(any(bSn < 0)){ # control for negative S groups
      for(i in 1:m){
        if(bSn[i] < 0){
          extra <- bSn[i]
          bSn[i] <- 0  # S group stops at 0
          bIn[i] <- bIn[i] + extra # remove extra that were already added to I
        }
      }
    }

    baseSIRout <- rbind(baseSIRout,c(time,bSn,bIn,bRn))

    bSo <- bSn
    bIo <- bIn
    bRo <- bRn
  }
  return(baseSIRout)
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