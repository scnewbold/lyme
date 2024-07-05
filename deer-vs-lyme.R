#===============================================================================
# deer-vs-lyme.R
#===============================================================================

#-------------------------------------------------------------------------------
# PRELIMINARIES:
#-------------------------------------------------------------------------------
if(TRUE){
  
  # Clear environment to start fresh
  rm(list=setdiff(ls(),lsf.str())) 
  
  # Grab script name for file handling:
  script.name <- basename(rstudioapi::getSourceEditorContext()$path) 
  script.name <- gsub(".R","",script.name)
  
  # Packages
  list.of.packages <- 
    c('dplyr',
      'extrafont',
      'googlesheets4',
      'lubridate',
      'MAd',
      'MASS',
      'metafor',
      'NlcOptim',
      'numDeriv',
      'openxlsx',
      'optimx',
      'rgl',
      'readxl',
      'robumeta',
      'stringr',
      'tidyverse',
      'xtable',
      'tikzDevice'
    )
  
  new.packages <- list.of.packages[!(list.of.packages %in% 
                                       installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  lapply(list.of.packages,function(x){library(x,character.only=TRUE)})
  
  # Clear console:
  cat('\014');
  
  # Clear all plots:
  try(dev.off(dev.list()["RStudioGD"]),silent=TRUE)
  try(dev.off(),silent=TRUE)
  
  # Define paths for file handling:
  this.dir <- dirname(parent.frame(2)$ofile) # source file dir
  setwd(this.dir)                            # set wd to source file dir
  code.path   <- getwd()                     # define code path
  output.path <- getwd()                     # define output path
  input.path  <- getwd()                      # define input path

  # Create output file in working directory:
  date.time     <- gsub(" ","_",Sys.time())
  date.time     <- gsub("-","_",date.time)
  date.time     <- gsub(":","_",date.time)
  out.file.name <- paste(output.path,'/',script.name,'-',date.time,'.out',sep='')
  outfile       <- file.create(out.file.name) 
  
  # WRITE SOURCE FILE TO OUTPUT FILE:
  {
    source.file.name <- paste(code.path,'/',script.name,'.R',sep='')
    
    # Read lines of source file:
    Rscript <- readLines(source.file.name)
    
    # Write lines of source file to output file:
    for(i in 1:length(Rscript)){cat('\n',Rscript[i],file=out.file.name,append=TRUE)}
    cat('\n\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R SCRIPT ABOVE                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|',file=out.file.name,append=TRUE)
    cat('\n| R OUTPUT BELOW                                                            |',file=out.file.name,append=TRUE)
    cat('\n|---------------------------------------------------------------------------|\n',file=out.file.name,append=TRUE)
    
  }
  
}

#-------------------------------------------------------------------------------
# FUNCTIONS:
#-------------------------------------------------------------------------------
{
  pause        <- function(sec){
    
    if(missing(sec)){
      readline(prompt='Paused. Press [enter] to continue.\n')
    }else{
      Sys.sleep(sec)
    }
    
  }
  
  tic          <- function(){
    start.time<<-proc.time()[3]
  }
  
  toc          <- function(){
    elapsed<-(proc.time()[3]-start.time)
    return(sprintf("elapsed time = %0.4f sec",elapsed))
  }

  eq.fn        <- function(tau){
    
    # Computes equilibrium densities of deer (D), adult ticks (A) and infected
    # ticks (II) given a hunting fee (tau)
   
    flag     <- 0
    adjuster <- .1
    
    xh    <- c()
    xh[1] <- 0
    done  <- 0
    i     <- 0
    while(done==0){
      
      i <- i + 1
    
      # Deer:
      D  <- max(0,K*(1-n.h*xh[i]*gam/r))
      
      # Ticks:
      F.L <- 1-exp(-eta.M*M-eta.O*O-eta.D*D)
      f.L <- (eta.M*M*kappa.M+eta.O*O*kappa.O+eta.D*D*kappa.D)/(eta.M*M+eta.O*O+eta.D*D)
      F.N <- 1-exp(-theta.M*M-theta.O*O-theta.D*D)
      f.N <- (theta.M*M*kappa.M+theta.O*O*kappa.O+theta.D*D*kappa.D)/(theta.M*M+theta.O*O+theta.D*D)
      F.A <- 1-exp(-delta.M*M-delta.O*O-delta.D*D)
      f.A <- (delta.M*M*kappa.M+delta.O*O*kappa.O+delta.D*D*kappa.D)/(delta.M*M+delta.O*O+delta.D*D)
      S   <- Beta*(phi*sigma.L*sigma.N*sigma.A*F.L*F.N*F.A-1); S[which(S<=0)]   <- 0
      L   <- phi*S/(1+S/Beta)                                ; L[which(L<=0)]   <- 0
      Nu  <- L *F.L*(1-f.L)*sigma.L                          ; Nu[which(Nu<=0)] <- 0
      Ni  <- L *F.L*f.L    *sigma.L                          ; Ni[which(Ni<=0)] <- 0
      Au  <- Nu*F.N*(1-f.N)*sigma.N                          ; Au[which(Au<=0)] <- 0
      Ai  <- (Nu*f.N + Ni)*F.N*sigma.N                       ; Ai[which(Ai<=0)] <- 0
      II  <- Ni + Ai
      
      # Hunting trip demand:
      # xxh <- adjuster * min(50,exp(a.h+b.h*gam*D-lambda.h*(p.h+tau+rho.h*II*Z))) + (1-adjuster)*xh
      
      xh[i+1] <- adjuster * exp(a.h+b.h*gam*D-lambda.h*(p.h+tau+rho.h*II*Z)) + (1-adjuster)*xh[i]
      
      # print(xh[i+1]); pause()
      
      # Converged?
      if((abs(xh[i+1]-xh[i])/xh[i+1])<0.00001){done <- 1}#else{xh <-xxh}
      
      if(i==100){
        # plot(xh);pause()
        xh[i+1] <- mean(xh[80:100])
        done <- 1
        flag <- 1
      }
    
    }
    
    return(list(D    = D,
                A    = Au+Ai,
                Ai   = Ai,
                Ni   = Ni,
                II   = II,
                xh   = xh[i+1],
                h    = xh[i+1]*gam*D,
                H    = xh[i+1]*gam*D*n.h,
                flag = flag))
     
  }
  
  CS.fn        <- function(tau){
    
    outs  <- eq.fn(tau)
    D     <- outs$D
    A     <- outs$A
    II    <- outs$II
    xh    <- outs$xh
    xk    <- exp(a.k-lambda.k*(p.k+rho.k*II*Z))
    CSh   <- n.h*xh/lambda.h + n.h*xh*tau
    CShh  <- n.h*xh/lambda.h 
    Rev   <- n.h*xh*tau
    CSk   <- n.k*xk/lambda.k
    CShk  <-  CSh+CSk
    return(list(CSh  = CSh,
                CSk  = CSk,
                CShk = CShk,
                Rev  = Rev,
                CShh = CShh))
    
  }
  
  CSh.fn       <- function(tau){
    
    outs <- eq.fn(tau)
    D    <- outs$D
    A    <- outs$A
    II   <- outs$II
    xh   <- outs$xh
    xk   <- exp(a.k-lambda.k*(p.k+rho.k*II*Z))
    CSh  <- n.h*xh/lambda.h + n.h*xh*tau
    CSk  <- n.k*xk/lambda.k
    CShk <- CSh+CSk
    return(CSh)
    
  }
  
  tau.h.opt    <- function(){
    round(optimize(CSh.fn,lower=-1000,upper=1000,maximum=TRUE,tol=0.01)$maximum,2)
  }
  
  tau.hk.opt   <- function(){
    round(optimize(CShk.fn,lower=-1000,upper=1000,maximum=TRUE,tol=0.01)$maximum,2)
  }
  
  B.opt        <- function(tau){
    
    kappa.M0 <- kappa.M
    B.range <- seq(1,100000,length=100)
    V <- matrix(0,length(B.range),1)
    j <- 0
    for(B in B.range){j <- j + 1
      kappa.M <<- kappa.M0 * 1/(psi*B*90)*(1-exp(-psi*B*90))
      outs    <- CS.fn(tau)
      V[j]    <- outs$CShk - price*B*zeta
    }
    B.opt <- B.range[which.max(V)]
    kappa.M <<- kappa.M0
    return(list(B.opt=B.opt,CS=max(V)))
    
  }
  
  tauB.h.opt   <- function(){
    
    tau <- 0
    B   <- 0
    
    done <- FALSE
    while(done==FALSE){
      
      tau2 <- tau.h.opt() 
      B2   <- B.opt(tau2)$B.opt
      if(abs((tau2-tau)/(tau+0.0001))<0.001 & abs((B2-B)/(B+0.0001))<0.001){
        done <- TRUE
      }else{
        tau <- tau2
        B   <- B2
      }
      
    }
    
    return(list(tau=tau,B=B))
    
  }
  
  tauB.hk.opt  <- function(){
    
    tau <- 0
    B   <- 0
    
    done <- FALSE
    while(done==FALSE){
      
      tau2 <- tau.hk.opt() 
      B2   <- B.opt(tau2)$B.opt
      if(abs((tau2-tau)/(tau+0.0001))<0.001 & abs((B2-B)/(B+0.0001))<0.001){
        done <- TRUE
      }else{
        tau <- tau2
        B   <- B2
      }
      
    }
    
    return(list(tau=tau,B=B))
    
  }
  
  CShk.fn      <- function(tau){
    
    outs <- eq.fn(tau)
    D    <- outs$D
    A    <- outs$A
    II   <- outs$II
    xh   <- outs$xh
    xk   <- exp(a.k-lambda.k*(p.k+rho.k*II*Z))
    CSh  <- n.h*xh/lambda.h + n.h*xh*tau
    CSk  <- n.k*xk/lambda.k
    CShk <- CSh+CSk
    return(CShk)
    
  }
  
  IvsD.fn      <- function(D){
    
    # Computes a relationship between equilibrium infected ticks (I) and deer
    # density (D)
    
    # D is a vector
    
    # Ticks:
    F.L <- 1-exp(-eta.M*M-eta.O*O-eta.D*D)
    f.L <- (eta.M*M*kappa.M+eta.O*O*kappa.O+eta.D*D*kappa.D)/(eta.M*M+eta.O*O+eta.D*D)
    F.N <- 1-exp(-theta.M*M-theta.O*O-theta.D*D)
    f.N <- (theta.M*M*kappa.M+theta.O*O*kappa.O+theta.D*D*kappa.D)/(theta.M*M+theta.O*O+theta.D*D)
    F.A <- 1-exp(-delta.M*M-delta.O*O-delta.D*D)
    f.A <- (delta.M*M*kappa.M+delta.O*O*kappa.O+delta.D*D*kappa.D)/(delta.M*M+delta.O*O+delta.D*D)
    S   <- Beta*(phi*sigma.L*sigma.N*sigma.A*F.L*F.N*F.A-1); S[which(S<=0)]   <- 0
    L   <- phi*S/(1+S/Beta)                                ; L[which(L<=0)]   <- 0
    Nu  <- L *F.L*(1-f.L)*sigma.L                          ; Nu[which(Nu<=0)] <- 0
    Ni  <- L *F.L*f.L    *sigma.L                          ; Ni[which(Ni<=0)] <- 0
    Au  <- Nu*F.N*(1-f.N)*sigma.N                          ; Au[which(Au<=0)] <- 0
    Ai  <- (Nu*f.N + Ni)*F.N*sigma.N                       ; Ai[which(Ai<=0)] <- 0
    II  <- Ni + Ai
    
    return(II)
    
  }
  
  tau.fn       <- function(){
    
    # Line search:
    if(TRUE){
      
      passes <- 3
      
      tau.lo <- min(tau.range)
      tau.hi <- max(tau.range)
      
      for(pass in 1:passes){
        
        tau.range <- seq(tau.lo,tau.hi,length=100)
        
        blank   <- matrix(0,length(tau.range),1)
        D.tau   <- blank
        A.tau   <- blank
        II.tau  <- blank
        xh.tau  <- blank
        xk.tau  <- blank
        CSh.tau <- blank
        CSk.tau <- blank
        CS.tau  <- blank
        
        i <- 0
        for(tau in tau.range){i <- i + 1
          outs       <- eq.fn(tau)
          D.tau[i]   <- outs$D
          A.tau[i]   <- outs$A
          II.tau[i]  <- outs$II
          xh.tau[i]  <- outs$xh
          xk.tau[i]  <- exp(a.k-lambda.k*(p.k+rho.k*II.tau[i]*Z))
          CSh.tau[i] <- n.h*xh.tau[i]/lambda.h + n.h*xh.tau[i]*tau
          CSk.tau[i] <- n.k*xk.tau[i]/lambda.k
          CS.tau[i]  <- CSh.tau[i]+CSk.tau[i]
        }
        
        D.0  <- D.tau[which(tau.range==0)]
        II.0 <- II.tau[which(tau.range==0)]
        CS.0 <- CS.tau[which(tau.range==0)]
        
        tau.star.h  <- tau.range[which(CSh.tau==max(CSh.tau))]
        tau.star.hk <- tau.range[which(CS.tau==max(CS.tau))]
        CS          <- max(CS.tau)
        
        tau.lo <- min(c(tau.star.h,tau.star.hk))
        tau.lo <- tau.range[max(which(tau.range<tau.lo))]
        tau.hi <- max(c(tau.star.h,tau.star.hk))
        tau.hi <- tau.range[min(which(tau.range>tau.hi))]
      
      }
      
    }

    # Bisection search:
    if(FALSE){
      
      # Search for tau.star.hk:
      tau.lo      <- min(tau.range)
      tau.hi      <- max(tau.range)
      dCSdtau.lo  <- CS.fn(tau.lo+1)$CS-CS.fn(tau.lo)$CS
      dCSdtau.hi  <- CS.fn(tau.hi+1)$CS-CS.fn(tau.hi)$CS
      done <- 0
      while(done==0){
        tau.md     <- (tau.lo+tau.hi)/2
        dCSdtau.md <- CS.fn(tau.md+1)$CS-CS.fn(tau.md)$CS
        if(dCSdtau.md>0){tau.lo <- tau.md; dCSdtau.lo <- dCSdtau.md}
        if(dCSdtau.md<0){tau.hi <- tau.md; dCSdtau.hi <- dCSdtau.md}
        if((tau.hi-tau.lo)<1){done<-1}
      }
      tau.star.hk <- tau.md
      CS <- CS.fn(tau.md)$CS
      
      # Search for tau.star.h:
      tau.lo      <- min(tau.range)
      tau.hi      <- max(tau.range)
      dCSdtau.lo  <- CS.fn(tau.lo+1)$CSh-CS.fn(tau.lo)$CSh
      dCSdtau.hi  <- CS.fn(tau.hi+1)$CSh-CS.fn(tau.hi)$CSh
      done <- 0
      while(done==0){
        tau.md     <- (tau.lo+tau.hi)/2
        dCSdtau.md <- CS.fn(tau.md+1)$CSh-CS.fn(tau.md)$CSh
        if(dCSdtau.md>0){tau.lo <- tau.md; dCSdtau.lo <- dCSdtau.md}
        if(dCSdtau.md<0){tau.hi <- tau.md; dCSdtau.hi <- dCSdtau.md}
        if((tau.hi-tau.lo)<1){done<-1}
      }
      tau.star.h <- tau.md
      
    }
        
    return(list( tau.star.h  = tau.star.h,
                 tau.star.hk = tau.star.hk,
                 CS          = CS) )
    
  }
  
  drawparms.fn <- function(shape){
    
    if(shape=='uniform'){
      
      parms <<- list()
      
      # exogenous parameters 'given'
      M       <<- runif(1,M_range[1]       ,M_range[3])      ; parms$M       <<- M
      O       <<- runif(1,O_range[1]       ,O_range[3])      ; parms$O       <<- O
      sigma.A <<- runif(1,sigma.A_range[1] ,sigma.A_range[3]); parms$sigma.A <<- sigma.A
      kappa.M <<- runif(1,kappa.M_range[1] ,kappa.M_range[3]); parms$kappa.M <<- kappa.M
      kappa.O <<- runif(1,kappa.O_range[1] ,kappa.O_range[3]); parms$kappa.O <<- kappa.O
      kappa.D <<- runif(1,kappa.D_range[1] ,kappa.D_range[3]); parms$kappa.D <<- kappa.D
      R       <<- runif(1,R_range[1]       ,R_range[3])      ; parms$R       <<- R
      r       <<- runif(1,r_range[1]       ,r_range[3])      ; parms$r       <<- r
      p.h     <<- runif(1,p.h_range[1]     ,p.h_range[3])    ; parms$p.h     <<- p.h
      gam     <<- runif(1,gam_range[1]     ,gam_range[3])    ; parms$gam     <<- gam
      n.k     <<- runif(1,n.k_range[1]     ,n.k_range[3])    ; parms$n.k     <<- n.k
      p.k     <<- runif(1,p.k_range[1]     ,p.k_range[3])    ; parms$p.k     <<- p.k
      Z       <<- runif(1,Z_range[1]       ,Z_range[3])      ; parms$Z       <<- Z

      # endogenous variables, assumed in eq for calibration
      D        <<- runif(1,D_range[1]        ,D_range[3]       ); parms$D        <<- D
      F.L      <<- runif(1,F.L_range[1]      ,F.L_range[3]     ); parms$F.L      <<- F.L
      F.N      <<- runif(1,F.N_range[1]      ,F.N_range[3]     ); parms$F.N      <<- F.N
      F.A      <<- runif(1,F.A_range[1]      ,F.A_range[3]     ); parms$F.A      <<- F.A
      BM.L     <<- runif(1,BM.L_range[1]     ,BM.L_range[3]    ); parms$BM.L     <<- BM.L  
      BM.N     <<- runif(1,BM.N_range[1]     ,BM.N_range[3]    ); parms$BM.N     <<- BM.N
      BM.A     <<- runif(1,BM.A_range[1]     ,BM.A_range[3]    ); parms$BM.A     <<- BM.A
      BO.L     <<- runif(1,BO.L_range[1]     ,BO.L_range[3]    ); parms$BO.L     <<- BO.L    
      BO.N     <<- runif(1,BO.N_range[1]     ,BO.N_range[3]    ); parms$BO.N     <<- BO.N
      BO.A     <<- runif(1,BO.A_range[1]     ,BO.A_range[3]    ); parms$BO.A     <<- BO.A
      BD.L     <<- runif(1,BD.L_range[1]     ,BD.L_range[3]    ); parms$BD.L     <<- BD.L
      BD.N     <<- runif(1,BD.N_range[1]     ,BD.N_range[3]    ); parms$BD.N     <<- BD.N
      BD.A     <<- runif(1,BD.A_range[1]     ,BD.A_range[3]    ); parms$BD.A     <<- BD.A
      e.ph     <<- runif(1,e.ph_range[1]     ,e.ph_range[3]    ); parms$e.ph     <<- e.ph
      e.gamD   <<- runif(1,e.gamD_range[1]   ,e.gamD_range[3]  ); parms$e.gamD   <<- e.gamD
      e.pk     <<- runif(1,e.pk_range[1]     ,e.pk_range[3]    ); parms$e.pk     <<- e.pk
      ell      <<- runif(1,ell_range[1]      ,ell_range[3]     ); parms$ell      <<- ell
      mr.ratio <<- runif(1,mr.ratio_range[1] ,mr.ratio_range[3]); parms$mr.ratio <<- mr.ratio
      x.h      <<- runif(1,x.h_range[1]      ,x.h_range[3]     ); parms$x.h      <<- x.h
      x.k      <<- runif(1,x.k_range[1]      ,x.k_range[3]     ); parms$x.k      <<- x.k
      
      # calibrate remaining exogenous parameters:
      calibrate.fn()
      
    }
    
    if(shape=='triangular'){
        
        parms <<- list()
        
        # exogenous parameters 'given'
        M       <<- rtri(1,M_range[1]       ,M_range[3]       ,M_range[2])      ; parms$M       <<- M
        O       <<- rtri(1,O_range[1]       ,O_range[3]       ,O_range[2])      ; parms$O       <<- O
        sigma.A <<- rtri(1,sigma.A_range[1] ,sigma.A_range[3] ,sigma.A_range[2]); parms$sigma.A <<- sigma.A
        kappa.M <<- rtri(1,kappa.M_range[1] ,kappa.M_range[3] ,kappa.M_range[2]); parms$kappa.M <<- kappa.M
        kappa.O <<- rtri(1,kappa.O_range[1] ,kappa.O_range[3] ,kappa.O_range[2]); parms$kappa.O <<- kappa.O
        kappa.D <<- rtri(1,kappa.D_range[1] ,kappa.D_range[3] ,kappa.D_range[2]); parms$kappa.D <<- kappa.D
        R       <<- rtri(1,R_range[1]       ,R_range[3]       ,R_range[2])      ; parms$R       <<- R
        r       <<- rtri(1,r_range[1]       ,r_range[3]       ,r_range[2])      ; parms$r       <<- r
        p.h     <<- rtri(1,p.h_range[1]     ,p.h_range[3]     ,p.h_range[2])    ; parms$p.h     <<- p.h
        gam     <<- rtri(1,gam_range[1]     ,gam_range[3]     ,gam_range[2])    ; parms$gam     <<- gam
        n.k     <<- rtri(1,n.k_range[1]     ,n.k_range[3]     ,n.k_range[2])    ; parms$n.k     <<- n.k
        p.k     <<- rtri(1,p.k_range[1]     ,p.k_range[3]     ,p.k_range[2])    ; parms$p.k     <<- p.k
        Z       <<- rtri(1,Z_range[1]       ,Z_range[3]       ,Z_range[2])      ; parms$Z       <<- Z
        
        # endogenous variables, assumed in eq for calibration
        D        <<- rtri(1,D_range[1]        ,D_range[3]        ,D_range[2]       ); parms$D        <<- D
        F.L      <<- rtri(1,F.L_range[1]      ,F.L_range[3]      ,F.L_range[2]     ); parms$F.L      <<- F.L
        F.N      <<- rtri(1,F.N_range[1]      ,F.N_range[3]      ,F.N_range[2]     ); parms$F.N      <<- F.N
        F.A      <<- rtri(1,F.A_range[1]      ,F.A_range[3]      ,F.A_range[2]     ); parms$F.A      <<- F.A
        BM.L     <<- rtri(1,BM.L_range[1]     ,BM.L_range[3]     ,BM.L_range[2]    ); parms$BM.L     <<- BM.L  
        BM.N     <<- rtri(1,BM.N_range[1]     ,BM.N_range[3]     ,BM.N_range[2]    ); parms$BM.N     <<- BM.N
        BM.A     <<- rtri(1,BM.A_range[1]     ,BM.A_range[3]     ,BM.A_range[2]    ); parms$BM.A     <<- BM.A
        BO.L     <<- rtri(1,BO.L_range[1]     ,BO.L_range[3]     ,BO.L_range[2]    ); parms$BO.L     <<- BO.L    
        BO.N     <<- rtri(1,BO.N_range[1]     ,BO.N_range[3]     ,BO.N_range[2]    ); parms$BO.N     <<- BO.N
        BO.A     <<- rtri(1,BO.A_range[1]     ,BO.A_range[3]     ,BO.A_range[2]    ); parms$BO.A     <<- BO.A
        BD.L     <<- rtri(1,BD.L_range[1]     ,BD.L_range[3]     ,BD.L_range[2]    ); parms$BD.L     <<- BD.L
        BD.N     <<- rtri(1,BD.N_range[1]     ,BD.N_range[3]     ,BD.N_range[2]    ); parms$BD.N     <<- BD.N
        BD.A     <<- rtri(1,BD.A_range[1]     ,BD.A_range[3]     ,BD.A_range[2]    ); parms$BD.A     <<- BD.A
        e.ph     <<- rtri(1,e.ph_range[1]     ,e.ph_range[3]     ,e.ph_range[2]    ); parms$e.ph     <<- e.ph
        e.gamD   <<- rtri(1,e.gamD_range[1]   ,e.gamD_range[3]   ,e.gamD_range[2]  ); parms$e.gamD   <<- e.gamD
        e.pk     <<- rtri(1,e.pk_range[1]     ,e.pk_range[3]     ,e.pk_range[2]    ); parms$e.pk     <<- e.pk
        ell      <<- rtri(1,ell_range[1]      ,ell_range[3]      ,ell_range[2]     ); parms$ell      <<- ell
        mr.ratio <<- rtri(1,mr.ratio_range[1] ,mr.ratio_range[3] ,mr.ratio_range[2]); parms$mr.ratio <<- mr.ratio
        x.h      <<- rtri(1,x.h_range[1]      ,x.h_range[3]      ,x.h_range[2]     ); parms$x.h      <<- x.h
        x.k      <<- rtri(1,x.k_range[1]      ,x.k_range[3]      ,x.k_range[2]     ); parms$x.k      <<- x.k
        
        
        # calibrate remaining exogenous parameters:
        calibrate.fn()
        
      }
    
    return(parms)
    
  }
  
  calibrate.fn <- function(){
    
    # Fixed model parameters
    eta.M    <<- -log(1-F.L)/(M+O*BO.L/BM.L+D*BD.L/BM.L)
    eta.O    <<- eta.M*BO.L/BM.L
    eta.D    <<- eta.M*BD.L/BM.L
    theta.M  <<- -log(1-F.N)/(M+O*BO.N/BM.N+D*BD.N/BM.N)
    theta.O  <<- theta.M*BO.N/BM.N
    theta.D  <<- theta.M*BD.N/BM.N
    delta.D  <<- -log(1-F.A)/(M*BM.A/BD.A+O*BO.A/BD.A+D)
    delta.M  <<- delta.D*BM.A/BD.A
    delta.O  <<- delta.D*BO.A/BD.A
    
    # Endogenous state variables
    L        <<- -BM.L/eta.M  *log(1-F.L)/F.L
    N        <<- -BM.N/theta.M*log(1-F.N)/F.N
    A        <<- BD.A/delta.D*(1+delta.M*M+delta.O*O+delta.D*D)
    S        <<- A*F.A*sigma.A
    
    # Fixed model parameters
    sigma.L  <<- N/(L*F.L)
    sigma.N  <<- A/(N*F.N)
    phi      <<- R/(sigma.L*sigma.N*sigma.A*F.L*F.N*F.A)
    Beta     <<- L*S/(phi*S-L)
    K        <<- D/(1-mr.ratio)
    n.h      <<- mr.ratio*r/(gam*x.h)
    b.h      <<- e.gamD/(gam*D)
    lambda.h <<- -e.ph/p.h
    lambda.k <<- -e.pk/p.k

    
    #------------------------
    # Dilution scenario:
    # D       <<- D       * 10
    # eta.D   <<- eta.D   * 100 (BD.L+, BM.L-)
    # theta.D <<- theta.D * 10  (BD.N+, BM.N-)
    #------------------------
    
    #------------------------
    # Amplification scenario:
    # lower phi to amplify
    # phi  <<- phi  * .6
    #------------------------
    
    
    #------------------------
    # No Lyme scenario:
    # cost of Lyme set to 0 
    # Z  <<- 0
    #------------------------
    
    # Endogenous state variables
    L    <<- phi*S/(1+S/Beta)
    Ni   <<- L*(1-exp(-eta.M*M-eta.O*O-eta.D*D))*(eta.M*M*kappa.M+eta.O*O*kappa.O+eta.D*D*kappa.D)/(eta.M*M+eta.O*O+eta.D*D)*sigma.L
    Nu   <<- L*(1-exp(-eta.M*M-eta.O*O-eta.D*D))*(eta.M*M*(1-kappa.M)+eta.O*O*(1-kappa.O)+eta.D*D*(1-kappa.D))/(eta.M*M+eta.O*O+eta.D*D)*sigma.L
    Au   <<- Nu*(1-exp(-theta.M*M-theta.O*O-theta.D*D))*(theta.M*M*(1-kappa.M)+theta.O*O*(1-kappa.O)+theta.D*D*(1-kappa.D))/(theta.M*M+theta.O*O+theta.D*D)*sigma.N
    Ai   <<- Nu*(1-exp(-theta.M*M-theta.O*O-theta.D*D))*(theta.M*M*kappa.M+theta.O*O*kappa.O+theta.D*D*kappa.D)/(theta.M*M+theta.O*O+theta.D*D)*sigma.N + 
             Ni*(1-exp(-theta.M*M-theta.O*O-theta.D*D))*sigma.N
    II   <<- Ni+Ai
    
    # Fixed model parameters
    rho.k <<- ell / (x.k*II) 
    rho.h <<- ell / (x.h*II)
    a.k   <<- log(x.k) + lambda.k*(p.k + rho.k*II*Z)
    a.h   <<- log(x.h) + lambda.h*(p.h + rho.h*II*Z) - b.h*gam*D 
    
  }
}

#-------------------------------------------------------------------------------
# MAIN PROGRAM:
#-------------------------------------------------------------------------------

TeXfigs <- TRUE

# SET GIVEN EXOGENOUS PARAMETERS AND ASSUMED EQ VALUES FOR ENDOGENOUS STATE 
# VARIABLES AND CALIBRATE REMAINING EXOGENOUS PARAMETERS:
{
  
  # Sources: https://docs.google.com/spreadsheets/d/1I1UNut_v1I0gkp9xBbSTBUYl8q-2hfBh2fHW5rVsUSQ/edit#gid=0
  
  parms       <- list()
  # exogenous parameters 'given'
  # ranges are:      lo    , primary, hi
  M_range       <- c(1e3   , 8e3    , 25e3)  ; M       <- M_range[2]      ; parms$M       <- M
  O_range       <- c(1e3   , 8e3    , 25e3)  ; O       <- O_range[2]      ; parms$O       <- O
  sigma.A_range <- c(0.2   , 0.5    , 0.7)   ; sigma.A <- sigma.A_range[2]; parms$sigma.A <- sigma.A
  kappa.M_range <- c(0.70  , 0.90   , 1.0)   ; kappa.M <- kappa.M_range[2]; parms$kappa.M <- kappa.M
  kappa.O_range <- c(0.10  , 0.30   , 0.50)  ; kappa.O <- kappa.O_range[2]; parms$kappa.O <- kappa.O
  kappa.D_range <- c(0.00  , 0.025  , 0.05)  ; kappa.D <- kappa.D_range[2]; parms$kappa.D <- kappa.D
  r_range       <- c(0.5   , 0.7    , 0.9)   ; r       <- r_range[2]      ; parms$r       <- r
  p.h_range     <- c(100   , 167.4  , 200)   ; p.h     <- p.h_range[2]    ; parms$p.h     <- p.h
  gam_range     <- c(.00161, .00322 , .00644); gam     <- gam_range[2]    ; parms$gam     <- gam
  n.k_range     <- c(5     , 22     , 100)   ; n.k     <- n.k_range[2]    ; parms$n.k     <- n.k
  p.k_range     <- c(50    , 85     , 100)   ; p.k     <- p.k_range[2]    ; parms$p.k     <- p.k
  Z_range       <- c(9342  , 14344  , 19347) ; Z       <- Z_range[2]      ; parms$Z       <- Z

  # endogenous variables, assumed in eq for calibration
  R_range        <- c(1.5  , 3      , 4.5)  ; R        <- R_range[2]       ; parms$R        <- R 
  D_range        <- c(40   , 80     , 120)  ; D        <- D_range[2]       ; parms$D        <- D
  F.L_range      <- c(0.60 , 0.75   , 0.90) ; F.L      <- F.L_range[2]     ; parms$F.L      <- F.L
  F.N_range      <- c(0.60 , 0.75   , 0.90) ; F.N      <- F.N_range[2]     ; parms$F.N      <- F.N
  F.A_range      <- c(0.60 , 0.75   , 0.90) ; F.A      <- F.A_range[2]     ; parms$F.A      <- F.A
  BM.L_range     <- c(1    , 10     , 30)   ; BM.L     <- BM.L_range[2]    ; parms$BM.L     <- BM.L
  BM.N_range     <- c(.2   , 1      , 5)    ; BM.N     <- BM.N_range[2]    ; parms$BM.N     <- BM.N
  BM.A_range     <- c(0    , 0      , 0)    ; BM.A     <- BM.A_range[2]    ; parms$BM.A     <- BM.A
  BO.L_range     <- c(1    , 10     , 30)   ; BO.L     <- BO.L_range[2]    ; parms$BO.L     <- BO.L
  BO.N_range     <- c(.2   , 1      , 5)    ; BO.N     <- BO.N_range[2]    ; parms$BO.N     <- BO.N
  BO.A_range     <- c(0    , 0      , 0)    ; BO.A     <- BO.A_range[2]    ; parms$BO.A     <- BO.A
  BD.L_range     <- c(200  , 400    , 800)  ; BD.L     <- BD.L_range[2]    ; parms$BD.L     <- BD.L
  BD.N_range     <- c(50   , 150    , 300)  ; BD.N     <- BD.N_range[2]    ; parms$BD.N     <- BD.N
  BD.A_range     <- c(150  , 300    , 600)  ; BD.A     <- BD.A_range[2]    ; parms$BD.A     <- BD.A
  e.ph_range     <- c(-1.5 , -1     , -.5)  ; e.ph     <- e.ph_range[2]    ; parms$e.ph     <- e.ph
  e.gamD_range   <- c(0.5  , 1      , 2.0)  ; e.gamD   <- e.gamD_range[2]  ; parms$e.gamD   <- e.gamD
  e.pk_range     <- c(-4.5 , -3     , -1.5) ; e.pk     <- e.pk_range[2]    ; parms$e.pk     <- e.pk
  ell_range      <- c(.005 , 0.0144 , 0.03) ; ell      <- ell_range[2]     ; parms$ell      <- ell
  mr.ratio_range <- c(0.1  , 0.50   , 0.9)  ; mr.ratio <- mr.ratio_range[2]; parms$mr.ratio <- mr.ratio
  x.h_range      <- c(5    , 9.5    , 15)   ; x.h      <- x.h_range[2]     ;  parms$x.h     <- x.h
  x.k_range      <- c(19   , 36.6   , 76)   ; x.k      <- x.k_range[2]     ; parms$x.k      <- x.k
  
  parms0 <- parms
  
  # calibrate remaining exogenous parameters:
  calibrate.fn()
  
}

tau.range <- seq(-600,600,1)

blank   <- matrix(0,length(tau.range),1)
D.tau   <- blank
T.tau   <- blank
II.tau  <- blank
xh.tau  <- blank
xk.tau  <- blank
CSh.tau <- blank
CSk.tau <- blank
CS.tau  <- blank
H.tau   <- blank

outs <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)

i <- 0
for(tau in tau.range){i <- i + 1
  outs       <- eq.fn(tau)
  D.tau[i]   <- outs$D  # Deer
  II.tau[i]  <- outs$II  # Infected ticks
  xh.tau[i]  <- outs$xh # Trips per hunter
  xk.tau[i]  <- exp(a.k-lambda.k*(p.k+rho.k*II.tau[i]*Z)) # trips per hiker
  CSh.tau[i] <- n.h*xh.tau[i]/lambda.h + n.h*xh.tau[i]*tau # hunters consumer surplus (+ fee revenues)
  CSk.tau[i] <- n.k*xk.tau[i]/lambda.k # hikers consumer surplus
  CS.tau[i]  <- CSh.tau[i]+CSk.tau[i]  # total surplus
  H.tau[i]   <- xh.tau[i]*gam*D.tau[i]*n.h
}

D.0  <- D.tau[which.min(tau.range^2)]
II.0 <- II.tau[which.min(tau.range^2)]
CS.0 <- CS.tau[which.min(tau.range^2)]

tau.star.h  <- optimize(CSh.fn ,lower=-1000,upper=1000,maximum=TRUE,tol=0.01)$maximum
tau.star.hk <- optimize(CShk.fn,lower=-1000,upper=1000,maximum=TRUE,tol=0.01)$maximum

# GRAPHS:

# panel a: D/D0 and I vs tau
{
  show <- which(D.tau>0)
  lft  <- -600
  rht  <- 600
  top  <- 2
  btm  <- 0
  
  tikz(paste(output.path,'/panel-a.tex',sep=''),
       width      = 4,
       height     = 4,
       pointsize  = 10,
       standAlone = TRUE)
  
  par(mar=c(5,4,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
  plot(tau.range[show],D.tau[show]/D.0 , 
       type     = 'l',
       col      = 'black',
       xlim     = c(lft,rht),
       ylim     = c(btm,top),
       main     = '',
       xlab     = '',
       ylab     = '',
       axes     = FALSE,
       family   = 'Bookman',
       cex      = .8,
       cex.axis = .8,
       cex.lab  = .8)
  lines(tau.range[show] , II.tau[show]/II.0 , col='black',lty=2)
  lines(tau.range , matrix(1,length(tau.range),1) , lty=3)
  axis(1,las=1,line=0,labels=seq(-600,600,100),at=seq(-600,600,100),cex.axis=.8) # Draw x axis
  axis(2,las=0,line=0,labels=seq(0,2,.25),at=seq(0,2,.25),cex.axis=.8) # Draw y axis
  mtext('Relative density'         ,side=2, line=2, cex=.8, las=0) # y-axis label
  mtext('Hunting trip fee, $\\tau$',side=1, line=2, cex=.8, las=1) # x-axis
  legend(lft+1*(rht-lft),btm+0.2*(top-btm),
         legend=c('Deer','Infected ticks'),col=c('black','black'),
         lty=c(1,2),xjust=1,cex=.8)
  
  dev.off()

  
}

# panel b: CS, CSh, and Csk vs tau
{
  show <- which(D.tau>0)
  lft  <- min(tau.range)
  rht  <- max(tau.range)
  top  <- max(CS.tau[show])
  btm  <- min(c(CSh.tau[show],CSk.tau[show],CS.tau[show]))

  tikz(paste(output.path,'/panel-b.tex',sep=''),
       width      = 4,
       height     = 4,
       pointsize  = 10,
       standAlone = TRUE)
  
  par(mar=c(5,4,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
  plot(tau.range[show],CS.tau[show], 
       type     = 'l',
       col      = 'black',
       xlim     = c(lft,rht),
       ylim     = c(btm,top),
       main     = '',
       xlab     = '',
       ylab     = '',
       axes     = FALSE,
       family   = 'Bookman',
       cex      = .8,
       cex.axis = .8,
       cex.lab  = .8)
  
  lines(tau.range[show],CSh.tau[show],col='red')
  lines(tau.range[show],CSk.tau[show],col='blue')
  
  lines(c(tau.star.hk,tau.star.hk),c(btm,max(CS.tau )),lty=3,col='black')
  lines(c(tau.star.h,tau.star.h)  ,c(btm,max(CSh.tau)),lty=3,col='red')
  
  legend(300,0,legend=c('$CS$','$CS_h$','$CS_k$'),
         col=c('black','red','blue'),lty=c(1,1,1),cex=.8)
  
  axis(1,las=1,line=0,labels=seq(-600,600,100),at=seq(-600,600,100),cex.axis=.8) # Draw x axis
  axis(2,las=0,line=0,cex.axis=.8) # Draw y axis
  mtext('Consumer surplus [\\$]'         ,side=2, line=2, cex=.8, las=0) # y-axis label
  mtext('$\\tau$, Hunting trip fee [\\$]',side=1, line=2, cex=.8, las=1) # x-axis
  
  dev.off()
  
}

# panel c: I vs D
{
  DD <- seq(0,200,1)
  II <- IvsD.fn(DD)
  
  D.0  <- D.tau[which.min(tau.range^2)]
  D.h  <- D.tau[which.min((tau.range-tau.star.h )^2)]
  D.hk <- D.tau[which.min((tau.range-tau.star.hk)^2)]
  
  plot(DD,II,type='l',col='black',ylim=c(0,max(II)*1.1))
  lines(c(K,K),      c(0,II[which.min((DD-K)^2)]),lty=1)
  lines(c(D.0,D.0),  c(0,II[which.min((DD-D.0)^2)]),lty=3)
  lines(c(D.hk,D.hk),c(0,II.tau[which.min((tau.range-tau.star.hk)^2)]),lty=3,col='green')
  lines(c(D.h,D.h),  c(0,II.tau[which.min((tau.range-tau.star.h )^2)]),lty=3,col='brown')
  
  lft  <- min(DD)
  rht  <- max(DD)
  top  <- max(II)*1.1
  btm  <- min(II)*0.9 * 0
  # top  <- 18000
  # btm  <- 8000
  
  tikz(paste(output.path,'/panel-c.tex',sep=''),
       width      = 4.5,
       height     = 6.5,
       pointsize  = 12,
       standAlone = TRUE)
  plot(DD,II,
       type     = 'l',
       col      = 'black',
       xlim     = c(lft,rht),
       ylim     = c(btm,top),
       main     = '',
       xlab     = '',
       ylab     = '',
       axes     = FALSE,
       family   = 'Bookman',
       cex      = 1.25,
       cex.axis = 1.25,
       cex.lab  = 1.25)
  lines(c(D.0,D.0),  c(0,II.tau[which.min(tau.range^2)]),lty=3)
  lines(c(D.h,D.h),  c(0,II.tau[which.min((tau.range-tau.star.h )^2)]),lty=3)
  lines(c(D.hk,D.hk),c(0,II.tau[which.min((tau.range-tau.star.hk)^2)]),lty=3)
  axis(1,las=1) # Draw x axis
  axis(2,las=0) # Draw y axis
  mtext('Infected tick density',side=2, line=2.5, cex.lab=1.2, las=0) # y-axis label
  mtext('Deer density' ,side=1, line=2.5, cex.lab=1.2, las=1) # x-axis
  dev.off()
  
}

# CREATE TABLE OF BENCHMARK RESULTS:
if(TRUE){
  # Open access
  D.open    <- eq.fn(0)$D
  H.open    <- eq.fn(0)$H
  Ni.open   <- eq.fn(0)$Ni
  Ai.open   <- eq.fn(0)$Ai
  II.open   <- eq.fn(0)$II
  CShh.open <- CS.fn(0)$CShh
  Rev.open  <- CS.fn(0)$Rev
  CShk.open <- CS.fn(0)$CShk
  
  
  # MSY
  H.tau <- matrix(0,length(tau.range),1)
  j <- 0
  for(tau in tau.range){j <- j + 1
    H.tau[j] <- eq.fn(tau)$H
  }
  tau.MSY  <- tau.range[which.max(H.tau)]
  
  D.MSY    <- eq.fn(tau.MSY)$D
  H.MSY    <- eq.fn(tau.MSY)$H
  Ni.MSY   <- eq.fn(tau.MSY)$Ni
  Ai.MSY   <- eq.fn(tau.MSY)$Ai
  II.MSY   <- eq.fn(tau.MSY)$II
  CShh.MSY <- CS.fn(tau.MSY)$CShh
  Rev.MSY  <- CS.fn(tau.MSY)$Rev
  CShk.MSY <- CS.fn(tau.MSY)$CShk
  
  # Max CSh + Rev
  outs       <- optimize(CSh.fn,c(0,200),maximum=TRUE)
  
  tau.CSh  <- outs$maximum
  D.CSh    <- eq.fn(tau.CSh)$D
  H.CSh    <- eq.fn(tau.CSh)$H
  Ni.CSh   <- eq.fn(tau.CSh)$Ni
  Ai.CSh   <- eq.fn(tau.CSh)$Ai
  II.CSh   <- eq.fn(tau.CSh)$II
  CShh.CSh <- CS.fn(tau.CSh)$CShh
  Rev.CSh  <- CS.fn(tau.CSh)$Rev
  CShk.CSh <- CS.fn(tau.CSh)$CShk
  
  
  # Max CShk 
  outs      <- optimize(CShk.fn,c(0,200),maximum=TRUE)
  
  tau.CShk  <- outs$maximum
  D.CShk    <- eq.fn(tau.CShk)$D
  H.CShk    <- eq.fn(tau.CShk)$H
  Ni.CShk   <- eq.fn(tau.CShk)$Ni
  Ai.CShk   <- eq.fn(tau.CShk)$Ai
  II.CShk   <- eq.fn(tau.CShk)$II
  CShh.CShk <- CS.fn(tau.CShk)$CShh
  Rev.CShk  <- CS.fn(tau.CShk)$Rev
  CShk.CShk <- CS.fn(tau.CShk)$CShk
  
  if(TRUE){
    cat('\014')
    cat('BENCHMARK RESULTS TABLE:\n\n')
    cat('Management regime            tau     D     H     I   CSh  CShk\n')
    cat('-------------------------- ----- ----- ----- ----- ----- -----\n')
    cat(sprintf('Open access                & %7.1f & %7.2f & %7.2f & %7.3e & %7.3e & %7.3e\\\\\n',0,D.open,H.open,Ni.open+Ai.open,CShh.open+Rev.open,CShk.open))
    cat(sprintf('Max harvest                & %7.1f & %7.2f & %7.2f & %7.3e & %7.3e & %7.3e\\\\\n',tau.MSY,D.MSY,H.MSY,Ni.MSY+Ai.MSY,CShh.MSY+Rev.MSY,CShk.MSY))
    cat(sprintf('Max hunter surplus         & %7.1f & %7.2f & %7.2f & %7.3e & %7.3e & %7.3e\\\\\n',tau.CSh,D.CSh,H.CSh,Ni.CSh+Ai.CSh,CShh.CSh+Rev.CSh,CShk.CSh))   
    cat(sprintf('Max hunter + hiker surplus & %7.1f & %7.2f & %7.2f & %7.3e & %7.3e & %7.3e\\\\\n',tau.CShk,D.CShk,H.CShk,Ni.CShk+Ai.CShk,CShh.CShk+Rev.CShk,CShk.CShk))
  }

  tab <- as.table(rbind(c(0,       D.open, H.open, Ni.open, Ai.open,  CShh.open, Rev.open, CShk.open), 
                        c(tau.MSY, D.MSY,  H.MSY,  Ni.MSY,  Ai.MSY,   CShh.MSY,  Rev.MSY,  CShk.MSY),
                        c(tau.CSh, D.CSh,  H.CSh,  Ni.CSh,  Ai.CSh,   CShh.CSh,  Rev.CSh,  CShk.CSh), 
                        c(tau.CShk,D.CShk, H.CShk, Ni.CShk, Ai.CShk,  CShh.CShk, Rev.CShk, CShk.CShk)))
  
  dimnames(tab) <- list(max = c("Open Access", "Max H", "Max CShR", "Max CShRk"),
                        " " = c("tau", "Deer", "Total Harvest", "Infected Nymphs", "Infected Adults", "CSh", "Rev", "CShRk"))
  tab
  tab <- round(tab, digits=1)
  cat('\nBenchmark results table:\n',file=out.file.name,append=TRUE)
  print(xtable(tab, digits=1),file=out.file.name,append=TRUE)
}

# BENCHMARK GRAPH:
if(TRUE){
  
  tikz(paste(output.path,'/fig-bench.tex',sep=''),
       width      = 5,
       height     = 5,
       pointsize  = 12,
       standAlone = TRUE)
  
  H.tau.0 <- H.tau*2000
  
  how  <- which(D.tau > 0)
  lft  <- min(tau.range) * 0 - 200
  rht  <- max(tau.range) * 0 + 300
  top  <- max(c(CSh.tau[show], CSk.tau[show], CS.tau[show])) * 1.1 * 0 + 80000
  btm  <- min(c(CSh.tau[show], CSk.tau[show], CS.tau[show])) * 0
  
  par(mar=c(3,3,0.5,3.1)) # bottom, left, top, right
  
  # Plot the main data
  plot(tau.range[show], CS.tau[show], 
       type     = 'l',
       col      = 'black',
       xlim     = c(lft, rht),
       ylim     = c(btm, top),
       main     = '',
       xlab     = '',
       ylab     = '',
       axes     = FALSE,
       family   = 'Bookman',
       cex      = .8,
       cex.axis = .8,
       cex.lab  = .8)
  
  lines(tau.range[show], CSk.tau[show], col='magenta')
  lines(tau.range[show], H.tau.0[show], col='blue')
  lines(tau.range[show], CSh.tau[show], col='red')
  
  tau.star.k <- tau.range[max(which(CSk.tau==max(CSk.tau)))]
  
  lines(c(tau.star.hk, tau.star.hk), c(btm, max(CS.tau)), lty=3, col='black')
  lines(c(tau.star.k , tau.star.k ), c(0,max(CSk.tau)),lty=3,col='magenta')
  lines(c(tau.MSY, tau.MSY), c(btm, max(H.tau.0)), lty=3, col='blue')
  lines(c(tau.star.h, tau.star.h), c(btm, max(CSh.tau)), lty=3, col='red')
  
  text(tau.star.k  ,0,'$\\tau_k$'    ,pos=2)
  text(tau.MSY     ,0,'$\\tau_{msy}$',pos=2)
  text(tau.star.h  ,0,'$\\tau_h$'    ,pos=4)
  text(tau.star.hk ,0,'$\\tau_{hk}$' ,pos=2)
  
  # Draw x-axis and y-axis
  axis(1, las=1, line=0, labels=seq(-600, 600, 100), at=seq(-600, 600, 100), cex.axis=.8) # Draw x axis
  axis(2, las=0, line=0, labels=seq(0,70,10), at=seq(0,70000,10000), cex.axis=.8) # Draw y axis
  
  # Add y-axis label
  mtext('Consumer surplus [1,000\\$/yr/mi$^2$]', side=2, line=2, cex=.8, las=0) # y-axis label
  mtext('$\\tau$, Hunting fee [\\$/trip]', side=1, line=2, cex=.8, las=1) # x-axis
  
  # Add another y-axis on the right side
  axis(4, ylim = c(btm, top), labels=seq(0,35,5), at=seq(0,70000,10000), las = 0, cex.axis = .8)
  mtext('Harvest [/yr/mi$^2$]', side = 4, line = 2, cex = .8, las = 0)
  
  # Add another curve line for the new y-axis
  lines(tau.range[show], H.tau.0[show], col='blue')
  
  # Add legend
  legend(-200,70e3, legend=c('Total surplus', 'Hiker surplus', 'Hunter surplus' , 'Harvest'),
         col=c('black', 'magenta', 'red','blue'), lty=c(1, 1, 1,1), cex=.8)
  
  dev.off()
}

# SENSITIVITY ANALYSIS -> PARAMETER PLOTS:
if(TRUE){
  
  Texfigs <- TRUE
  
  # (m/r) hunting mortality / deer intrinsic growth ratio [endogenous]
  {
    mr.ratio0       <- mr.ratio
    mr.ratio.range  <- seq(mr.ratio_range[1],mr.ratio_range[3],length=100)
    tau.h.mr.ratio  <- matrix(0,length(mr.ratio.range),1)
    tau.hk.mr.ratio <- matrix(0,length(mr.ratio.range),1)
    i <- 0
    for(mr.ratio in mr.ratio.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on m/r =',sprintf('%-.2f\n',mr.ratio))
      tau.h.mr.ratio[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.mr.ratio[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.mr.ratio,tau.h.mr.ratio))*1.05
    bot <- min(c(tau.hk.mr.ratio,tau.h.mr.ratio))*0.95
    plot( mr.ratio.range,tau.hk.mr.ratio,type='l',col='black',ylim=c(bot,top))
    lines(mr.ratio.range,tau.h.mr.ratio ,lty=2   ,col='black')
    lines(c(mr.ratio0,mr.ratio0),c(0,tau.star.h),lty=3)
    mr.ratio <- mr.ratio0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-mr.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( mr.ratio.range,tau.hk.mr.ratio,type='l',col='black',ylim=c(bot,top),main='$m/r$')
      lines(mr.ratio.range,tau.h.mr.ratio ,lty=2   ,col='black')
      lines(c(mr.ratio0,mr.ratio0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (n.k) hiker density [exogenous]
  {
    n.k0       <- n.k
    n.k.range  <- seq(n.k_range[1],n.k_range[3],length=100)
    tau.h.n.k  <- matrix(0,length(n.k.range),1)
    tau.hk.n.k <- matrix(0,length(n.k.range),1)
    i <- 0
    for(n.k in n.k.range){i <- i + 1
      cat('\014');cat('Working on n.k =',sprintf('%-.2f\n',n.k))
      tau.h.n.k[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.n.k[i] <- optimize(CShk.fn,lower=-1000,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.n.k,tau.h.n.k))*1.05
    bot <- min(c(tau.hk.n.k,tau.h.n.k))*0.95
    plot( n.k.range,tau.hk.n.k ,type='l',col='black',ylim=c(bot,top))
    lines(n.k.range,tau.h.n.k  ,lty=2   ,col='black')
    lines(c(n.k0,n.k0),c(bot,tau.star.h),lty=3)
    n.k <- n.k0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-nk.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( n.k.range,tau.hk.n.k,type='l',col='black',ylim=c(bot,top),main='$n_k$')
      lines(n.k.range,tau.h.n.k ,lty=2   ,col='black')
      lines(c(n.k0,n.k0),c(bot,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (R) tick intrinsic growth rate per generation [endogenous]
  {
    R0       <- R
    R.range  <- seq(R_range[1],R_range[3],length=100)
    tau.h.R  <- matrix(0,length(R.range),1)
    tau.hk.R <- matrix(0,length(R.range),1)
    i <- 0
    for(R in R.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on R =',sprintf('%-.2f\n',R))
      tau.h.R[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.R[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.R,tau.h.R))*1.05
    bot <- min(c(tau.hk.R,tau.h.R))*0.95
    plot( R.range,tau.hk.R,type='l',col='black',ylim=c(bot,top))
    lines(R.range,tau.h.R ,lty=2   ,col='black')
    lines(c(R0,R0),c(0,tau.star.h),lty=3)
    R <- R0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-R.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( R.range,tau.hk.R,type='l',col='black',ylim=c(bot,top),main='$R$')
      lines(R.range,tau.h.R ,lty=2   ,col='black')
      lines(c(R0,R0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (Z) Lyme disease cost [exogenous]
  {
    Z0       <- Z
    Z.range  <- seq(Z_range[1],Z_range[3],length=100)
    tau.h.Z  <- matrix(0,length(Z.range),1)
    tau.hk.Z <- matrix(0,length(Z.range),1)
    i <- 0
    for(Z in Z.range){i <- i + 1
    cat('\014');cat('Working on Z =',sprintf('%-.2f\n',Z))
    tau.h.Z[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    tau.hk.Z[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.Z,tau.h.Z))*1.05
    bot <- min(c(tau.hk.Z,tau.h.Z))*0.95
    plot( Z.range,tau.hk.Z,type='l',col='black',ylim=c(bot,top))
    lines(Z.range,tau.h.Z ,lty=2   ,col='black')
    lines(c(Z0,Z0),c(0,tau.star.h),lty=3)
    Z <- Z0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-Z.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( Z.range,tau.hk.Z,type='l',col='black',ylim=c(bot,top),main='$Z$')
      lines(Z.range,tau.h.Z ,lty=2   ,col='black')
      lines(c(Z0,Z0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (e.gamD) hunting demand expected catch elasticity [endogenous]
  {
    e.gamD0       <- e.gamD
    e.gamD.range  <- seq(e.gamD_range[1],e.gamD_range[3],length=100)
    tau.h.e.gamD  <- matrix(0,length(e.gamD.range),1)
    tau.hk.e.gamD <- matrix(0,length(e.gamD.range),1)
    i <- 0
    for(e.gamD in e.gamD.range){i <- i + 1
      cat('\014');cat('Working on e.gamD =',sprintf('%-.2f\n',e.gamD))
      calibrate.fn()
      tau.h.e.gamD[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.e.gamD[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.e.gamD,tau.h.e.gamD))*1.05
    bot <- min(c(tau.hk.e.gamD,tau.h.e.gamD))*0.95
    plot( e.gamD.range,tau.hk.e.gamD,type='l',col='black',ylim=c(bot,top))
    lines(e.gamD.range,tau.h.e.gamD ,lty=2   ,col='black')
    lines(c(e.gamD0,e.gamD0),c(0,tau.star.h),lty=3)
    e.gamD <- e.gamD0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-egamD.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( e.gamD.range,tau.hk.e.gamD,type='l',col='black',ylim=c(bot,top),main='$\\varepsilon_{\\gamma D}$')
      lines(e.gamD.range,tau.h.e.gamD ,lty=2   ,col='black')
      lines(c(e.gamD0,e.gamD0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (O) other host abundance [exogenous]
  {
    O0       <- O
    O.range  <- seq(O_range[1],O_range[3],length=100)
    tau.h.O  <- matrix(0,length(O.range),1)
    tau.hk.O <- matrix(0,length(O.range),1)
    i <- 0
    for(O in O.range){i <- i + 1
    cat('\014');cat('Working on O =',sprintf('%-.1f\n',O))
    tau.h.O[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    tau.hk.O[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.O,tau.h.O))*1.05
    bot <- min(c(tau.hk.O,tau.h.O))*0.95
    plot( O.range,tau.hk.O,type='l',col='black',ylim=c(bot,top))
    lines(O.range,tau.h.O ,lty=2   ,col='black')
    lines(c(O0,O0),c(0,tau.star.h),lty=3)
    O <- O0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-O.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( O.range,tau.hk.O,type='l',col='black',ylim=c(bot,top),main='$O$')
      lines(O.range,tau.h.O ,lty=2   ,col='black')
      lines(c(O0,O0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (p.k) price per hiking trip [exogenous]
  {
    p.k0       <- p.k
    p.k.range  <- seq(p.k_range[1],p.k_range[3],length=100)
    tau.h.p.k  <- matrix(0,length(p.k.range),1)
    tau.hk.p.k <- matrix(0,length(p.k.range),1)
    i <- 0
    for(p.k in p.k.range){i <- i + 1
      cat('\014');cat('Working on p.k =',sprintf('%-.1f\n',p.k))
      tau.h.p.k[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.p.k[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.p.k,tau.h.p.k))*1.05
    bot <- min(c(tau.hk.p.k,tau.h.p.k))*0.95
    plot( p.k.range,tau.hk.p.k,type='l',col='black',ylim=c(bot,top))
    lines(p.k.range,tau.h.p.k ,lty=2   ,col='black')
    lines(c(p.k0,p.k0),c(0,tau.star.h),lty=3)
    p.k <- p.k0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-pk.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( p.k.range,tau.hk.p.k,type='l',col='black',ylim=c(bot,top),main='$p_k$')
      lines(p.k.range,tau.h.p.k ,lty=2   ,col='black')
      lines(c(p.k0,p.k0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (ell) hiker annual risk of Lyme [endogenous]
  {
    ell0       <- ell
    ell.range  <- seq(ell_range[1],ell_range[3],length=100)
    tau.h.ell  <- matrix(0,length(ell.range),1)
    tau.hk.ell <- matrix(0,length(ell.range),1)
    i <- 0
    for(ell in ell.range){i <- i + 1
      cat('\014');cat('Working on ell =',sprintf('%-.1f\n',ell))
      calibrate.fn()
      tau.h.ell[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.ell[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.ell,tau.h.ell))*1.05
    bot <- min(c(tau.hk.ell,tau.h.ell))*0.95
    plot( ell.range,tau.hk.ell,type='l',col='black',ylim=c(bot,top))
    lines(ell.range,tau.h.ell ,lty=2   ,col='black')
    lines(c(ell0,ell0),c(0,tau.star.h),lty=3)
    ell <- ell0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-ell.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( ell.range,tau.hk.ell,type='l',col='black',ylim=c(bot,top),main='$\\ell$')
      lines(ell.range,tau.h.ell ,lty=2   ,col='black')
      lines(c(ell0,ell0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (sigma.A) fed adult tick survival rate [exogenous]
  {
    sigma.A0       <- sigma.A
    sigma.A.range  <- seq(sigma.A_range[1],sigma.A_range[3],length=100)
    tau.h.sigma.A  <- matrix(0,length(sigma.A.range),1)
    tau.hk.sigma.A <- matrix(0,length(sigma.A.range),1)
    i <- 0
    for(sigma.A in sigma.A.range){i <- i + 1
      cat('\014');cat('Working on sigma.A =',sprintf('%-.1f\n',sigma.A))
      tau.h.sigma.A[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.sigma.A[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.sigma.A,tau.h.sigma.A))*1.05
    bot <- min(c(tau.hk.sigma.A,tau.h.sigma.A))*0.95
    plot( sigma.A.range,tau.hk.sigma.A,type='l',col='black',ylim=c(bot,top))
    lines(sigma.A.range,tau.h.sigma.A ,lty=2   ,col='black')
    lines(c(sigma.A0,sigma.A0),c(0,tau.star.h),lty=3)
    sigma.A <- sigma.A0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-sigmaA.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( sigma.A.range,tau.hk.sigma.A,type='l',col='black',ylim=c(bot,top),main='$\\sigma_A$')
      lines(sigma.A.range,tau.h.sigma.A ,lty=2   ,col='black')
      lines(c(sigma.A0,sigma.A0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (e.ph) hunting demand price elasticity [endogenous]
  {
    e.ph0       <- e.ph
    e.ph.range  <- seq(e.ph_range[1],e.ph_range[3],length=100)
    tau.h.e.ph  <- matrix(0,length(e.ph.range),1)
    tau.hk.e.ph <- matrix(0,length(e.ph.range),1)
    i <- 0
    for(e.ph in e.ph.range){i <- i + 1
    calibrate.fn()
      cat('\014');cat('Working on e.ph =',sprintf('%-.2f\n',e.ph))
      tau.h.e.ph[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.e.ph[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.e.ph,tau.h.e.ph))*1.05
    bot <- min(c(tau.hk.e.ph,tau.h.e.ph))*0.95
    plot( e.ph.range,tau.hk.e.ph,type='l',col='black',ylim=c(bot,top))
    lines(e.ph.range,tau.h.e.ph ,lty=2   ,col='black')
    lines(c(e.ph0,e.ph0),c(0,tau.star.h),lty=3)
    e.ph <- e.ph0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-eph.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( e.ph.range,tau.hk.e.ph,type='l',col='black',ylim=c(bot,top),main='$\\varepsilon_{p_h}$')
      lines(e.ph.range,tau.h.e.ph ,lty=2   ,col='black')
      lines(c(e.ph0,e.ph0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (r) deer intrinsic growth rate [exogenous]
  {
    r0       <- r
    r.range  <- seq(r_range[1],r_range[3],length=100)
    tau.h.r  <- matrix(0,length(r.range),1)
    tau.hk.r <- matrix(0,length(r.range),1)
    i <- 0
    for(r in r.range){i <- i + 1
      cat('\014');cat('Working on r =',sprintf('%-.1f\n',r))
      tau.h.r[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.r[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.r,tau.h.r))*1.05
    bot <- min(c(tau.hk.r,tau.h.r))*0.95
    plot( r.range,tau.hk.r,type='l',col='black',ylim=c(bot,top))
    lines(r.range,tau.h.r ,lty=2   ,col='black')
    lines(c(r0,r0),c(0,tau.star.h),lty=3)
    r <- r0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-rr.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( r.range,tau.hk.r,type='l',col='black',ylim=c(bot,top),main='$r$')
      lines(r.range,tau.h.r ,lty=2   ,col='black')
      lines(c(r0,r0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (F.A) host finding success of adult ticks [endogenous]
  {
    F.A0       <- F.A
    F.A.range  <- seq(F.A_range[1],F.A_range[3],length=100)
    tau.h.F.A  <- matrix(0,length(F.A.range),1)
    tau.hk.F.A <- matrix(0,length(F.A.range),1)
    i <- 0
    for(F.A in F.A.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on F.A =',sprintf('%-.2f\n',F.A))
      tau.h.F.A[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.F.A[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.F.A,tau.h.F.A))*1.05
    bot <- min(c(tau.hk.F.A,tau.h.F.A))*0.95
    plot( F.A.range,tau.hk.F.A,type='l',col='black',ylim=c(bot,top))
    lines(F.A.range,tau.h.F.A ,lty=2   ,col='black')
    lines(c(F.A0,F.A0),c(0,tau.star.h),lty=3)
    F.A <- F.A0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-FA.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( F.A.range,tau.hk.F.A,type='l',col='black',ylim=c(bot,top),main='$F_A$')
      lines(F.A.range,tau.h.F.A ,lty=2   ,col='black')
      lines(c(F.A0,F.A0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BM.N) body burden of nymphs on mice [endogenous]
  {
    BM.N0       <- BM.N
    BM.N.range  <- seq(BM.N_range[1],BM.N_range[3],length=100)
    tau.h.BM.N  <- matrix(0,length(BM.N.range),1)
    tau.hk.BM.N <- matrix(0,length(BM.N.range),1)
    i <- 0
    for(BM.N in BM.N.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BM.N =',sprintf('%-.2f\n',BM.N))
      tau.h.BM.N[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BM.N[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BM.N,tau.h.BM.N))*1.05
    bot <- min(c(tau.hk.BM.N,tau.h.BM.N))*0.95
    plot( BM.N.range,tau.hk.BM.N,type='l',col='black',ylim=c(bot,top))
    lines(BM.N.range,tau.h.BM.N ,lty=2   ,col='black')
    lines(c(BM.N0,BM.N0),c(0,tau.star.h),lty=3)
    BM.N <- BM.N0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BMN.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BM.N.range,tau.hk.BM.N,type='l',col='black',ylim=c(bot,top),main='$B_N^M$')
      lines(BM.N.range,tau.h.BM.N ,lty=2   ,col='black')
      lines(c(BM.N0,BM.N0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (D) deer density [endogenous]
  {
    D0       <- D
    D.range  <- seq(D_range[1],D_range[3],length=100)
    tau.h.D  <- matrix(0,length(D.range),1)
    tau.hk.D <- matrix(0,length(D.range),1)
    i <- 0
    for(D in D.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on D =',sprintf('%-.2f\n',D))
      tau.h.D[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.D[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.D,tau.h.D))*1.05
    bot <- min(c(tau.hk.D,tau.h.D))*0.95
    plot( D.range,tau.hk.D,type='l',col='black',ylim=c(bot,top))
    lines(D.range,tau.h.D ,lty=2   ,col='black')
    lines(c(D0,D0),c(0,tau.star.h),lty=3)
    D <- D0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-D.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( D.range,tau.hk.D,type='l',col='black',ylim=c(bot,top),main='$D$')
      lines(D.range,tau.h.D ,lty=2   ,col='black')
      lines(c(D0,D0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BO.N) body burden of nymphs on other hosts [endogenous]
  {
    BO.N0       <- BO.N
    BO.N.range  <- seq(BO.N_range[1],BO.N_range[3],length=100)
    tau.h.BO.N  <- matrix(0,length(BO.N.range),1)
    tau.hk.BO.N <- matrix(0,length(BO.N.range),1)
    i <- 0
    for(BO.N in BO.N.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BO.N =',sprintf('%-.2f\n',BO.N))
      tau.h.BO.N[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BO.N[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BO.N,tau.h.BO.N))*1.05
    bot <- min(c(tau.hk.BO.N,tau.h.BO.N))*0.95
    plot( BO.N.range,tau.hk.BO.N,type='l',col='black',ylim=c(bot,top))
    lines(BO.N.range,tau.h.BO.N ,lty=2   ,col='black')
    lines(c(BO.N0,BO.N0),c(0,tau.star.h),lty=3)
    BO.N <- BO.N0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BON.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BO.N.range,tau.hk.BO.N,type='l',col='black',ylim=c(bot,top),main='$B_N^O$')
      lines(BO.N.range,tau.h.BO.N ,lty=2   ,col='black')
      lines(c(BO.N0,BO.N0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (F.N) host finding success of nymphal ticks [endogenous]
  {
    F.N0       <- F.N
    F.N.range  <- seq(F.N_range[1],F.N_range[3],length=100)
    tau.h.F.N  <- matrix(0,length(F.N.range),1)
    tau.hk.F.N <- matrix(0,length(F.N.range),1)
    i <- 0
    for(F.N in F.N.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on F.N =',sprintf('%-.2f\n',F.N))
      tau.h.F.N[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.F.N[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.F.N,tau.h.F.N))*1.05
    bot <- min(c(tau.hk.F.N,tau.h.F.N))*0.95
    plot( F.N.range,tau.hk.F.N,type='l',col='black',ylim=c(bot,top))
    lines(F.N.range,tau.h.F.N ,lty=2   ,col='black')
    lines(c(F.N0,F.N0),c(0,tau.star.h),lty=3)
    F.N <- F.N0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-FN.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( F.N.range,tau.hk.F.N,type='l',col='black',ylim=c(bot,top),main='$F_N$')
      lines(F.N.range,tau.h.F.N ,lty=2   ,col='black')
      lines(c(F.N0,F.N0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (kappa.O) other host competence [exogenous]
  {
    kappa.O0       <- kappa.O
    kappa.O.range  <- seq(kappa.O_range[1],kappa.O_range[3],length=100)
    tau.h.kappa.O  <- matrix(0,length(kappa.O.range),1)
    tau.hk.kappa.O <- matrix(0,length(kappa.O.range),1)
    i <- 0
    for(kappa.O in kappa.O.range){i <- i + 1
      cat('\014');cat('Working on kappa.O =',sprintf('%-.2f\n',kappa.O))
      tau.h.kappa.O[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.kappa.O[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.kappa.O,tau.h.kappa.O))*1.05
    bot <- min(c(tau.hk.kappa.O,tau.h.kappa.O))*0.95
    plot( kappa.O.range,tau.hk.kappa.O,type='l',col='black',ylim=c(bot,top))
    lines(kappa.O.range,tau.h.kappa.O ,lty=2   ,col='black')
    lines(c(kappa.O0,kappa.O0),c(0,tau.star.h),lty=3)
    kappa.O <- kappa.O0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-kappaO.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( kappa.O.range,tau.hk.kappa.O,type='l',col='black',ylim=c(bot,top),main='$\\kappa_O$')
      lines(kappa.O.range,tau.h.kappa.O ,lty=2   ,col='black')
      lines(c(kappa.O0,kappa.O0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (gamma) deer hunting catchability coefficient [exogenous]
  {
    gam0       <- gam
    gam.range  <- seq(gam_range[1],gam_range[3],length=100)
    tau.h.gam  <- matrix(0,length(gam.range),1)
    tau.hk.gam <- matrix(0,length(gam.range),1)
    i <- 0
    for(gam in gam.range){i <- i + 1
      cat('\014');cat('Working on gam =',sprintf('%-.2f\n',gam))
      tau.h.gam[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.gam[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.gam,tau.h.gam))*1.05
    bot <- min(c(tau.hk.gam,tau.h.gam))*0.95
    plot( gam.range,tau.hk.gam,type='l',col='black',ylim=c(bot,top))
    lines(gam.range,tau.h.gam ,lty=2   ,col='black')
    lines(c(gam0,gam0),c(0,tau.star.h),lty=3)
    gam <- gam0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-gam.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( gam.range,tau.hk.gam,type='l',col='black',ylim=c(bot,top),main='$\\gamma$')
      lines(gam.range,tau.h.gam ,lty=2   ,col='black')
      lines(c(gam0,gam0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }

  # (M) mice abundance [exogenous]
  {
    M0       <- M
    M.range  <- seq(M_range[1],M_range[3],length=100)
    tau.h.M  <- matrix(0,length(M.range),1)
    tau.hk.M <- matrix(0,length(M.range),1)
    i <- 0
    for(M in M.range){i <- i + 1
      cat('\014');cat('Working on M =',sprintf('%-.1f\n',M))
      tau.h.M[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.M[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.M,tau.h.M))*1.05
    bot <- min(c(tau.hk.M,tau.h.M))*0.95
    plot( M.range,tau.hk.M,type='l',col='black',ylim=c(bot,top))
    lines(M.range,tau.h.M ,lty=2   ,col='black')
    lines(c(O0,O0),c(0,tau.star.h),lty=3)
    M <- M0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-M.tex',sep=''),
        width      = 2.1,
        height     = 2.1,
        pointsize  = 12,
        standAlone = TRUE)
        par(mar=c(2,2,1,1)) # bottom, left, top, right
        plot( M.range,tau.hk.M,type='l',col='black',ylim=c(bot,top),main='$M$')
        lines(M.range,tau.h.M ,lty=2   ,col='black')
        lines(c(M0,M0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BD.N) body burden of nymphs on deer [endogenous]
  {
    BD.N0       <- BD.N
    BD.N.range  <- seq(BD.N_range[1],BD.N_range[3],length=100)
    tau.h.BD.N  <- matrix(0,length(BD.N.range),1)
    tau.hk.BD.N <- matrix(0,length(BD.N.range),1)
    i <- 0
    for(BD.N in BD.N.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BD.N =',sprintf('%-.2f\n',BD.N))
      tau.h.BD.N[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BD.N[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BD.N,tau.h.BD.N))*1.05
    bot <- min(c(tau.hk.BD.N,tau.h.BD.N))*0.95
    plot( BD.N.range,tau.hk.BD.N,type='l',col='black',ylim=c(bot,top))
    lines(BD.N.range,tau.h.BD.N ,lty=2   ,col='black')
    lines(c(BD.N0,BD.N0),c(0,tau.star.h),lty=3)
    BD.N <- BD.N0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BDN.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BD.N.range,tau.hk.BD.N,type='l',col='black',ylim=c(bot,top),main='$B_N^D$')
      lines(BD.N.range,tau.h.BD.N ,lty=2   ,col='black')
      lines(c(BD.N0,BD.N0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }

  # (kappa.M) mice competence [exogenous]
  {
    kappa.M0       <- kappa.M
    kappa.M.range  <- seq(kappa.M_range[1],kappa.M_range[3],length=100)
    tau.h.kappa.M  <- matrix(0,length(kappa.M.range),1)
    tau.hk.kappa.M <- matrix(0,length(kappa.M.range),1)
    i <- 0
    for(kappa.M in kappa.M.range){i <- i + 1
      cat('\014');cat('Working on kappa.M =',sprintf('%-.2f\n',kappa.M))
      tau.h.kappa.M[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.kappa.M[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.kappa.M,tau.h.kappa.M))*1.05
    bot <- min(c(tau.hk.kappa.M,tau.h.kappa.M))*0.95
    plot( kappa.M.range,tau.hk.kappa.M,type='l',col='black',ylim=c(bot,top))
    lines(kappa.M.range,tau.h.kappa.M ,lty=2   ,col='black')
    lines(c(kappa.M0,kappa.M0),c(0,tau.star.h),lty=3)
    kappa.M <- kappa.M0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-kappaM.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( kappa.M.range,tau.hk.kappa.M,type='l',col='black',ylim=c(bot,top),main='$\\kappa_M$')
      lines(kappa.M.range,tau.h.kappa.M ,lty=2   ,col='black')
      lines(c(kappa.M0,kappa.M0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BD.A) body burden of adults on deer [endogenous]
  {
    BD.A0       <- BD.A
    BD.A.range  <- seq(BD.A_range[1],BD.A_range[3],length=100)
    tau.h.BD.A  <- matrix(0,length(BD.A.range),1)
    tau.hk.BD.A <- matrix(0,length(BD.A.range),1)
    i <- 0
    for(BD.A in BD.A.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BD.A =',sprintf('%-.2f\n',BD.A))
      tau.h.BD.A[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BD.A[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BD.A,tau.h.BD.A))*1.05
    bot <- min(c(tau.hk.BD.A,tau.h.BD.A))*0.95
    plot( BD.A.range,tau.hk.BD.A,type='l',col='black',ylim=c(bot,top))
    lines(BD.A.range,tau.h.BD.A ,lty=2   ,col='black')
    lines(c(BD.A0,BD.A0),c(0,tau.star.h),lty=3)
    BD.A <- BD.A0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BDA.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BD.A.range,tau.hk.BD.A,type='l',col='black',ylim=c(bot,top),main='$B_A^D$')
      lines(BD.A.range,tau.h.BD.A ,lty=2   ,col='black')
      lines(c(BD.A0,BD.A0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BM.L) body burden of larvae on mice [endogenous]
  {
    BM.L0       <- BM.L
    BM.L.range  <- seq(BM.L_range[1],BM.L_range[3],length=100)
    tau.h.BM.L  <- matrix(0,length(BM.L.range),1)
    tau.hk.BM.L <- matrix(0,length(BM.L.range),1)
    i <- 0
    for(BM.L in BM.L.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BM.L =',sprintf('%-.2f\n',BM.L))
      tau.h.BM.L[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BM.L[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BM.L,tau.h.BM.L))*1.05
    bot <- min(c(tau.hk.BM.L,tau.h.BM.L))*0.95
    plot( BM.L.range,tau.hk.BM.L,type='l',col='black',ylim=c(bot,top))
    lines(BM.L.range,tau.h.BM.L ,lty=2   ,col='black')
    lines(c(BM.L0,BM.L0),c(0,tau.star.h),lty=3)
    BM.L <- BM.L0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BML.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BM.L.range,tau.hk.BM.L,type='l',col='black',ylim=c(bot,top),main='$B_L^M$')
      lines(BM.L.range,tau.h.BM.L ,lty=2   ,col='black')
      lines(c(BM.L0,BM.L0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (F.L) host finding success of larval ticks [endogenous]
  {
    F.L0       <- F.L
    F.L.range  <- seq(F.L_range[1],F.L_range[3],length=100)
    tau.h.F.L  <- matrix(0,length(F.L.range),1)
    tau.hk.F.L <- matrix(0,length(F.L.range),1)
    i <- 0
    for(F.L in F.L.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on F.L =',sprintf('%-.2f\n',F.L))
      tau.h.F.L[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.F.L[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.F.L,tau.h.F.L))*1.05
    bot <- min(c(tau.hk.F.L,tau.h.F.L))*0.95
    plot( F.L.range,tau.hk.F.L,type='l',col='black',ylim=c(bot,top))
    lines(F.L.range,tau.h.F.L ,lty=2   ,col='black')
    lines(c(F.L0,F.L0),c(0,tau.star.h),lty=3)
    F.L <- F.L0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-FL.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( F.L.range,tau.hk.F.L,type='l',col='black',ylim=c(bot,top),main='$F_L$')
      lines(F.L.range,tau.h.F.L ,lty=2   ,col='black')
      lines(c(F.L0,F.L0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BD.L) body burden of larvae on deer [endogenous]
  {
    BD.L0       <- BD.L
    BD.L.range  <- seq(BD.L_range[1],BD.L_range[3],length=100)
    tau.h.BD.L  <- matrix(0,length(BD.L.range),1)
    tau.hk.BD.L <- matrix(0,length(BD.L.range),1)
    i <- 0
    for(BD.L in BD.L.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BD.L =',sprintf('%-.2f\n',BD.L))
      tau.h.BD.L[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BD.L[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BD.L,tau.h.BD.L))*1.05
    bot <- min(c(tau.hk.BD.L,tau.h.BD.L))*0.95
    plot( BD.L.range,tau.hk.BD.L,type='l',col='black',ylim=c(bot,top))
    lines(BD.L.range,tau.h.BD.L ,lty=2   ,col='black')
    lines(c(BD.L0,BD.L0),c(0,tau.star.h),lty=3)
    BD.L <- BD.L0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BDL.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BD.L.range,tau.hk.BD.L,type='l',col='black',ylim=c(bot,top),main='$B_L^D$')
      lines(BD.L.range,tau.h.BD.L ,lty=2   ,col='black')
      lines(c(BD.L0,BD.L0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (p.h) price per hunting trip [exogenous]
  {
    p.h0       <- p.h
    p.h.range  <- seq(p.h_range[1],p.h_range[3],length=100)
    tau.h.p.h  <- matrix(0,length(p.h.range),1)
    tau.hk.p.h <- matrix(0,length(p.h.range),1)
    i <- 0
    for(p.h in p.h.range){i <- i + 1
      cat('\014');cat('Working on p.h =',sprintf('%-.1f\n',p.h))
      tau.h.p.h[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.p.h[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.p.h,tau.h.p.h))*1.05
    bot <- min(c(tau.hk.p.h,tau.h.p.h))*0.95
    plot( p.h.range,tau.hk.p.h,type='l',col='black',ylim=c(bot,top))
    lines(p.h.range,tau.h.p.h ,lty=2   ,col='black')
    lines(c(p.h0,p.h0),c(0,tau.star.h),lty=3)
    p.h <- p.h0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-ph.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( p.h.range,tau.hk.p.h,type='l',col='black',ylim=c(bot,top),main='$p_h$')
      lines(p.h.range,tau.h.p.h ,lty=2   ,col='black')
      lines(c(p.h0,p.h0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (kappa.D) deer competence [exogenous]
  {
    kappa.D0       <- kappa.D
    kappa.D.range  <- seq(kappa.D_range[1],kappa.D_range[3],length=100)
    tau.h.kappa.D  <- matrix(0,length(kappa.D.range),1)
    tau.hk.kappa.D <- matrix(0,length(kappa.D.range),1)
    i <- 0
    for(kappa.D in kappa.D.range){i <- i + 1
      cat('\014');cat('Working on kappa.D =',sprintf('%-.2f\n',kappa.D))
      tau.h.kappa.D[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.kappa.D[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.kappa.D,tau.h.kappa.D))*1.05
    bot <- min(c(tau.hk.kappa.D,tau.h.kappa.D))*0.95
    plot( kappa.D.range,tau.hk.kappa.D,type='l',col='black',ylim=c(bot,top))
    lines(kappa.D.range,tau.h.kappa.D ,lty=2   ,col='black')
    lines(c(kappa.D0,kappa.D0),c(0,tau.star.h),lty=3)
    kappa.D <- kappa.D0
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-kappaD.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( kappa.D.range,tau.hk.kappa.D,type='l',col='black',ylim=c(bot,top),main='$\\kappa_D$')
      lines(kappa.D.range,tau.h.kappa.D ,lty=2   ,col='black')
      lines(c(kappa.D0,kappa.D0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (x.h) demand for hunting trips [endogenous]
  {
    x.h0       <- x.h
    x.h.range  <- seq(x.h_range[1],x.h_range[3],length=100)
    tau.h.x.h  <- matrix(0,length(x.h.range),1)
    tau.hk.x.h <- matrix(0,length(x.h.range),1)
    i <- 0
    for(x.h in x.h.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on x.h =',sprintf('%-.2f\n',x.h))
      tau.h.x.h[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.x.h[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.x.h,tau.h.x.h))*1.05
    bot <- min(c(tau.hk.x.h,tau.h.x.h))*0.95
    plot( x.h.range,tau.hk.x.h,type='l',col='black',ylim=c(bot,top))
    lines(x.h.range,tau.h.x.h ,lty=2   ,col='black')
    lines(c(x.h0,x.h0),c(0,tau.star.h),lty=3)
    x.h <- x.h0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-xh.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( x.h.range,tau.hk.x.h,type='l',col='black',ylim=c(bot,top),main='$x_h$')
      lines(x.h.range,tau.h.x.h ,lty=2   ,col='black')
      lines(c(x.h0,x.h0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (x.k) demand for hiking trips [endogenous]
  {
    x.k0       <- x.k
    x.k.range  <- seq(x.k_range[1],x.k_range[3],length=100)
    tau.h.x.k  <- matrix(0,length(x.k.range),1)
    tau.hk.x.k <- matrix(0,length(x.k.range),1)
    i <- 0
    for(x.k in x.k.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on x.k =',sprintf('%-.2f\n',x.k))
      tau.h.x.k[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.x.k[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.x.k,tau.h.x.k))*1.05
    bot <- min(c(tau.hk.x.k,tau.h.x.k))*0.95
    plot( x.k.range,tau.hk.x.k,type='l',col='black',ylim=c(bot,top))
    lines(x.k.range,tau.h.x.k ,lty=2   ,col='black')
    lines(c(x.k0,x.k0),c(0,tau.star.h),lty=3)
    x.k <- x.k0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-xk.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( x.k.range,tau.hk.x.k,type='l',col='black',ylim=c(bot,top),main='$x_k$')
      lines(x.k.range,tau.h.x.k ,lty=2   ,col='black')
      lines(c(x.k0,x.k0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BO.L) body burden of larvae on other hosts [endogenous]
  {
    BO.L0       <- BO.L
    BO.L.range  <- seq(BO.L_range[1],BO.L_range[3],length=100)
    tau.h.BO.L  <- matrix(0,length(BO.L.range),1)
    tau.hk.BO.L <- matrix(0,length(BO.L.range),1)
    i <- 0
    for(BO.L in BO.L.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BO.L =',sprintf('%-.2f\n',BO.L))
      tau.h.BO.L[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BO.L[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BO.L,tau.h.BO.L))*1.05
    bot <- min(c(tau.hk.BO.L,tau.h.BO.L))*0.95
    plot( BO.L.range,tau.hk.BO.L,type='l',col='black',ylim=c(bot,top))
    lines(BO.L.range,tau.h.BO.L ,lty=2   ,col='black')
    lines(c(BO.L0,BO.L0),c(0,tau.star.h),lty=3)
    BO.L <- BO.L0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BOL.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BO.L.range,tau.hk.BO.L,type='l',col='black',ylim=c(bot,top),main='$B_L^O$')
      lines(BO.L.range,tau.h.BO.L ,lty=2   ,col='black')
      lines(c(BO.L0,BO.L0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (e.pk) hiking demand price elasticity [endogenous]
  {
    e.pk0       <- e.pk
    e.pk.range  <- seq(e.pk_range[1],e.pk_range[3],length=100)
    tau.h.e.pk  <- matrix(0,length(e.pk.range),1)
    tau.hk.e.pk <- matrix(0,length(e.pk.range),1)
    i <- 0
    for(e.pk in e.pk.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on e.pk =',sprintf('%-.2f\n',e.pk))
      tau.h.e.pk[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.e.pk[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.e.pk,tau.h.e.pk))*1.05
    bot <- min(c(tau.hk.e.pk,tau.h.e.pk))*0.95
    plot( e.pk.range,tau.hk.e.pk,type='l',col='black',ylim=c(bot,top))
    lines(e.pk.range,tau.h.e.pk ,lty=2   ,col='black')
    lines(c(e.pk0,e.pk0),c(0,tau.star.h),lty=3)
    e.pk <- e.pk0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-epk.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( e.pk.range,tau.hk.e.pk,type='l',col='black',ylim=c(bot,top),main='$\\varepsilon_{p_k}$')
      lines(e.pk.range,tau.h.e.pk ,lty=2   ,col='black')
      lines(c(e.pk0,e.pk0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BO.A) body burden of adults on other hosts [endogenous]
  {
    BO.A0       <- BO.A
    BO.A.range  <- seq(BO.A_range[1],BO.A_range[3],length=100)
    tau.h.BO.A  <- matrix(0,length(BO.A.range),1)
    tau.hk.BO.A <- matrix(0,length(BO.A.range),1)
    i <- 0
    for(BO.A in BO.A.range){i <- i + 1
      calibrate.fn()
      cat('\014');cat('Working on BO.A =',sprintf('%-.2f\n',BO.A))
      tau.h.BO.A[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
      tau.hk.BO.A[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BO.A,tau.h.BO.A))*1.05
    bot <- min(c(tau.hk.BO.A,tau.h.BO.A))*0.95
    plot( BO.A.range,tau.hk.BO.A,type='l',col='black',ylim=c(bot,top))
    lines(BO.A.range,tau.h.BO.A ,lty=2   ,col='black')
    lines(c(BO.A0,BO.A0),c(0,tau.star.h),lty=3)
    BO.A <- BO.A0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BOA.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BO.A.range,tau.hk.BO.A,type='l',col='black',ylim=c(bot,top),main='$B_A^O$')
      lines(BO.A.range,tau.h.BO.A ,lty=2   ,col='black')
      lines(c(BO.A0,BO.A0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }
  
  # (BM.A) body burden of adults on mice [endogenous]
  {
    BM.A0       <- BM.A
    BM.A.range  <- seq(BM.A_range[1],BM.A_range[3],length=100)
    tau.h.BM.A  <- matrix(0,length(BM.A.range),1)
    tau.hk.BM.A <- matrix(0,length(BM.A.range),1)
    i <- 0
    for(BM.A in BM.A.range){i <- i + 1
    calibrate.fn()
    cat('\014');cat('Working on BM.A =',sprintf('%-.2f\n',BM.A))
    tau.h.BM.A[i]  <- optimize(CSh.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    tau.hk.BM.A[i] <- optimize(CShk.fn,lower=0,upper=1000,maximum=TRUE,tol=0.01)$maximum
    }
    top <- max(c(tau.hk.BM.A,tau.h.BM.A))*1.05
    bot <- min(c(tau.hk.BM.A,tau.h.BM.A))*0.95
    plot( BM.A.range,tau.hk.BM.A,type='l',col='black',ylim=c(bot,top))
    lines(BM.A.range,tau.h.BM.A ,lty=2   ,col='black')
    lines(c(BM.A0,BM.A0),c(0,tau.star.h),lty=3)
    BM.A <- BM.A0
    calibrate.fn()
    
    if(TeXfigs){
      
      tikz(paste(output.path,'/fig2-BMA.tex',sep=''),
           width      = 2.1,
           height     = 2.1,
           pointsize  = 12,
           standAlone = TRUE)
      par(mar=c(2,2,1,1)) # bottom, left, top, right
      plot( BM.A.range,tau.hk.BM.A,type='l',col='black',ylim=c(bot,top),main='$B_A^M$')
      lines(BM.A.range,tau.h.BM.A ,lty=2   ,col='black')
      lines(c(BM.A0,BM.A0),c(0,tau.star.h),lty=3)
      dev.off()
      
    }
    
  }

}

# SENSITIVITY ANALYSIS -> TORNADO PLOT:
if(TRUE){
  
  tor.tau    <- matrix(0,33,2) # 33 parameters or endogenous calibration variables
  parm.names <- c()
  teX.names  <- c()
  
  # tornado bars
  {
    Dtau.0 <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
    
    # exogenous parameters (except R):
    if(TRUE){
      
      parm.names[1] <- "M" # exogenous
      teX.names[1]  <- "$M$"
      M0 <- M
      M <- min(M_range); tor.tau[1,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      M <- max(M_range); tor.tau[1,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      M <- M0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[2] <- "O" # exogenous
      teX.names[2]  <- "$O$"
      O0 <- O
      O <- min(O_range); tor.tau[2,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      O <- max(O_range); tor.tau[2,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      O <- O0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[3] <- "sigma.A" # exogenous
      teX.names[3]  <- "$\\sigma_A$"
      sigma.A0 <- sigma.A
      sigma.A <- min(sigma.A_range); tor.tau[3,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      sigma.A <- max(sigma.A_range); tor.tau[3,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      sigma.A <- sigma.A0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[4] <- "kappa.M" # exogenous
      teX.names[4]  <- "$\\kappa_M$"
      kappa.M0 <- kappa.M
      kappa.M <- min(kappa.M_range); tor.tau[4,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.M <- max(kappa.M_range); tor.tau[4,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.M <- kappa.M0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[5] <- "kappa.O" # exogenous
      teX.names[5] <- "$\\kappa_O$"
      kappa.O0 <- kappa.O
      kappa.O <- min(kappa.O_range); tor.tau[5,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.O <- max(kappa.O_range); tor.tau[5,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.O <- kappa.O0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[6] <- "kappa.D" # exogenous
      teX.names[6] <- "$\\kappa_D$"
      kappa.D0 <- kappa.D
      kappa.D <- min(kappa.D_range); tor.tau[6,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.D <- max(kappa.D_range); tor.tau[6,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      kappa.D <- kappa.D0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}

      parm.names[7] <- "r" # exogenous
      teX.names[7]  <- "$r$"
      r0 <- r
      r <- min(r_range); tor.tau[7,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      r <- max(r_range); tor.tau[7,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      r <- r0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[8] <- "p.h" # exogenous
      teX.names[8]  <- "$p_h$"
      p.h0 <- p.h
      p.h <- min(p.h_range); tor.tau[8,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      p.h <- max(p.h_range); tor.tau[8,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      p.h <- p.h0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[9] <- "gamma" # exogenous
      teX.names[9]  <- "$\\gamma$"
      gam0 <- gam
      gam <- min(gam_range); tor.tau[9,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      gam <- max(gam_range); tor.tau[9,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      gam <- gam0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}

      parm.names[10] <- "n.k" # exogenous
      teX.names[10] <- "$n_k$"
      n.k0 <- n.k
      n.k <- min(n.k_range); tor.tau[10,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      n.k <- max(n.k_range); tor.tau[10,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      n.k <- n.k0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[11] <- "p.k" # exogenous
      teX.names[11]  <- "$p_k$"
      p.k0 <- p.k
      p.k <- min(p.k_range); tor.tau[11,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      p.k <- max(p.k_range); tor.tau[11,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      p.k <- p.k0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[12] <- "Z" # exogenous
      teX.names[12]  <- "$Z$"
      Z0 <- Z
      Z <- min(Z_range); tor.tau[12,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      Z <- max(Z_range); tor.tau[12,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      Z <- Z0
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}

    }
    
    # endogenous calibration variables:
    if(TRUE){
      
      parm.names[13] <- "R" 
      teX.names[13]  <- "($R$)"
      R0 <- R
      R <- min(R_range); calibrate.fn(); tor.tau[13,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      R <- max(R_range); calibrate.fn(); tor.tau[13,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      R <- R0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[14] <- "D" # endogenous calibration variable
      teX.names[14]  <- "($D$)"
      D0 <- D
      D <- min(D_range); calibrate.fn(); tor.tau[14,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      D <- max(D_range); calibrate.fn(); tor.tau[14,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      D <- D0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[15] <- "F.L" # endogenous (depends on D)
      teX.names[15]  <- "($F_L$)"
      F.L0 <- F.L
      F.L <- min(F.L_range); calibrate.fn(); tor.tau[15,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.L <- max(F.L_range); calibrate.fn(); tor.tau[15,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.L <- F.L0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[16] <- "F.N" # endogenous (depends on D)
      teX.names[16]  <- "($F_N$)"
      F.N0 <- F.N
      F.N <- min(F.N_range); calibrate.fn(); tor.tau[16,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.N <- max(F.N_range); calibrate.fn(); tor.tau[16,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.N <- F.N0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[17] <- "F.A" # endogenous (depends on D)
      teX.names[17]  <- "($F_A$)"
      F.A0 <- F.A
      F.A <- min(F.A_range); calibrate.fn(); tor.tau[17,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.A <- max(F.A_range); calibrate.fn(); tor.tau[17,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      F.A <- F.A0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[18] <- "BM.L" # endogenous
      teX.names[18]  <- "($B_L^M$)"
      BM.L0 <- BM.L
      BM.L <- min(BM.L_range); calibrate.fn(); tor.tau[18,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.L <- max(BM.L_range); calibrate.fn(); tor.tau[18,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.L <- BM.L0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[19] <- "BM.N" # endogenous
      teX.names[19]  <- "($B_N^M$)"
      BM.N0 <- BM.N
      BM.N <- min(BM.N_range); calibrate.fn(); tor.tau[19,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.N <- max(BM.N_range); calibrate.fn(); tor.tau[19,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.N <- BM.N0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[20] <- "BM.A" # endogenous
      teX.names[20]  <- "($B_A^M$)"
      BM.A0 <- BM.A
      BM.A <- min(BM.A_range); calibrate.fn(); tor.tau[20,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.A <- max(BM.A_range); calibrate.fn(); tor.tau[20,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BM.A <- BM.A0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[21] <- "BO.L" # endogenous
      teX.names[21]  <- "($B_L^O$)"
      BO.L0 <- BO.L
      BO.L <- min(BO.L_range); calibrate.fn(); tor.tau[21,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.L <- max(BO.L_range); calibrate.fn(); tor.tau[21,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.L <- BO.L0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[22] <- "BO.N" # endogenous
      teX.names[22]  <- "($B_N^O$)"
      BO.N0 <- BO.N
      BO.N <- min(BO.N_range); calibrate.fn(); tor.tau[22,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.N <- max(BO.N_range); calibrate.fn(); tor.tau[22,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.N <- BO.N0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[23] <- "BO.A" # endogenous
      teX.names[23]  <- "($B_A^O$)"
      BO.A0 <- BO.A
      BO.A <- min(BO.A_range); calibrate.fn(); tor.tau[23,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.A <- max(BO.A_range); calibrate.fn(); tor.tau[23,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BO.A <- BO.A0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[24] <- "BD.L" # endogenous
      teX.names[24]  <- "($B_L^D$)"
      BD.L0 <- BD.L
      BD.L <- min(BD.L_range); calibrate.fn(); tor.tau[24,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.L <- max(BD.L_range); calibrate.fn(); tor.tau[24,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.L <- BD.L0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[25] <- "BD.N" # endogenous
      teX.names[25]  <- "($B_N^D$)"
      BD.N0 <- BD.N
      BD.N <- min(BD.N_range); calibrate.fn(); tor.tau[25,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.N <- max(BD.N_range); calibrate.fn(); tor.tau[25,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.N <- BD.N0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[26] <- "BD.A" # endogenous
      teX.names[26]  <- "($B_A^D$)"
      BD.A0 <- BD.A
      BD.A <- min(BD.A_range); calibrate.fn(); tor.tau[26,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.A <- max(BD.A_range); calibrate.fn(); tor.tau[26,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      BD.A <- BD.A0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[27] <- "e.ph" # endogenous
      teX.names[27]  <- "($\\varepsilon_{p_h}$)"
      e.ph0 <- e.ph
      e.ph <- min(e.ph_range); calibrate.fn(); tor.tau[27,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.ph <- max(e.ph_range); calibrate.fn(); tor.tau[27,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.ph <- e.ph0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[28] <- "e.gamD" # endogenous
      teX.names[28]  <- "($\\varepsilon_{\\gamma D}$)"
      e.gamD0 <- e.gamD
      e.gamD  <- min(e.gamD_range); calibrate.fn(); tor.tau[28,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.gamD  <- max(e.gamD_range); calibrate.fn(); tor.tau[28,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.gamD  <- e.gamD0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}

      parm.names[29] <- "e.pk" # endogenous
      teX.names[29]  <- "($\\varepsilon_{p_k}$)"
      e.pk0 <- e.pk
      e.pk <- min(e.pk_range); calibrate.fn(); tor.tau[29,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.pk <- max(e.pk_range); calibrate.fn(); tor.tau[29,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      e.pk <- e.pk0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[30] <- "ell" # endogenous
      teX.names[30]  <- "($\\ell$)"
      ell0 <- ell
      ell <- min(ell_range); calibrate.fn(); tor.tau[30,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      ell <- max(ell_range); calibrate.fn(); tor.tau[30,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      ell <- ell0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[31] <- "mr.ratio" # endogenous
      teX.names[31]  <- "($m/r$)"
      mr.ratio0 <- mr.ratio
      mr.ratio <- min(mr.ratio_range); calibrate.fn(); tor.tau[31,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      mr.ratio <- max(mr.ratio_range); calibrate.fn(); tor.tau[31,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      mr.ratio <- mr.ratio0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[32] <- "x.h" # endogenous
      teX.names[32]  <- "($x_h$)"
      x.h0 <- x.h
      x.h <- min(x.h_range); calibrate.fn(); tor.tau[32,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      x.h <- max(x.h_range); calibrate.fn(); tor.tau[32,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      x.h <- x.h0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}
      
      parm.names[33] <- "x.k" # endogenous
      teX.names[33]  <- "($x_k$)"
      x.k0 <- x.k
      x.k <- min(x.k_range); calibrate.fn(); tor.tau[33,1] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      x.k <- max(x.k_range); calibrate.fn(); tor.tau[33,2] <- (tau.hk.opt()-tau.h.opt())/tau.h.opt()
      x.k <- x.k0
      calibrate.fn()
      if(((tau.hk.opt()-tau.h.opt())/tau.h.opt()) != Dtau.0){pause()}

    }
    
    order      <- sort(abs(tor.tau[,2]-tor.tau[,1]),decreasing=FALSE,index.return=TRUE)$ix
    tor.tau    <- tor.tau[order,]
    tor.names  <- parm.names[order]
    torX.names <- teX.names[order]
  }
  
  # plot
  {
    tikz(paste(output.path,'/tornado.tex',sep=''),
         width      = 4,
         height     = 5,
         pointsize  = 12,
         standAlone = TRUE)
    
    par(mfrow=c(1,1))
    par(mar=c(3,0.5,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
    plot(c(tor.tau[1,1],tor.tau[1,2]),c(33,33),type='l',col='white',
         ylim=c(1,33),
         xlim=c(-.65,0),
         xaxt="n",xlab="tornado widths",ylab="",axes=FALSE)
    box(lwd=1)
    for (k in seq(33, 1, -1)) {
      if (tor.tau[k, 1] <= Dtau.0) {
        # Fill the left side of the bar with a different color
        rect(tor.tau[k, 1], k - 0.4, Dtau.0   , k + 0.4, col = '#999999', border = NA)
        rect(Dtau.0   , k - 0.4, tor.tau[k, 2], k + 0.4, col = '#CCCCCC', border = NA)
      }else{
        rect(tor.tau[k, 2], k - 0.4, Dtau.0, k + 0.4, col = '#CCCCCC', border = NA)
        rect(Dtau.0   , k - 0.4, tor.tau[k, 1], k + 0.4, col = '#999999', border = NA)
      }
      text(min(tor.tau[k,]),k,torX.names[k],cex=.65,pos=2)
      # mtext(torX.names[k], side = 2, line = 2, at = k, las = 2, cex=.65)
    }
    mtext(c(-.6,-.5,-.4,-.3,-.2,-.1,0), side = 1, at = c(-.6,-.5,-.4,-.3,-.2,-.1,0), las = 1, cex=.8)
    mtext('$(\\tau_{hk}-\\tau_h)/\\tau_h$', side = 1, at = -.25 , las = 1 , cex = .8, line = 1.2 )
    lines(c(Dtau.0,Dtau.0),c(1,33),lty=3)
    
    dev.off()
  }
  
}

# pause()

# UNCERTAINTY ANALYSIS -> MONTE CARLO SIMULATION:
if(TRUE){

  set.seed(1234)  
  MC <- 5000
  tau.h.MC   <- matrix(0,MC,1)
  tau.hk.MC  <- matrix(0,MC,1)
  II.MC      <- matrix(0,MC,1)
  B.h.MC     <- matrix(0,MC,1)
  B.hk.MC    <- matrix(0,MC,1)
  B.star0.MC <- matrix(0,MC,1)
  B.star.MC  <- matrix(0,MC,1)
  dIdD.MC    <- matrix(0,MC,1)
  CSh.MC     <- matrix(0,MC,1)
  CShk.MC    <- matrix(0,MC,1)
  CSB0.MC    <- matrix(0,MC,1)
  CSB.MC     <- matrix(0,MC,1)
  flag.MC    <- matrix(0,MC,1)
  parms.MC   <- list()
  for(mc in 1:MC){

    # Random parameter draw:
    parms          <- drawparms.fn('uniform')
    parms.MC[[mc]] <- parms
    calibrate.fn()
    
    # dI/dD at prevailing D:
    DD <- c(D,D*1.05)
    II <- IvsD.fn(DD)
    dIdD.MC[mc] <- (II[2]-II[1])/(DD[2]-DD[1])
    
     # Optimal hunting fees:
    tau.h.MC[mc]  <- tau.h.opt()
    tau.hk.MC[mc] <- tau.hk.opt()
    CSh.MC[mc]    <- CS.fn(tau.h.MC[mc])$CSh
    CShk.MC[mc]   <- CS.fn(tau.hk.MC[mc])$CShk
    II.MC[mc]     <- eq.fn(tau.hk.MC[mc])$II

    outs <- eq.fn(tau.h.MC[mc])
    flag.MC[mc] <- outs$flag

    # Running histograms:
    if(FALSE){
      par(mfrow=c(1,3))
      hist(dIdD.MC[1:mc],main='dI/dD')
      hist(tau.h.MC[1:mc],main='tau.h')
      hist((tau.hk.MC[1:mc]-tau.h.MC[1:mc])/tau.h.MC[1:mc],main='(tau.hk-tau.h)/tau.h')

    }
    if(TRUE){
      if(floor(mc/1)==(mc/1)){
        
        cat('\014',sprintf('Completed %-.0f of %-0.f Monte Carlo reps\n',mc,MC))
        
        par(mfrow = c(1,2),        # 2x1 layout
            oma   = c(0, 0, 0, 0), # rows of text at the outer bottom, left, top, right
            mar   = c(4, 3, 1, 1), # rows of text at ticks and to separate plots
            mgp   = c(2, 1, 0) )   # axis label at 2 rows distance, tick labels at 1 row
        
        hist1 <- hist(tau.h.MC[1:mc],breaks=seq(min(tau.h.MC[1:mc]),max(tau.h.MC[1:mc])+20,20),plot=FALSE)
        plot(hist1$mids,hist1$counts/mc,type='l',xlim=c(-500,1000),ylab="",yaxt="n",xlab="tauh")
        lines(c(0,0),c(0,1))
        lines(c(tau.star.h,tau.star.h),c(0,1),col='red')

        hist2 <- hist(tau.hk.MC[1:mc],breaks=seq(min(tau.hk.MC[1:mc]),max(tau.hk.MC[1:mc])+20,20),plot=FALSE)
        plot(hist2$mids,hist2$counts/mc,type='l',xlim=c(-500,1000),ylab="",yaxt="n",xlab="tau_hk")
        lines(c(0,0),c(0,1))
        lines(c(tau.star.hk,tau.star.hk),c(0,1),col='red')

        pause(.05)
        
      }
    }
    
    # TeX figure w histograms:
    if(mc == MC){
      
      # tau.h frequency distribution:
      {
        tikz(paste(output.path,'/fig-tauh-mc.tex',sep=''),
             width      = 6,
             height     = 3,
             pointsize  = 12,
             standAlone = TRUE)
        
        par(mfrow=c(1,1),
            mar=c(2,.25,.50,.25), # [bottom,left,top,right]
            oma=c(0,1,1,1))
        
        hist1 <- hist(tau.h.MC[1:mc],breaks=seq(min(tau.h.MC[1:mc]),max(tau.h.MC[1:mc])+20,20),plot=FALSE)
        plot(hist1$mids,hist1$counts/mc,type='l',xlim=c(-500,1000),ylab="",yaxt="n",xlab="")
        lines(c(0,0),c(0,1))
        lines(c(tau.star.h,tau.star.h),c(0,1),col='red')
        mtext("Hunter surplus-maximizing fee",side=3,line=.5,cex=1)
        mtext(sprintf("{\\small %-.1f}",tau.star.h),side=1,at=tau.star.h,line=.3,cex=1)
        
        dev.off()
      }
      
      # tau.hk frequency distribution:
      {
        tikz(paste(output.path,'/fig-tauhk-mc.tex',sep=''),
             width      = 6,
             height     = 3,
             pointsize  = 12,
             standAlone = TRUE)
        
        par(mfrow=c(1,1),
            mar=c(2,.25,.50,.25), # [bottom,left,top,right]
            oma=c(0,1,1,1))
        
        hist1 <- hist(tau.hk.MC[1:mc],breaks=seq(min(tau.hk.MC[1:mc]),max(tau.hk.MC[1:mc])+20,20),plot=FALSE)
        plot(hist1$mids,hist1$counts/mc,type='l',xlim=c(-500,1000),ylab="",yaxt="n",xlab="")
        lines(c(0,0),c(0,1))
        lines(c(tau.star.hk,tau.star.hk),c(0,1),col='red')
        mtext("Hunter+hiker surplus-maximizing fee",side=3,line=.5,cex=1)
        mtext(sprintf("{\\small %-.1f}",tau.star.hk),side=1,at=tau.star.hk,line=.3,cex=1)
        
        dev.off()
      }
      
      tikz(paste(output.path,'/MC-figs.tex',sep=''),
           width      = 7,
           height     = 4,
           pointsize  = 12,
           standAlone = TRUE)
      
      par(mfrow=c(2,3),
          mar=c(4,4,1,0.5), # [bottom,left,top,right]
          oma=c(0,1,1,1))

      hist(dIdD.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$\\partial I/\\partial D$',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      hist(tau.h.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$\\tau_h$',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      hist((tau.hk.MC-tau.h.MC)/tau.h.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$(\\tau_{hk}-\\tau_h)/\\tau_h$',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      hist(B.star0.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$B^{\\star}(\\tau=0)$',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      hist(B.star.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$B^{\\star}(\\tau=\\tau^{\\star})$',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      hist(CSB0.MC/CShk.MC,main='',xlab='',ylab='',ylim=c(0,MC))
      mtext('$CS$ ratio',side=3,line=.5,cex=.7)
      box(lty=1, col="black")
      
      dev.off()
      
    }
    
  }
  
  if(TeXfigs==TRUE){
    # fig a: dI/dD histogram
    {  
      tikz(paste(output.path,'/fig-mc-a.tex',sep=''),
         width      = 4,
         height     = 3,
         pointsize  = 12,
         standAlone = TRUE)
      par(mar=c(3.5,3.5,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
      hist(dIdD.MC,
           main     = '',
           xlab     = '',
           ylab     = '',
           axes     = FALSE,
           family   = 'Bookman',
           cex      = 1.25,
           cex.axis = 1.25,
           cex.lab  = 1.25)
      axis(1,las=1) # Draw x axis
      axis(2,las=0) # Draw y axis
      mtext('Frequency',side=2, line=2.5, cex.lab=1.2, las=0) # y-axis label
      mtext('$I_D$' ,side=1, line=2.5, cex.lab=1.2, las=1) # x-axis
      dev.off()
    }
    
    # fig b: tau.h histogram
    {  
      tikz(paste(output.path,'/fig-mc-b.tex',sep=''),
          width      = 4,
          height     = 3,
          pointsize  = 12,
          standAlone = TRUE)
      par(mar=c(3.5,3.5,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
      hist(tau.h.MC,
           main     = '',
           xlab     = '',
           ylab     = '',
           axes     = FALSE,
           family   = 'Bookman',
           cex      = 1.25,
           cex.axis = 1.25,
           cex.lab  = 1.25)
      axis(1,las=1) # Draw x axis
      axis(2,las=0) # Draw y axis
      mtext('Frequency',side=2, line=2.5, cex.lab=1.2, las=0) # y-axis label
      mtext('$\\tau_h$' ,side=1, line=2.5, cex.lab=1.2, las=1) # x-axis
      dev.off()
    }
    
    # fig c: tau.hk histogram
    {  
      tikz(paste(output.path,'/fig-mc-c.tex',sep=''),
         width      = 4,
         height     = 3,
         pointsize  = 12,
         standAlone = TRUE)
      par(mar=c(3.5,3.5,0.5,0.5)) # figure border whitespace [bottom,left,top,right]
      hist(tau.hk.MC,
           main     = '',
           xlab     = '',
           ylab     = '',
           axes     = FALSE,
           family   = 'Bookman',
           cex      = 1.25,
           cex.axis = 1.25,
           cex.lab  = 1.25)
      axis(1,las=1) # Draw x axis
      axis(2,las=0) # Draw y axis
      mtext('Frequency',side=2, line=2.5, cex.lab=1.2, las=0) # y-axis label
      mtext('$\\tau_{hk}-\\tau_h$' ,side=1, line=2.5, cex.lab=1.2, las=1) # x-axis
      dev.off()
    }

  }
  
  # Print summary stats to output file:
  cat('\n\nMonte Carlo summary stats:\n',file=out.file.name,append=TRUE)
  cat(sprintf('St dev of tau.h outcomes  = %-.4f\n',sd(tau.h.MC)),file=out.file.name,append=TRUE)
  cat(sprintf('St dev of tau.hk outcomes = %-.4f\n',sd(tau.hk.MC)),file=out.file.name,append=TRUE)
  cat(sprintf('Percent of tau.h outcomes larger than benchmark tau.h = %-.4f\n',mean(tau.h.MC>tau.star.h)),file=out.file.name,append=TRUE)
  cat(sprintf('Percent of tau.hk outcomes larger than benchmark tau.hk = %-.4f\n',mean(tau.hk.MC>tau.star.hk)),file=out.file.name,append=TRUE)
  cat(sprintf('Percent of tau.hk outcomes less than zero = %-.4f\n',mean(tau.hk.MC<0)),file=out.file.name,append=TRUE)
  cat(sprintf('Percent of outcomes with tau.hk less than tau.h = %-.4f\n',mean(tau.hk.MC<tau.h.MC)),file=out.file.name,append=TRUE)
  
}

