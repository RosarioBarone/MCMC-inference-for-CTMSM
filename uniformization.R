
##---------------------------------------------------------------
Uniformization=function(SaT,RateM,abs,Lam,U,InvU,death)
{
  id=SaT[7]
  BgnSt=SaT[1];EndSt=SaT[2];Tm=SaT[4]-SaT[3]
  Fstate=integer()
  Ftime=double()
  FEndSt=EndSt
  FSaT=SaT[4]
  if (EndSt==abs & death==1){
  prob.abs=(round(( U%*%diag(exp(Tm*Lam))%*%InvU)[BgnSt,]*RateM[,abs],digits=10))
  EndSt=sample(nrow(RateM),size=1,prob=prob.abs)
  Fstate=abs
  Ftime=SaT[4]
  FEndSt=integer()
  FSaT=double()}

  ## Diagonalization of rate matrix
  PrbBgnEnd <- (U%*%diag(exp(Tm*Lam))%*%InvU)[BgnSt,EndSt]

  ## Determine max diagonal entry and construct transition matrix
  nSt <- nrow(RateM)
  Mx <- max(-diag(RateM))
  TransM <- diag(rep(1,times=nSt))+RateM/Mx
  TransMn <- TransM
  ## TransMn is n'th power of TransM
  ##------------------------------------------------------------
  ## Simulate number of jumps

  rU <- runif(1)
  ifelse(BgnSt==EndSt,cum <- dpois(0,Mx*Tm)/PrbBgnEnd,cum <- 0)
  notExceed <- TRUE
  if (cum>rU) { notExceed <- FALSE }
  nJmp <- 0
  while (notExceed){
    nJmp <- nJmp+1
    prb <- dpois(nJmp,Mx*Tm)*TransMn[BgnSt,EndSt]/PrbBgnEnd
    cum <- cum+prb
    if (cum>rU) notExceed <- FALSE
    ## Update transition matrices
    ## TransArr holds n'th power of TransM
    if (nJmp==1) TransArr <- array(TransM,c(nSt,nSt,1))
    if (nJmp!=1) TransArr <- array(c(TransArr,TransMn),c(nSt,nSt,nJmp))
    TransMn <- TransMn%*%TransM
  }


#  cat("nJmp:",nJmp,"\n")
  ##--------------------------------------------------------

  ## if (nJmp==0)
  if (nJmp==0){
    Path <- list()
    if (SaT[5]==0 & SaT[6]==0){
    Path$St <- integer()
    Path$Tm <- double()}
    if (SaT[5]==0 & SaT[6]==1){
     Path$St <- c(FEndSt,Fstate)
    Path$Tm <- c(FSaT,Ftime)}
  if (SaT[5]==1 & SaT[6]==0)
    {Path$St=c(BgnSt)
    Path$Tm=c(0)}
 if (SaT[5]==1 & SaT[6]==1)
  {Path$St=c(BgnSt,FEndSt,Fstate)
    Path$Tm=c(0,FSaT,Ftime)}
    Path$id=rep(id,length(Path$Tm))
    return(rbind(Path$St,Path$Tm,Path$id))
  }



  ## if (nJmp==1)
  if (nJmp==1){
    ## if virtual jump: done
    if (BgnSt==EndSt){
      Path <- list()
       if (SaT[5]==0 & SaT[6]==0){
      Path$St <- c(integer(),Fstate)
      Path$Tm <- c(double(),Ftime)
      }
       if (SaT[5]==0 & SaT[6]==1){
      Path$St <- c(FEndSt,Fstate)
      Path$Tm <- c(FSaT,Ftime)}

    if (SaT[5]==1 & SaT[6]==0){
      Path$St <- c(BgnSt,Fstate)
      Path$Tm <- c(0,Ftime)}



      if (SaT[5]==1 & SaT[6]==1){
      Path$St <- c(BgnSt,FEndSt,Fstate)
      Path$Tm <- c(0,FSaT,Ftime)}


      Path$id=rep(id,length(Path$Tm))
    return(rbind(Path$St,Path$Tm,Path$id))


    }


    if (BgnSt!=EndSt){
      Path <- list()

if (SaT[5]==0 & SaT[6]==0){
      Path$St <- c(FEndSt,Fstate)
      Path$Tm <- c(SaT[3]+Tm*runif(1),Ftime)}

if (SaT[5]==0 & SaT[6]==1){
      Path$St <- c(EndSt,FEndSt,Fstate)
      Path$Tm <- c(SaT[3]+Tm*runif(1),FSaT,Ftime)}

 if (SaT[5]==1 & SaT[6]==0){
      Path$St <- c(BgnSt,FEndSt,Fstate)
      Path$Tm <- c(0,SaT[3]+Tm*runif(1),Ftime)}


  if (SaT[5]==1 & SaT[6]==1){
      Path$St <- c(BgnSt,EndSt,FEndSt,Fstate)
      Path$Tm <- c(0,SaT[3]+Tm*runif(1),FSaT,Ftime)}

Path$id=rep(id,length(Path$Tm))
    return(rbind(Path$St,Path$Tm,Path$id))

    }



  }


  ## Case (nJmp >= 2):
  ## Simulate jumping times
  JmpTmV <- Tm*sort(runif(nJmp))
  ## Simulate states (last state always EndSt)
  JmpStV <- rep(0,(nJmp-1))
  Prb1 <- TransM[BgnSt,]
  for (jump in 1:(nJmp-1)){
    Prb2Mat <- TransArr[,,(nJmp-jump)]
    Prb2 <- Prb2Mat[,EndSt]
    JmpStV[jump] <-
      sample(1:nSt,size=1,replace=TRUE,round(Prb1*Prb2,digits=20))
    Prb1 <- TransM[JmpStV[jump],]
  }

  JmpStV <- c(JmpStV,EndSt)
  ## Remove virtual substitutions
  TrueSub <- c(BgnSt,JmpStV[1:(nJmp-1)])!=c(JmpStV[1:nJmp])
  Path <- list()
   if (SaT[5]==0 & SaT[6]==0){
  Path$St <- c(JmpStV[TrueSub],Fstate)
  Path$Tm <- c(SaT[3]+JmpTmV[TrueSub],Ftime)}

if (SaT[5]==0 & SaT[6]==1){
  Path$St <- c(JmpStV[TrueSub],FEndSt,Fstate)
  Path$Tm <- c(SaT[3]+JmpTmV[TrueSub],FSaT,Ftime)}

if (SaT[5]==1 & SaT[6]==0){
  Path$St <- c(BgnSt,JmpStV[TrueSub],Fstate)
  Path$Tm <- c(0,SaT[3]+JmpTmV[TrueSub],Ftime)}

if (SaT[5]==1 & SaT[6]==1){
  Path$St <- c(BgnSt,JmpStV[TrueSub],FEndSt,Fstate)
  Path$Tm <- c(0,SaT[3]+JmpTmV[TrueSub],FSaT,Ftime)}

Path$id=rep(id,length(Path$Tm))
    return(rbind(Path$St,Path$Tm,Path$id))

}




