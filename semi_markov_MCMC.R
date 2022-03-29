gibbs.msm.SM=function(data,NMCMC,RateM.init,alpha.dir,f,g,abs,death,SM){

alpha.dir=matrix(as.numeric(RateM.init>0),ncol=ncol(RateM.init))


Obs = data[,1];States = data[,2];Time = data[,3]

Nstates=length(unique(States))
Nobs=n=length(unique(data$id))
na.states=nrow(RateM.init)-(abs>0)# number of non absorbing state

Obs=as.factor(Obs); levels(Obs)=1:Nobs
Obs=as.numeric(Obs)

dataset = cbind(Obs,States,Time)

ObsStPanel=list();ObsTmPanel= list()
for (i in 1:Nobs){
ObsStPanel[[i]] = dataset[Obs==i,2]
ObsTmPanel[[i]] = dataset[Obs==i,3]
}
ni=sapply(ObsStPanel,FUN=length)

A=list()
for(i in 1:Nobs){
A[[i]]=cbind(
ObsStPanel[[i]][-ni[i]],ObsStPanel[[i]][-1],
ObsTmPanel[[i]][-ni[i]],ObsTmPanel[[i]][-1],
as.numeric((1:(ni[i]-1)==1)),
as.numeric((1:(ni[i]-1)==(ni[i]-1))),
rep(i,ni[i]-1))
}


D=matrix(ncol=7,nrow=sum((ni-1)))
l=1
for(i in 1:Nobs){
for (j in 1:(ni[i]-1)) {D[l,]=A[[i]][j,]
l=l+1}
}

RateM=RateM.init
gamma=-diag(RateM)
alpha=rep(1,length(gamma))

Eigen=eigen(RateM);Lam <- Eigen$values
U <- Eigen$vectors;InvU <- solve(U)


ZtL=parApply(cl,D,FUN=Uniformization,MAR=1,RateM=RateM, Lam=Lam,U=U,InvU=InvU,abs=abs,death=death)



P=matrix(0,nrow=nrow(RateM),ncol=ncol(RateM))

param=matrix(ncol=(ncol(RateM)*(ncol(RateM)+1)+1),nrow=NMCMC)


for(iter in 1:NMCMC){

############### GENERAZIONE TRAIETTORIE ##########################

###############ZtL è l'attuale lista con la sequenza di traiettorie ma ricalcolare alcuni elementi di Zt per cola di alpha
if(SM==1){
Zt=matrix(unlist(ZtL),ncol=3,byrow=TRUE)# prima collona =stati, seconda colonna=tempi,terza colonna=obs
a=unlist(tapply(Zt[,2],INDEX=Zt[,3],FUN=diff))#soggiorni
b=unlist(tapply(Zt[,1],INDEX=Zt[,3],FUN=tail,n=-1))#stato finale
c=unlist(tapply(Zt[,1],INDEX=Zt[,3],FUN=head,n=-1))#stato iniziale
d=a^alpha[c];
e=(c!=b) #non censure
gg=unlist(tapply(Zt[,3],INDEX=Zt[,3],FUN=head,n=-1))#ossevazioni
Zt=cbind(c,b,a,d,e,gg)

TRAN.seq=array(table(Zt[,1],Zt[,2],Zt[,6]),dim=c(3,3,Nobs))
diag.seq=apply(TRAN.seq,FUN=diag,MAR=3)
n.seq=t(apply(TRAN.seq,FUN=sum,MAR=c(1,3))-diag.seq)



SOJ.seq      =t(tapply(Zt[,3],INDEX=list(Zt[,1],Zt[,6]),FUN=sum))
SOJ.seq[is.na(SOJ.seq)]=0

SOJ.alpha.seq=t(tapply(Zt[,4],INDEX=list(Zt[,1],Zt[,6]),FUN=sum))
SOJ.alpha.seq[is.na(SOJ.alpha.seq)]=0

LSOJ.nocen.seq=t(tapply(log(Zt[,3])*Zt[,5],INDEX=list(Zt[,1],Zt[,6]),FUN=sum))
LSOJ.nocen.seq[is.na(LSOJ.nocen.seq)]=0
}
##################propongo le nuove traiettorie#################à

Eigen=eigen(RateM);Lam <- Eigen$values
U <- Eigen$vectors;InvU <- solve(U)

ZtL.prop=parApply(cl,D,FUN=Uniformization,MAR=1,RateM=RateM, Lam=Lam,U=U,InvU=InvU,abs=abs,death=death)
if(SM==0){ZtL=ZtL.prop;acc=Nobs}
else{
Zt.prop=matrix(unlist(ZtL.prop),ncol=3,byrow=TRUE)



a=unlist(tapply(Zt.prop[,2],INDEX=Zt.prop[,3],FUN=diff))
b=unlist(tapply(Zt.prop[,1],INDEX=Zt.prop[,3],FUN=tail,n=-1))
c=unlist(tapply(Zt.prop[,1],INDEX=Zt.prop[,3],FUN=head,n=-1))
d=a^alpha[c]
e=(c!=b)
gg=unlist(tapply(Zt.prop[,3],INDEX=Zt.prop[,3],FUN=head,n=-1))
Zt.prop=cbind(c,b,a,d,e,gg)



TRAN.seq.prop=array(table(Zt.prop[,1],Zt.prop[,2],Zt.prop[,6]),dim=c(3,3,Nobs))
diag.seq.prop=apply(TRAN.seq.prop,FUN=diag,MAR=3)
n.seq.prop=t(apply(TRAN.seq.prop,FUN=sum,MAR=c(1,3))-diag.seq.prop)



SOJ.seq.prop      =t(tapply(Zt.prop[,3],INDEX=list(Zt.prop[,1],Zt.prop[,6]),FUN=sum))
SOJ.seq.prop[is.na(SOJ.seq.prop)]=0

SOJ.alpha.seq.prop=t(tapply(Zt.prop[,4],INDEX=list(Zt.prop[,1],Zt.prop[,6]),FUN=sum))
SOJ.alpha.seq.prop[is.na(SOJ.alpha.seq.prop)]=0

LSOJ.nocen.seq.prop=t(tapply(log(Zt.prop[,3])*Zt.prop[,5],INDEX=list(Zt.prop[,1],Zt.prop[,6]),FUN=sum))
LSOJ.nocen.seq.prop[is.na(LSOJ.nocen.seq.prop)]=0



lnum.seq=-(gamma^alpha)*t(SOJ.alpha.seq.prop)+(log(alpha)+alpha*log(gamma))*t(n.seq.prop)+(alpha-1)*t(LSOJ.nocen.seq.prop)+log(gamma)*t(n.seq)     -gamma*t(SOJ.seq)
lden.seq=-(gamma^alpha)*t(SOJ.alpha.seq)     +(log(alpha)+alpha*log(gamma))*t(n.seq)     +(alpha-1)*t(LSOJ.nocen.seq)     +log(gamma)*t(n.seq.prop)-gamma*t(SOJ.seq.prop)

lnum.seq=colSums(lnum.seq[-abs,])
lden.seq=colSums(lden.seq[-abs,])


acc=runif(Nobs)<exp(lnum.seq-lden.seq) #per stati non assorbenti NON VA BENE


for (i in 1:Nobs) {if(acc[i]) {
indici=which(D[,7]==i)
for(l in indici) {ZtL[[l]]=ZtL.prop[[l]]}}}
}

#####################GENERAZIONE PARAMETRI##########################à
Zt=matrix(unlist(ZtL),ncol=3,byrow=TRUE)
a=unlist(tapply(Zt[,2],INDEX=Zt[,3],FUN=diff))
b=unlist(tapply(Zt[,1],INDEX=Zt[,3],FUN=tail,n=-1))
c=unlist(tapply(Zt[,1],INDEX=Zt[,3],FUN=head,n=-1))
d=a^alpha[c]
e=(c!=b)
gg=unlist(tapply(Zt[,3],INDEX=Zt[,3],FUN=head,n=-1))
Zt=cbind(c,b,a,d,e,gg)
TRAN=table(Zt[,1],Zt[,2])
SOJ=tapply(Zt[,3],INDEX=Zt[,1],FUN=sum)
SOJ.alpha=tapply(Zt[,4],INDEX=Zt[,1],FUN=sum)
LSOJ.nocen=tapply(log(Zt[,3])*Zt[,5],INDEX=Zt[,1],FUN=sum)
n.observed=rowSums(TRAN)-diag(TRAN)

#parametrizzazione alpha eta
eta=rgamma(na.states,shape=n.observed[1:na.states]+f,rate=SOJ.alpha[1:na.states]+g)
for ( s in 1:na.states) P[s,-s]=rdirichlet(1,alpha=(TRAN[s,-s]+alpha.dir[s,-s]))


if (SM==1){
for (s in 1:na.states){
alpha.prop=rnorm(1,alpha[s],0.25)
if (alpha.prop<0) alpha.prop=alpha[s]
SOJ.alpha.prop=tapply(a^alpha.prop,INDEX=Zt[,1],FUN=sum)
lnum=-eta[s]*SOJ.alpha.prop[s]+  n.observed[s]*(log(alpha.prop))+alpha.prop*LSOJ.nocen[s]-log(alpha.prop)+dnorm(log(alpha.prop),0,1,log=TRUE)
lden=-eta[s]*SOJ.alpha[s]     +  n.observed[s]*(log(alpha[s]))   +alpha[s]*LSOJ.nocen[s]-log(alpha[s])+dnorm(log(alpha[s]),0,1,log=TRUE)
if(runif(1)<exp(lnum-lden)) alpha[s]=alpha.prop}}

for (s in 1:na.states) {
    gamma[s]=eta[s]^(1/alpha[s])
    RateM[s,]=P[s,]*gamma[s]}
diag(RateM)=-gamma


param[iter,]=c(t(RateM),alpha,sum(acc)/Nobs)
print(iter)
#print(c(iter,gamma,alpha,sum(acc)/Nobs))

}

return(param)}
