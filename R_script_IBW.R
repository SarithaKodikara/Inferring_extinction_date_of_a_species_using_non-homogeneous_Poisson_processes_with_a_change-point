source("DBDA2E-utilities.R")
require(rjags)               # Must have previously installed package rjags.
library(coda)
set.seed(1234)
fileNameRoot="Jags-ExampleScript" # For output file names.
sighting<- c(1,2,3,4,5,7,8,9,10,11,12,13,14,16,17,19,20,23,24,26,27,28,29,32,33,34,35,36,37,38,39,40,41,42,44,45,46,47,49,51,52,53,54,55,58,61,62,65,69,70,71,72,74,75,76,77,79,84,85,88,89,90,91,102,107,108,109)
cat<- c(0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,1,0,1,1,1,0,0,1,1,1,1,0,1,1,0,1,1,0,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1)
Tt<-113
t_n<-42
N<-length(sighting)
y<-c()

dataList = list(    # Put the information into a lis
  d=sighting,
  y=y,
  t_n=t_n,
  N=N,
  cat=cat,
  Tt=Tt)

# Define the model:
modelpois = paste0("
                        data {
                        C <- 10000000000000 # JAGS does not warn if too small!
                        
                        
                        for (i in 1:N) {
                          ones[i] <- 1
                        }
                        
                        
                        }
                        
                        model {
                        
                        
                        for(i in 1:N){
                          y[i]<- ifelse(cat[i]<1,
                          ifelse(tau<=t_n,(10^(-100)), ifelse(tau<=Tt, ((a/(sigma^a))*d[i]^(a-1)*exp(-1*((tau/sigma)^a))),((a/(sigma^a))*d[i]^(a-1)*exp(-1*((Tt/sigma)^a))))),
                          ifelse(tau<=t_n,(10^(-100)), ifelse(d[i]<=tau, ((a_U1/(sigma_U1^a_U1))*d[i]^(a_U1-1)*exp(-1*((tau/sigma_U1)^a_U1))) ,  ((a_U2/(sigma_U2^a_U2))*d[i]^(a_U2-1)*exp(-1*((Tt/sigma_U2)^a_U2- (tau/sigma_U2)^a_U2))))))/C
                          
                         
                          
                          ones[i]~ dbern( y[i] )
                        }
                        
                        tau ~ dexp(0.005) 
                        sigma~dunif(0,1000)
                        a~dunif(0,1000)
                        sigma_U1~dunif(0,1000)
                        a_U1~dunif(0,1000)
                        sigma_U2~dunif(0,1000)
                        a_U2~dunif(0,1000)
                        }
                        
                        ") # close quote for modelString



writeLines( modelpois , con="model_pois.txt" )

# Run the chains:
jagsmodelpois = jags.model( file="model_pois.txt" , data=dataList, n.chains=4 , n.adapt=10000)
update( jagsmodelpois , n.iter=10000 )
codaSamplespois= coda.samples( jagsmodelpois, variable.names=c("tau","a","sigma","a_U1","sigma_U1","a_U2","sigma_U2") ,n.iter=130000 , thin = 13)
save( codaSamplespois , file=paste0(fileNameRoot,"Mcmc_pois.Rdata") )


diagMCMC( codaObject=codaSamplespois  , parName="tau" )
diagMCMC( codaObject=codaSamplespois , parName="a" )
diagMCMC( codaObject=codaSamplespois , parName="sigma" )
diagMCMC( codaObject=codaSamplespois , parName="a_U1" )
diagMCMC( codaObject=codaSamplespois , parName="sigma_U1" )
diagMCMC( codaObject=codaSamplespois , parName="a_U2" )
diagMCMC( codaObject=codaSamplespois , parName="sigma_U2" )

ma <- matrix(nrow=10000, ncol=2)
colnames(ma) <- c("tau","tau_E")

prob=list(ma, ma, ma, ma)
for(i in 1:4){
  for(j in 1:10000){
    prob[[i]][j,1] = codaSamplespois[[i]][j,7]
    prob[[i]][j,2] = codaSamplespois[[i]][j,7]+1897
  }
}

newCodaSamples= mcmc.list(as.mcmc(prob[[1]][,]),as.mcmc(prob[[2]][,]),as.mcmc(prob[[3]][,]),as.mcmc(prob[[4]][,])) 
save( newCodaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") ) 


openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost(codaSamplespois[,"tau"], main="" ,xlim=c(40,80), ylim=c(0,0.1), xlab=bquote(tau[E]),showCurve=TRUE, cenTend=FALSE)


openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost(newCodaSamples[,"tau_E"], main="" ,xlim=c(1897+40,1897+80), ylim=c(0,0.1), xlab=bquote(tau[E]),showCurve=TRUE, cenTend=FALSE)

openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamplespois[,"a"], main="" ,xlim=c(0,2), xlab=bquote(alpha[c]),showCurve=TRUE)


openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamplespois[,"a_U1"], main="" ,xlim=c(.5,7), xlab=bquote(alpha[u1]),showCurve=TRUE)

openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamplespois[,"a_U2"], main="" ,xlim=c(0,6), xlab=bquote(alpha[u2]),showCurve=TRUE)


