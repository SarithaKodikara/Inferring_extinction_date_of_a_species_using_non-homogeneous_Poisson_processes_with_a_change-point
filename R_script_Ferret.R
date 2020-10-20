source("DBDA2E-utilities.R")
require(rjags)               # Must have previously installed package rjags.
library(coda)
set.seed(1234)
fileNameRoot="Jags-ExampleScript" # For output file names.
sighting<-c(6/12, 7/12, 8/12, 10/12, 1+5/12, 1+6/12, 1+7/12,1+8/12, 1+9/12, 1+10/12, 2+6/12, 2+7/12, 3+5/12, 3+8/12, 3+10/12, 4+5/12, 4+9/12, 4+10/12, 5+6/12, 7+6/12, 9+9/12, 9+10/12, 10+2/12, 10+3/12, 10+7/12, 11+7/12, 12+7/12, 12+9/12)
d<-sighting*12
Tt<-229
t_n<-153
N<-length(d)
y<-c()

dataList = list(    # Put the information into a lis
  d=d,
  y=y,
  t_n=t_n,
  N=N,
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
                          
                          y[i]<-ifelse(tau<=t_n,(10^(-100)), ifelse(tau<=Tt, ((a/(sigma^a))*d[i]^(a-1)*exp(-1*((tau/sigma)^a))),
                          ((a/(sigma^a))*d[i]^(a-1)*exp(-1*((Tt/sigma)^a)))))
                         
                         
                          
                          ones[i]~ dbern( y[i] )
                        }
                        
                        tau ~ dexp(0.005) 
                        sigma~dunif(0,1000)
                        a~dunif(0,1000)
                        
                        }
                        
                        ") # close quote for modelString



writeLines( modelpois , con="model_pois.txt" )

# Run the chains:
jagsmodelpois = jags.model( file="model_pois.txt" , data=dataList, n.chains=4 , n.adapt=10000)
update( jagsmodelpois , n.iter=10000)
codaSamplespois= coda.samples( jagsmodelpois, variable.names=c("tau","a","sigma") ,n.iter=3600000 , thin =360)
save( codaSamplespois , file=paste0(fileNameRoot,"Mcmc_pois.Rdata") )


diagMCMC( codaObject=codaSamplespois  , parName="tau")
diagMCMC( codaObject=codaSamplespois , parName="a" )
diagMCMC( codaObject=codaSamplespois , parName="sigma" )

ma <- matrix(nrow=10000, ncol=2)
colnames(ma) <- c("tau","tau_E")

prob=list(ma, ma, ma, ma)
for(i in 1:4){
  for(j in 1:10000){
    prob[[i]][j,1] = codaSamplespois[[i]][j,3]
    prob[[i]][j,2] = codaSamplespois[[i]][j,3]/12+1972
  }
}

newCodaSamples= mcmc.list(as.mcmc(prob[[1]][,]),as.mcmc(prob[[2]][,]),as.mcmc(prob[[3]][,]),as.mcmc(prob[[4]][,])) 
save( newCodaSamples , file=paste0(fileNameRoot,"Mcmc.Rdata") ) 


openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost(codaSamplespois[,"tau"], main="" ,xlim=c(40,80), ylim=c(0,0.1), xlab=bquote(tau[E]),showCurve=TRUE, cenTend=FALSE)


openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost(newCodaSamples[,"tau_E"], main="" ,xlim=c(1984,1995) , ylim=c(0,0.8), xlab=bquote(tau[E]),showCurve=TRUE, cenTend=FALSE)

openGraph(height=5,width=7)
par( mar=c(3.5,0.5,2.5,0.5) , mgp=c(2.25,0.7,0) )
plotPost( codaSamplespois[,"a"], main="" ,xlim=c(0,2), xlab=bquote(alpha[c]),showCurve=TRUE)




