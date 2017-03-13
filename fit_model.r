#--------------------------------------------------#
#--------------- read in data ---------------------#
#--------------------------------------------------#

# 0 = non-detection
# 1 = detection of a female without young
# 2 = detection of a female with one young
# 3 = detection of a female with two youngs
# 4 = detection of a female with unknown number of youngs
datahcr <- read.table("deerupdatedentree2y.txt",header=F) # last column is censored index 
head(datahcr)

# for simplicity, remove censored individuals
ind <- which(datahcr[,ncol(datahcr)]==-1)
datahcr <- datahcr[-ind,] 
datahcr <- datahcr[,-ncol(datahcr)] # delete last column
head(datahcr)

# number of individuals 
n <- dim(datahcr)[[1]] 

# number of capture occasions
K <- dim(datahcr)[[2]] 

# compute date of first capture
e <- NULL
for (i in 1:n){
	temp <- 1:K
	e <- c(e,min(temp[datahcr[i,]!=0]))
}

#----------------------------------------------#
#--------------- fit model --------------------#
#----------------------------------------------#

# compute number of states dropping structural 0s (typically before first capture)
lgr = NULL
for (i in 1:n){lgr = c(lgr,length(e[i]:K))}
lgr = 1:sum(lgr)
index = matrix(NA,nrow=n,ncol=K)
ind=1
for (i in 1:n){
	for (j in e[i]:K){
		index[i,j] = lgr[ind]
		ind = ind+1
	}
}
index

# calculate age
age <- matrix(NA,nrow=nrow(datahcr),ncol=ncol(datahcr))
for (i in 1:nrow(age)){
	for (j in e[i]:ncol(age)){
		if ((j >= e[i]) & (j < e[i]+6)) age[i,j] <- 1
		if (j >= e[i]+6) age[i,j] <- 2
	}
}

# write model
sink("model.bug")
cat("

# state-space formulation of multievent capture-recapture models
# O. Gimenez, March 2017.

# OBSERVATIONS
# 0 = non-detection
# 1 = detection of a female without young
# 2 = detection of a female with one young
# 3 = detection of a female with two youngs
# 4 = detection of a female with unknown number of youngs

# STATES
# NB is for non-breeder
# B1 is for breeder with 1
# B2 is for breeder with 2
# D dead

# PARAMETERS
# phiNB  survival prob. of non-breeding individuals  / by age
# phiB  survival prob. of breeding individuals / by age
# pNB  detection prob. of NB individuals 
# pB  detection prob. of B individuals 
# psiNBB1 transition prob. from NB to B1
# psiNBB2 transition prob. from NB to B2
# psiB1NB transition prob. from B1 to NB
# psiB1B2 transition prob. from B1 to B2
# psiB2NB transition prob. from B2 to NB
# psiB2B1 transition prob. from B2 to B1
# piNB prob. of being in initial state NB
# piB1 prob. of being in initial state B1

model
{

	################################
	# STATE-SPACE MODEL LIKELIHOOD #
	################################
	
	# probabilities for each initial state
		px0[1] <- pi1 # prob. of being in initial state NB
		px0[2] <- pi2 # prob. of being in initial state B1
		px0[3] <- 1-pi1-pi2 # prob. of being in initial state B2
		px0[4] <- 0 # prob. of being in initial state dead
	
	# define probabilities of O(t) given X(t)
		po[1,1] <- 1-pNB
		po[1,2] <- pNB * delta
		po[1,3] <- 0
		po[1,4] <- 0
		po[1,5] <- pNB * (1-delta)

		po[2,1] <- 1-pB
		po[2,2] <- 0
		po[2,3] <- pB * deltap
		po[2,4] <- 0
		po[2,5] <- pB * (1-deltap)

		po[3,1] <- 1-pB
		po[3,2] <- 0
		po[3,3] <- 0
		po[3,4] <- pB * deltap
		po[3,5] <- pB * (1-deltap)

		po[4,1] <- 1
		po[4,2] <- 0
		po[4,3] <- 0
		po[4,4] <- 0
		po[4,5] <- 0

  # define probabilities of O(t) given X(t)
		po.init[1,1] <- 0
		po.init[1,2] <- delta
		po.init[1,3] <- 0
		po.init[1,4] <- 0
		po.init[1,5] <- 1-delta

		po.init[2,1] <- 0
		po.init[2,2] <- 0
		po.init[2,3] <- deltap
		po.init[2,4] <- 0
		po.init[2,5] <- 1-deltap

		po.init[3,1] <- 0
		po.init[3,2] <- 0
		po.init[3,3] <- 0
		po.init[3,4] <- deltap
		po.init[3,5] <- 1-deltap

		po.init[4,1] <- 1
		po.init[4,2] <- 0
		po.init[4,3] <- 0
		po.init[4,4] <- 0
		po.init[4,5] <- 0

	for (i in 1:N)  # for each indiv
	{
	
	# estimated probabilities of initial states are the proportions in each state at first capture occasion
	alive[i,First[i]] ~ dcat(px0[1:4])
    mydata[i,First[i]] ~ dcat(po.init[alive[i,First[i]],1:5])

		for (j in (First[i]+1):Years)  # loop over time

		{

		# define probabilities of X(t+1) given X(t) / 1er indice = ligne
		
		# if ageij = 1, then phi = phi-prime-aged
		# if ageij = 2, then phi = phi-old
		phiNB[i,j-1] <- phiNBpa * equals(age[i,j-1],1) + phiNBold * equals(age[i,j-1],2)
		phiB[i,j-1] <- phiBpa * equals(age[i,j-1],1) + phiBold * equals(age[i,j-1],2) 
		
		# probabilities of observations given states and states given states
		px[1,i,j-1,1] <- phiNB[i,j-1] * 1/(1+exp(alpha[1,1])+exp(alpha[1,2]))
		px[1,i,j-1,2] <- phiNB[i,j-1] * exp(alpha[1,1])/(1+exp(alpha[1,1])+exp(alpha[1,2]))
		px[1,i,j-1,3] <- phiNB[i,j-1] * exp(alpha[1,2])/(1+exp(alpha[1,1])+exp(alpha[1,2]))
		px[1,i,j-1,4] <- 1-phiNB[i,j-1]
		
		px[2,i,j-1,1] <- phiB[i,j-1] * exp(alpha[2,1])/(1+exp(alpha[2,1])+exp(alpha[2,2]))
		px[2,i,j-1,2] <- phiB[i,j-1] * 1/(1+exp(alpha[2,1])+exp(alpha[2,2]))
		px[2,i,j-1,3] <- phiB[i,j-1] * exp(alpha[2,2])/(1+exp(alpha[2,1])+exp(alpha[2,2]))
		px[2,i,j-1,4] <- 1-phiB[i,j-1] 
		
		px[3,i,j-1,1] <- phiB[i,j-1] * exp(alpha[3,1])/(1+exp(alpha[3,1])+exp(alpha[3,2]))
		px[3,i,j-1,2] <- phiB[i,j-1] * exp(alpha[3,2])/(1+exp(alpha[3,1])+exp(alpha[3,2]))
		px[3,i,j-1,3] <- phiB[i,j-1] * 1/(1+exp(alpha[3,1])+exp(alpha[3,2]))
		px[3,i,j-1,4] <- 1-phiB[i,j-1] 

		px[4,i,j-1,1] <- 0
		px[4,i,j-1,2] <- 0
		px[4,i,j-1,3] <- 0
		px[4,i,j-1,4] <- 1

		## STATE EQUATIONS ##
		# draw X(t) given X(t-1)
		alive[i,j] ~ dcat(px[alive[i,j-1],i,j-1,1:4])

		## OBSERVATION EQUATIONS ##
		# draw O(t) given X(t)
		mydata[i,j] ~ dcat(po[alive[i,j],1:5])

		}
	
	}

##########
# PRIORS #
##########

pNB ~ dunif(0, 1) # non-breeder detectability
pB ~ dunif(0, 1) # breeder detectability

phiBpa ~ dunif(0, 1) # breeder survival prime-aged
phiBold ~ dunif(0, 1) # breeder survival old
phiNBpa ~ dunif(0, 1) # non-breeder survival prime-aged
phiNBold ~ dunif(0, 1) # non-breeder survival old

delta ~ dunif(0, 1) # assignement
deltap ~ dunif(0, 1) # assignement

pi1 ~ dunif(0, 1) # prop
pi2 ~ dunif(0, 1) # prop

# transition probabilites - multinomial logit
for (i in 1:3){
	for (j in 1:2){
		alpha[i,j] ~ dnorm(0,0.1)
			      }
			  }
	
for (i in 1:N){
	for (j in First[i]:Years){
		lrs[indicateur[i,j]] <- alive[i,j]
		                     }
	          }
	
}

",fill=TRUE)
sink()

# data
mydatax <- list(N=n,Years=K,mydata=as.matrix(datahcr+1),First=e,indicateur=index,age=as.matrix(age)) 

# inits for hidden states
mask = (datahcr==4)
alive = datahcr
alive[mask] <- rbinom(1,1,0.5)+2 # assign uncertain to B1 or B2 with odds 50:50
for (i in 1:n) {
	for (j in 1:K) {
		if (j > e[i] & datahcr[i,j]==0) {alive[i,j] = which(rmultinom(1, 1, c(1/3,1/3,1/3))==1)}
		if (j < e[i]) {alive[i,j] = NA}
	}
}
alive <- as.matrix(alive)

# list of inits
init1 <- list(pNB=0.5,pB=0.5,phiNBpa=0.5,alive=alive,delta=0.5,deltap=0.5,pi1=0.6,pi2=0.2)
init2 <- list(pNB=0.8,pB=0.8,phiNBpa=0.8,alive=alive,delta=0.8,deltap=0.8,pi1=0.6,pi2=0.2)
inits <- list(init1,init2)

# specify the parameters to be monitored
parameters <- c("phiBpa","phiBold","phiNBpa","phiNBold","pi1","pi2","pB","pNB","delta","deltap","alpha","lrs")

# run model (warning: takes more than 2 hours)
library(rjags)
start<-as.POSIXlt(Sys.time())
jmodel <- jags.model("model.bug", mydatax, inits, n.chains = 2,n.adapt = 5000)
jsample <- coda.samples(jmodel, parameters, n.iter=20000, thin = 1)
end <-as.POSIXlt(Sys.time())
duration = end-start
save(jsample,jmodel,duration,file='fit_deer.Rdata')

#----------------------------------------------#
#--------------- post-process -----------------#
#---------------    results   -----------------#
#----------------------------------------------#

# proportion param
pi1 <- c(jsample[[1]][,'pi1'],jsample[[2]][,'pi1'])
mean(pi1)
sd(pi1)

pi2 <- c(jsample[[1]][,'pi2'],jsample[[2]][,'pi2'])
mean(pi2)
sd(pi2)

# survival parameters
phiNBpa <- c(jsample[[1]][,'phiNBpa'],jsample[[2]][,'phiNBpa'])
mean(phiNBpa)
sd(phiNBpa)

phiNBold <- c(jsample[[1]][,'phiNBold'],jsample[[2]][,'phiNBold'])
mean(phiNBold)
sd(phiNBold)

phiBpa <- c(jsample[[1]][,'phiBpa'],jsample[[2]][,'phiBpa'])
mean(phiBpa)
sd(phiBpa)

phiBold <- c(jsample[[1]][,'phiBold'],jsample[[2]][,'phiBold'])
mean(phiBold)
sd(phiBold)

# observation parameters
pNB <- c(jsample[[1]][,'pNB'],jsample[[2]][,'pNB'])
mean(pNB)
sd(pNB)

pB <- c(jsample[[1]][,'pB'],jsample[[2]][,'pB'])
mean(pB)
sd(pB)

delta <- c(jsample[[1]][,'delta'],jsample[[2]][,'delta'])
mean(delta)
sd(delta)

deltap <- c(jsample[[1]][,'deltap'],jsample[[2]][,'deltap'])
mean(deltap)
sd(deltap)

# back-transform transition probabilities estimates
a12 <- c(jsample[[1]][,'alpha[1,1]'],jsample[[2]][,'alpha[1,1]'])
a13 <- c(jsample[[1]][,'alpha[1,2]'],jsample[[2]][,'alpha[1,2]'])
a11 <- rep(0,length(a12)) ## ref
a1<-cbind(a11,a12,a13)
# psiNBNB psiNBB1 psiNBB2
apply(exp(a1)/apply(exp(a1),1,sum),2,mean)
a21 <- c(jsample[[1]][,'alpha[2,1]'],jsample[[2]][,'alpha[2,1]'])
a23 <- c(jsample[[1]][,'alpha[2,2]'],jsample[[2]][,'alpha[2,2]'])
a22 <- rep(0,length(a21)) ## ref
a2<-cbind(a21,a22,a23)
# psiB1NB psiB1B1 psiB1B2
apply(exp(a2)/apply(exp(a2),1,sum),2,mean)
a31 <- c(jsample[[1]][,'alpha[3,1]'],jsample[[2]][,'alpha[3,1]'])
a33 <- c(jsample[[1]][,'alpha[3,2]'],jsample[[2]][,'alpha[3,2]'])
a32 <- rep(0,length(a31)) ## ref
a3<-cbind(a31,a32,a33)
# psiB2NB psiB2B1 psiB2B2
apply(exp(a3)/apply(exp(a3),1,sum),2,mean)

# get transition parameter standard deviations
apply(exp(a1)/apply(exp(a1),1,sum),2,sd)
apply(exp(a2)/apply(exp(a2),1,sum),2,sd)
apply(exp(a3)/apply(exp(a3),1,sum),2,sd)

# states estimates
res <- rbind(jsample[[1]],jsample[[2]])
dim(res) # 40000 x 4499
lrs <- res[,9:4491]
dim(lrs) # 40000 x 4483
nrowarray = dim(lrs)[1]
lrsind = array(NA,c(nrowarray,n,K))
for (k in 1:nrowarray){
	ind=1
	for (i in 1:n){
		for (j in e[i]:K){
			lrsind[k,i,j] = lrs[k,index[i,j]]
			ind = ind+1
		}
	}
}
dim(lrsind)
# [1] 40000   211    35
save(lrsind,file='states.Rdata')
