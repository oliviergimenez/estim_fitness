#--------------------- read in data

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

# load estimated states
load('states.Rdata')

#--------------------- compute LRS à la McGraw&Caswell

# function to calculate individual fitness Ã  la McGraw & Caswell 1996
# INPUT: a vector of reproductive events (nb of offsprings/recruits) per year (e.g.)
# OUTPUT: individual fitness calculated as dominant eigenvalue of individual Leslie matrix
# USE: see two examples by McGraw&Caswell
# tit<-c(3.5, 5, 5, 6, 4.5)*2
# ind_fitness(tit)
# sparrowhawks <- c(0, 0, 0, 1.5, 0, 0, 1, 1.5 )*2
# ind_fitness(sparrowhawks)
ind_fitness <- function(y) {
nb.age.classes <- length(y) # get number of age classes
if (nb.age.classes==1) return(0)
if (nb.age.classes>1){
# build up individual Leslie matrix
ind.leslie.matrix <- matrix(0,nrow=nb.age.classes,ncol=nb.age.classes)
ind.leslie.matrix[1,] <- y/2
for (i in 2:nb.age.classes){ind.leslie.matrix[i,i-1] <- 1}       
#return(ind.leslie.matrix)
ind.fitness <- eigen(ind.leslie.matrix,only.values = T)$values # calculate eigenvalues
ind.fitness <- max(Re(ind.fitness)) # get real part if complex number and dominant eigenvalue
return(ind.fitness)}
}

nbMCMC = 40000
res <- matrix(NA,nrow=nbMCMC,ncol=n)
for (k in 1:nbMCMC){
	for (i in 1:n){
		# look for date of death
		DD = min(which(lrsind[k,i,e[i]:K]==4))
		if (DD!=Inf) {
			if (DD==2) res[k,i] <- ind_fitness(lrsind[k,i,e[i]]-1) # si la mort a lieu tout de suite aprÃ¨s premiÃ¨re capture
			if (DD>2) res[k,i] <- ind_fitness(lrsind[k,i,e[i]:(e[i]+DD-2)]-1) # si la mort a lieu au moins une occasion aprÃ¨s premiÃ¨re capture
		}
		if (DD==Inf) res[k,i] <- ind_fitness(lrsind[k,i,e[i]:K]-1)
	}
}
t(res)
# warnings are for individuals that did not die over the study period
# DD is Inf in this situation
dim(res) # nbMCMC x nbindiv
save(res,file='indfitdeer.Rdata')

# plot distribution of (posterior mean) individual fitness
ppi <- 300
png("indfit.png", width=6*ppi, height=6*ppi, res=ppi)
plot(density(apply(res,2,mean)),xlab='',ylab='',main='Fitness a la McGraw and Caswell',xlim=c(0,2.6))
dev.off()

#--------------------- compute standard LRS, à la Clutton-Brock

lrstot <- matrix(0,nrow=n,ncol=dim(lrsind)[1])
for (kk in 1:dim(lrsind)[1]){
  for (i in 1:n){
    for (j in e[i]:K){
      if (lrsind[kk,i,j] == 1) lrstot[i,kk] <- lrstot[i,kk] + 0
      if (lrsind[kk,i,j] == 2) lrstot[i,kk] <- lrstot[i,kk] + 1
      if (lrsind[kk,i,j] == 3) lrstot[i,kk] <- lrstot[i,kk] + 2
      if (lrsind[kk,i,j] == 4) lrstot[i,kk] <- lrstot[i,kk] + 0
    }
  }
}		

# compute LRS, credible intervals and standard deviation
lrsfinal <- apply(lrstot,1,median)
lrs.lb <- apply(lrstot,1,quantile,probs=2.5/100)
lrs.ub <- apply(lrstot,1,quantile,probs=97.5/100)
lrs.sd <- apply(lrstot,1,sd)

save(lrstot,lrsfinal,file='lrstot.Rdata')

ppi <- 300
png("lrs.png", width=6*ppi, height=6*ppi, res=ppi)
hist(lrsfinal,ylim=c(0,100),xlab="Lifetime reproductive success", ylab="Number of individuals",main='Fitness a la Clutton-Brock')
dev.off()



#--------------------- compute delifing, à la Coulson

# R code to implement Coulson's 2006 delifing method
# Gimenez, August 2012
# Estimating Individual Contributions to Population Growth: Evolutionary Fitness in Ecological Time
# T. Coulson, T.G. Benton, P. Lundberg, S.R.X. Dall, B.E. Kendall, J-M. Gaillard


# R equivalent to diag(x,k) in matlab
subdiag <-function(vec, size, offset=0){ 
M<-matrix(0,size,size)
M[row(M)-offset==col(M)]<-vec
return(M)}

# zero machine
eps = 2.2204e-016

# R equivalent to diag(x,k) in matlab
delife <-function(datax){ 

dat <- datax

# misc calculations
# create sorted list of unique names
names=sort(unique(dat[,1]))
length(names)
# create some useful constants
first_year = min(dat[,2])
last_year = max(dat[,2])
num_years = last_year-first_year+1
max_age = max(dat[,4])
number_individuals = length(names)

#--------------------------------------------------------
# STEP 1: calculate the annual population projection (transition) matrices, annual vectors of nb of individuals and wt per capita growth rate
#--------------------------------------------------------

# creates a three dimensional matrix for annual matrices
annual_mean_matrices = array(0,c(max(dat[,4]),max(dat[,4]),(max(dat[,2])-min(dat[,2])+1)))
annual_trans_vectors = matrix(0,max(dat[,4]),(max(dat[,2])-min(dat[,2])+1))
lambda <- rep(0,num_years)

for (years in first_year:last_year){ 

    # select data for each year
    dat2=dat; dat2=dat2[dat2[,2]==years,]
    
    if (is.matrix(dat2)!=T) next
    
    # calculate population age structure, n, and age-specific survival rates, s and recruitment rates, r 
    n <- rep(0,length=max_age)
    s <- rep(0,length=max_age)
    r <- rep(0,length=max_age)
    for (j in 1:max_age){ # this is 'data' and not 'dat' to ensure I always get the same sized matrix 'mat'
        n[j]=sum(dat2[,4]==j); s[j]=sum(dat2[,5]*(dat2[,4]==j)); r[j]=sum(dat2[,6]*(dat2[,4]==j))
    }
    # generate the transition matrix 
    n=n+(n==0)*eps # assign almost-0 to 0 values to avoid NAs
    s=s/n; s=s[1:length(s)-1]; mat=subdiag(s,length(s)+1,1); mat[1,]=r/n # fill in matrix
    n1 = mat%*%n
    lambda[years-min(dat[,2])+1] = sum(n1)/sum(n) # per capita growth rate wt
    annual_mean_matrices[,,(years-first_year+1)] = mat # pop projection matrix
    annual_trans_vectors[,(years-first_year+1)] = n # population vector
}


#--------------------------------------------------------
# STEP 2: remove each individual from the data set in turn (by deleting one row at a
# time) and recalculate the demographic matrix, population vector and population growth
#--------------------------------------------------------

indiv_matrices = array(0,c(max_age,max_age,num_years,number_individuals))
indiv_vectors = array(0,c(max_age,num_years,number_individuals))
lambda_ind = matrix(0,number_individuals,num_years)

# loop through animal
for (animal in 1:number_individuals){ 
    data1 = dat; data1=data1[data1[,1]!=names[animal],] # remove current individual i
#    loop through years for each animal
    for (years in first_year:last_year){
    	dat2=data1; dat2=subset(dat2,dat2[,2]==years)  # select appropriate subset of data for individual i and year t
        for (j in 1:max_age){
        	n[j]=sum(dat2[,4]==j); s[j]=sum(dat2[,5]*(dat2[,4]==j)); r[j]=sum(dat2[,6]*(dat2[,4]==j))
        }
        n=n+(n==0)*eps
        s=s/n; s=s[1:length(s)-1]; mat=subdiag(s,length(s)+1,1); mat[1,]=r/n # fill in matrix
        n1 = mat%*%n
        lambda_ind[animal,years-first_year+1] = sum(n1)/sum(n) # growth rate wt(-i) (indiv i removed)
        indiv_matrices[,,years-first_year+1,animal] = mat # pop matrix At(-i) (indiv i removed)
        indiv_vectors[,years-first_year+1,animal] = n # vector Nt(-i) (indiv i removed)
    }
}
# note: mean of the wt(-i)s for individuals within the population in year t is equal to wt

#--------------------------------------------------------
# STEP 3: calculate individual contributions to population growth
# pi = wt- wt(-i)
# a positive value of p_i reflects an individual that outperformed the mean individual 
# contribution to population growth and a negative value reflects an individual that underperformed.
#--------------------------------------------------------

lambda_new <- matrix(0,nrow=number_individuals,ncol=length(lambda))
for (i in 1:number_individuals){
    lambda_new[i,] = lambda - lambda_ind[i,] # indiv contrib of indiv i = p_i
}
lambda_ind = lambda_new

#--------------------------------------------------------
# STEP 4: estimate how sensitive population growth is to each individual
# The difference between At - At(-i) = termed St(i) produces a perturbation matrix for each individual
# these values describe how each individual contributed to population growth through survival and recruitment
# A negative element in St(i) means that the individual performed less well
# than average individuals of that age, while a positive value means the individual
# performed better than the average.
#--------------------------------------------------------

new <- array(0,c(max_age,max_age,num_years,length(names)))
for (i in 1:length(names)){
    for (j in 1:num_years){
        new[,,j,i] = annual_mean_matrices[,,j] - indiv_matrices[,,j,i] # the St(i)
    }
}

indiv_matrices = new

#--------------------------------------------------------
# STEP 5: store results in dat
#--------------------------------------------------------

dat <- cbind(datax,0,0,0)

for (i in 1:dim(dat)[1]){
    dat[i,7] = lambda_ind[dat[i,1]==names,dat[i,2]-first_year+1] # individual contrib pt(-i) to growth rate
    age_col = dat[i,4]
    age_row = age_col+1
    if (age_row==max_age+1) age_row=max_age
    year = dat[i,2]-first_year+1
    id = dat[i,1]==names
    dat[i,8] = indiv_matrices[1,age_col,year,id] # contrib through survival?
    dat[i,9] = indiv_matrices[age_row,age_col,year,id] # contrib through recruitment?
}

# moyenne par individu

return(dat)

}

year <- 1977:2011

delife.tot <- matrix(0,nrow=n,ncol=dim(lrsind)[1])

for (kk in 1:dim(lrsind)[1]){

#----- FORMAT DATA

# data must be in following format
# col 1 = id
# col 2 = year
# col 3 = sex
# col 4 = age
# col 5 = survival between t and t+1
# col 6 = recruitment between t and t+1

data.format <- NULL
ex <- lrsind[kk,,]

  for (i in 1:n){
    for (j in e[i]:K){
    	id.tp <- i
    	year.tp <- year[j]
    	sex.tp <- 1 # female only
    	age.tp <- j - e[i] + 1   	    	
    	if (ex[i,j] != 4)
    	{
    		surv.tp <- 1 
    		recruit.tp <- ex[i,j]-1
    		}
    	else 
    	{
    		surv.tp <- 0
    		recruit.tp <- 0
    		}
    	temp <- c(id.tp,year.tp,sex.tp,age.tp,surv.tp,recruit.tp)
    	data.format <- rbind(data.format,temp)
       }
  }

#----- DELIFE!

resbig <- delife(data.format)
res <- resbig[,c(1,2,7)]
id <- sort(unique(res[,1]))
pfinal <- rep(0,length(id))
ind <- 1
for (i in id)
{
pfinal[ind] <- mean(res[res[,1]==i,3])
ind <- ind+1
     }

delife.tot[,kk] <- pfinal

}		

save(delife.tot,resbig,file='delife.Rdata')

ppi <- 300
png("delifing.png", width=6*ppi, height=6*ppi, res=ppi)
plot(density(apply(delife.tot,1,mean)),xlab='',ylab='',main='Lifetime individual contribution (LIC)')
dev.off()

# all figures in a single one (Figure 1 in the paper)
ppi <- 300
png("fig1.png", width=6*ppi, height=6*ppi, res=ppi)
par(mfrow=c(3,1))
hist(lrsfinal,ylim=c(0,150),xlab="", ylab="",main=expression('a) Lifetime reproductive success (LRS)'))
plot(density(apply(res,2,mean)),xlab='',ylab='',main=expression(paste('b) Individual fitness (',lambda[ind],')',sep="")),xlim=c(0,2.6))
plot(density(apply(delife.tot,1,mean)),xlab='',ylab='',main=expression('c) Lifetime individual contribution (LIC)'))
dev.off()
