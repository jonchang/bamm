#A -1.5
#B  1

root <- -0.19421
dl <- -1.5 - root
dr <- 1 - root

rtmean <- 0
lfmean <- 0.574164

rtvar <- 5
lfvar <- 3.18513

dl
dr

dnorm(dl, lfmean, sqrt(lfvar), log=T) + dnorm(dr, rtmean,sqrt(rtvar), log=T)

 
 


# Generate whale data

library(mvtnorm)
v <- read.tree('whaletree.tre')
vx <- vcv.phylo(v) * 0.1
trait <- as.vector(rmvnorm(1, sigma=vx))
names(trait) <- rownames(vx)
trait[53:87] <- trait[53:87] - 10
dff <- data.frame(names(trait), trait) 
write.table(dff, file = "whaletrait", sep="\t", quote=F, row.names=F, col.names=F)



tx <- trait
tx[53:87] <- tx[53:87] - 5

mc <- read.csv('mcmc_out.txt') 
mc <- mc[1000:nrow(mc), ]
table(mc$N_shifts) / nrow(mc)

xx <- read.csv("event_data.txt") 
xx <- xx[xx$betainit == 0, ]

x <- read.csv("traitmodeldata.txt")


