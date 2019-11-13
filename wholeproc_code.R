load("dati_bbstest.RData")
source("modes.functions.r")
library(ks); library(ClusterR);library(snow); library(snowfall); library(ggplot2)
rm(hgridbs)
#####################
# #phase 1: variable selection
ngiri <- 1000
nvar <- 3
p <- rep(1/NCOL(datib),NCOL(datib))

keep <- matrix(NA,nrow=ngiri,ncol=nvar)
poss <- matrix(NA,nrow=ngiri,ncol=1)


for (i in 1:ngiri) {
  keep[i,] <- sample(1:ncol(datib),size=nvar,prob=p)
}

set.seed(12345)
sub <- matrix(sample(1:nrow(datib),nrow(datibs), rep=F), ncol=1)


sfInit(parallel=TRUE,cpus=3)
sfExport("datib","datibs","keep","poss","ngiri", "sub")
sfLibrary(ks)


poss <- sfSapply(1:ngiri, function(x) {
  hpib <- Hns(datib[sub,keep[x,]])
  poss[x] <- kde.test(datib[sub,keep[x,]],datibs[,keep[x,]],H1=hpib,H2=hpib)$pvalue
})



whichplow <- which(poss<0.01)
tot.examined <- table(keep)
tot.relevant <- table(keep[whichplow,])  

prop.relevant <- tot.relevant/tot.examined
barplot(prop.relevant)
save.image("finale.Rdata")

#####################
########################################################################
# phase 2: background density estimation

var.sel<- c(7,23)
datib.red <- datib[,var.sel] 
datib.red <- scale(datib.red)

datibs.red <- datibs[,var.sel]
datibs.red <- scale(datibs.red, attributes(datib.red)$`scaled:center`, attributes(datib.red)$`scaled:scale`)

plot(kde(datib.red))
plot(kde(datibs.red))
plot(kde(datitest.red))


clb <- kms(datib.red, tol.clust=0.05)

clb$nclust
clb$H
plot(kde(datib.red, H=clb$H))
save.image("finale.Rdata")


################################################
#phase 3: hbs selection
hgridbs <- c(0.010, 0.015, 0.035, 0.060, 0.070, 0.080, 0.090, 0.100, 0.120, 0.140, 0.160, 0.180,
 0.200, 0.250, 0.300, 0.350, 0.400, 0.450, 0.500)


fm <- nc <- numeric(length(hgridbs))
d<- 2
temp<- list()
print(c(date()))
for (i in 1:(length(hgridbs)-1)) {
  temp[[i]] <- kms(x = datibs.red, y=rbind(datib.red, datibs.red), H=diag(hgridbs[i],d))
  fm[i]<- external_validation(temp[[i]]$lab[1:20000], clb$lab, method = "fowlkes_mallows_index")
  nc[i] <- temp[[i]]$nclust
  print(c(i,date()))
  print(table(temp[[i]]$lab))
  print(fm[i])
  save.image("finale.Rdata")}

fm[length(hgridbs)] <- 1
nc[length(hgridbs)] <- 1

whichhopt <- which(fm==max(fm[nc>1]))
save.image("finale.Rdata")


# phase 5: testing significance of the signal mode
out.mode2        = modetest.fun3(datX=data.matrix(datibs.red), datY=data.matrix(datitest.red), 
                                 bw=hgridbs[whichhopt]^0.5, modes=matrix(temp[[whichhopt]]$mode[,2],ncol=2), alpha=0.0001, nboot = 1000)

modes      = out.mode2$modes
CI         = out.mode2$CIlead
CI
allCI      = out.mode2$allCI
save.image("finale.Rdata")

#pdf("test.pdf")
par(mar=c(3,2,1,1))
eigenport.fun(allci=allCI)
#dev.off()

#################################################
# phase 6: final estimation and classification
#################################################
# plot of hat f_bs
cont=c(10,20,30,44,60,80)
plot(kde(datibs.red, H=diag(hgridbs[whichhopt],2)), approx.cont=F, cont=cont)

#classification of X_bs
table(temp[[whichhopt]]$lab[20001:30000], label)
external_validation(temp[[whichhopt]]$lab[20001:30000], label, method = "fowlkes_mallows_index")

#classification of test data
clbstest <- kms(x = datibs.red, y=datitest.red, H=diag(hgridbs[whichhopt],2))
table(clbstest$lab, labeltest)
external_validation(clbstest$lab, labeltest, method = "fowlkes_mallows_index")
save.image("finale.Rdata")

