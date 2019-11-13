#code by Perone-Pacifico available at https://sites.google.com/a/uniroma1.it/marcoperonepacifico/R-code
#associated with the paper "Non-parametric inference for density modes" (Genovese et al., JRSSB, 2016)


library(matrixcalc)
library(mvtnorm)
library(scatterplot3d)


# kernel weights
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   n-dim vector of weights

hatfi.fun  = function(pt,dat,bw) {
   D       = ncol(dat)
   sqdis   = rowSums(t(pt-t(dat))^2) # square distance of point from data
   return( (2*pi*bw^2)^(-D/2)*exp(-sqdis/(2*bw^2)) )
   }

# kernel density estimate
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   kernel density estimate (scalar)

fhat.fun    = function(pt,dat,bw) {
   mean(hatfi.fun(pt,dat,bw))
   }

# mean shift
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   mean shift (D-dim vector)

ms.fun     = function(pt,dat,bw) {
   wi      = hatfi.fun(pt,dat,bw)
   colMeans(wi*dat)/(mean(wi))
   }

# mean shift sequence
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
#   eps  threshold for stopping rule
# output:
#   mode (D-dim vector)

msiter.fun = function(pt,dat,bw, eps=(bw^2)*10^(-7)) {
   start   = pt
   end     = ms.fun(start,dat,bw)
   while(sum(abs(end-start)) > eps){
      start=end
      end  = ms.fun(start,dat,bw)
   }
   return(end)
   }


# first and second derivatives of density estimate
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   D + D^2 dim vector: first D entries are gradient, last D^2 hessian

deriv.fun  = function(pt,dat,bw) {
   D       = ncol(dat)
   wi      = hatfi.fun(pt,dat,bw)
   fhatx   = mean(wi)
   msx     = colMeans(wi*dat)/(fhatx)
   grad    = (pt - msx)/(bw^2)
   crosspr = matrix(0, nrow=D, ncol=D)
   for(j in 1:D){
      crosspr[,j]= colMeans(wi*dat[,j]*dat)/(bw^4*fhatx)
      }
   hessian = crosspr - (msx %*% t(msx))/(bw^4) - diag(1/bw^2, nrow=D)
   
   return(c(grad, hessian))
   }


# spectral decomposition of hessian at one point
# input:
#   pt   one point in the space (D-dim vector)
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   eigenvalues of the hessian

spectral.fun= function(pt,dat,bw) {
   D       = ncol(dat)
   deriv   = deriv.fun(pt,dat,bw)
   hessian = matrix(deriv[-(1:D)], ncol=D)
   sort(eigen(hessian, symmetric=TRUE)$values, decreasing=TRUE)
}



# cluster detected modes (up to some digits)
# input:
#    modemat    matrix, each row is one detected mode
#    digits     number of digits 

modecl.fun = function(modemat,digits){
     X = round(modemat,digits)
     k = ncol(X)
     X = unique(X)
     return(matrix(X, ncol=ncol(modemat)))
     }

# # Hessian eigenvalues estimation
### THIS FUNCTION APPLIES TO A MATRIX
# input:
#   pts  matrix of critical points
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   matrix of eigenvalues

eigenest.fun   = function(pts,dat,bw){
   eigenval    = t(apply(X=pts, MARGIN=1, FUN=spectral.fun, dat, bw))
   return(matrix(eigenval, ncol=ncol(dat)))
}

# # gradient estimation
### THIS FUNCTION APPLIES TO A MATRIX
# input:
#   pts  matrix of critical points
#   dat  data points (nXD data matrix)
#   bw   bandwidth (scalar)
# output:
#   matrix of gradients

gradest.fun   = function(pts,dat,bw){
   grad       = t(apply(X=pts, MARGIN=1, FUN=deriv.fun, dat, bw)[1:ncol(dat),])
   return(matrix(grad, ncol=ncol(dat)))
}


# compute elementary symmetric polynomials
# input: vector of eigenvalues
# output: elementary simmetric polynomials

esp.fun = function(lam){
     out      = rep(0, length(lam))
     for (i in 1:length(lam)){
        subs  = matrix(combn(lam,i), nrow=i)
        out[i]= sum(apply(subs,2,prod))
        }
     return(out)
}



# roots of symmetric polynomial
# input: vector of elementary symmetric polynomials
# output: roots

roots.fun  = function(pol){
     signs = rep(c(-1,1), length=length(pol))
     coef  = c(rev(signs*pol),1)
     out   = Re(polyroot(coef))
     return(sort(out, decreasing=TRUE))
}


# bootstrap eigenvalues at one critical point
# input:
#   pt       critical point (D-dim vector)
#   dat      data points (nXD data matrix) to resample from
#   nboot    number of bootstrap replications 
#   nsample  size of bootstrap samples (default is nrow(dat))
#   bw       bandwidth (scalar)
# output:
#   eig      matrix with eigenvalues at candidate critical points
#            each row is one bootstrap replication

nonparboot.fun   = function(pt, dat, nboot, bw){
     eig         = matrix(data= 0, nrow=nboot, ncol=ncol(dat))
     for (i in 1:nboot){
        bootX    = sample(x=1:nrow(dat), size=nrow(dat), replace=TRUE)
        bootX    = matrix(dat[bootX,], ncol=ncol(dat)) #qui era ncol(X)
        eig[i,]  = spectral.fun(pt=pt,dat=bootX,bw=bw)
     }
     return(eig)
}


# find 1-alpha bootstrap confidence rectangle for all eigenvalues at one critical point
# input:
#   pt       critical point (D-dim vector)
#   dat      data points (nXD data matrix) to resample from
#   nboot    number of bootstrap replications 
#   bw       bandwidth (scalar)
#   alpha    for confidence interval
# output:    2-dim vector with extremes of CI of the leading eigenvalue

confint.fun = function(pt, dat, nboot, bw, alpha){
     shat   = esp.fun(spectral.fun(pt,dat,bw))
     booteig= nonparboot.fun(pt=pt, dat=dat, nboot=nboot, bw=bw)
     bootesp= t(apply(X=booteig, MARGIN=1, FUN=esp.fun))
     absdiff= abs(shat-t(bootesp))
     infdist= apply(X=absdiff, MARGIN=2, FUN=max)
     crit   = quantile(infdist, 1-alpha)
     keep   = (1:nboot)[infdist<= crit]
     #range(booteig[keep,1])
     out    = t(apply(X=matrix(booteig[keep,], ncol=ncol(dat)), MARGIN=2, FUN=range))
     return(out)
}


# perform testing procedure
# input:
#   datX     data points for mean shift (D-col matrix)
#   datY     data points for eigenvalues confid interv (D-col matrix)
#   bw       bandwidth (scalar)
#   nboot    number of bootstrap replications 
#   alpha    1-confidence level
#   digits   number of digits for modes fusion
# output:
#   modes    D-col matrix, each row is one mode detected from mean shift
#   CIlead   2-col matrix, each row is the confid interv relative of leading eigenv at corresponding mode
#   allCI    list, each entry correspond to a mode, the entry is a Dx2 matrix 
#            with conf interv of all eigenvalues at that mode

modetest.fun = function(datX, datY, bw, nboot=1000, alpha=0.05, digits=5){
     critical= matrix(t(apply(X=datX, MARGIN=1, FUN=msiter.fun, dat=datX, bw=bw)), ncol=ncol(datX))
     modes   = modecl.fun(modemat=critical, digits=digits)
     allCI   = list()
     CIlead  = matrix(0, nrow=nrow(modes), ncol=2)
     for (mode in 1:nrow(modes)){
          allCI[[mode]]= confint.fun(pt=modes[mode,], dat=datY, nboot=nboot, bw=bw, alpha=alpha/nrow(modes))
          CIlead[mode,]= allCI[[mode]][1,]
     }
     return(list(modes=modes,  CIlead=CIlead, allCI=allCI))
}


modetest.fun2 = function(datX, datY, bw, nboot=1000, alpha=0.05, out.ms){
  #critical= out.ms$end.points
  #modes   = modecl.fun(modemat=critical, digits=digits)
  modes <- matrix(out.ms$mode,ncol=ncol(datX))
  allCI   = list()
  CIlead  = matrix(0, nrow=nrow(modes), ncol=2)
  for (mode in 1:nrow(modes)){
    allCI[[mode]]= confint.fun(pt=modes[mode,], dat=datY, nboot=nboot, bw=bw, alpha=alpha/nrow(modes))
    CIlead[mode,]= allCI[[mode]][1,]
  }
  return(list(modes=modes,  CIlead=CIlead, allCI=allCI))
}

modetest.fun3 = function(datX, datY, bw, nboot=1000, alpha=0.05, modes){
  #critical= out.ms$end.points
  #modes   = modecl.fun(modemat=critical, digits=digits)
  #modes <- matrix(out.ms$mode,ncol=ncol(datX))
  allCI   = list()
  CIlead  = matrix(0, nrow=nrow(modes), ncol=2)
  for (mode in 1:nrow(modes)){
    allCI[[mode]]= confint.fun(pt=modes[mode,], dat=datY, nboot=nboot, bw=bw, alpha=alpha/nrow(modes))
    CIlead[mode,]= allCI[[mode]][1,]
  }
  return(list(modes=modes,  CIlead=CIlead, allCI=allCI))
}
# multibandtest
# input:
#   datX     data points for mean shift (D-col matrix)
#   datY     data points for eigenvalues confid interv (D-col matrix)
#   bwgrid   grid of possible bandwidths (vector)
#   nboot    number of bootstrap replications 
#   alpha    1-confidence level
#   digits   number of digits for modes fusion
# output:    list, each element is result of modetest.fun at corresponding bandwidth
   
multibandtest.fun = function(datX, datY, bwgrid, nboot=1000, alpha=0.05, digits=5){
     ngrid   = length(bwgrid)
     testres = list()
     for (i in 1:ngrid){
          testres[[i]] =modetest.fun(datX=datX, datY=datY, bw=bwgrid[i], nboot=nboot, alpha=alpha, digits=digits)
     }
     return(testres)
}



# plot confidence interval of leading eigenvalue at all candidate modes
# input:
#   rect     2-col matrix, each row is the confid interv at one mode
# output:    plot

CIplot.fun  = function(rect){
     bound  = max(abs(rect))
     k      = nrow(rect)
     plot(c(-bound, bound),c(0,k+1),type="n", xlab=expression(lambda[1]), ylab="Modes", yaxp= c(1,k,max(1,k-1)), las=1, )
     # plot(c(-bound, bound),c(0,k+1),type="n",yaxt="n",xlab=expression(lambda[1]),ylab="")
     # text(x=rep(-bound, k), y=0.2 + (1:k), labels=paste("M", 1:k, sep=""))
     for(i in 1:nrow(rect)){
          segments(rect[i,1],i,rect[i,2],i,lwd=5)
          segments(-bound,i,bound,i,lwd=.5)
         }
    abline(v=0,col="red")
}

# plot conf interv of all eigenvalues at all detected modes (even non signifiant)
# input:
#   allci  list of matrices (output allCI of modetest.fun)
#   isoscale TRUE if want the same scale on all plots, FALSE otherwise
# output   a series of plots

eigenport.fun= function(allci, isoscale=TRUE){
     bound = max(abs(unlist(allci)))
     for (mode in 1:length(allci)){
          rect =allci[[mode]]
          if (isoscale==FALSE) bound=max(abs(rect))
          k    = nrow(rect)
          #plot(c(-bound, bound),c(0,k+1),type="n", xlab=paste("Mode", mode), ylab=expression(lambda), 
             #  yaxp= c(1,k,max(1,k-1)), las=1, )
plot(c(-bound, bound),c(0,k+1),type="n",yaxp= c(1,k,max(1,k-1)), las=1 )
title(xlab=expression(lambda),ylab="",mgp=c(2,1,1))     

for(i in 1:nrow(rect)){
               segments(rect[i,1],i,rect[i,2],i,lwd=5)
               segments(-bound,i,bound,i,lwd=.5)
          }
         abline(v=0,col="red")
     }
}


# plot points and other stuff
# input:
#   dat     (nxD) matrix, each row is a point to plot
#   stuff   D-col matrix containing other points to be plotted (for instance detected modes)
#   isoscale TRUE if want the same scale on all dimensions, FALSE otherwise (needed only for D=2,3)
#   bw      bandwidth for density estimation (needed only for D=1)
# output:    plot

pointsplot.fun = function(dat, stuff, isoscale=TRUE, bw=0){
     D     = ncol(dat)
     if(D==1) {
          grid=seq(min(dat), max(dat), length=100)
          dens= sapply(X=grid, FUN=fhat.fun, dat=dat, bw=bw)
          plot(grid, dens, ylim=c(0, max(dens)*1.1), type="l", main="", xlab="", ylab="", yaxt="n")
          abline(v=stuff, col="red")
          text(x=stuff[,1], y=0, labels=1:nrow(stuff), col="blue")
     }
     if(D==2) {
          ran   = apply(X=dat, MARGIN=2, FUN=range)
          midran= (ran[2,] + ran[1,])/2
          maxran= ran[2,] - ran[1,]
          if(isoscale==TRUE) maxran= max(ran[2,] - ran[1,])
          plran = rbind(midran - 0.6*maxran, midran + 0.6*maxran)
          plot(dat, pch=19, xlim=plran[,1], ylim=plran[,2], xlab="", ylab="")
          points(stuff, pch=19, col="red")
          text(x=stuff[,1], y=stuff[,2], labels=1:nrow(stuff), col="yellow")

     }
     if(D==3) {
          ran   = apply(X=dat, MARGIN=2, FUN=range)
          midran= (ran[2,] + ran[1,])/2
          maxran= ran[2,] - ran[1,]
          if(isoscale==TRUE) maxran= max(ran[2,] - ran[1,])
          plran = rbind(midran - 0.6*maxran, midran + 0.6*maxran)
          pl1   = scatterplot3d(dat, pch=19,  xlim=plran[,1], ylim=plran[,2], zlim=plran[,3],
               xlab="",ylab="", zlab="") 
          pl1$points3d(stuff, pch=19, col="red")
          text(pl1$xyz.convert(stuff), labels=1:nrow(stuff), col="yellow")
     }
     if(D>3) print(paste(D,"-dimensional data - cannot plot", sep=""))
}

