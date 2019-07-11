#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  7 15:09:38 2018

@author: Tianyan
"""

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects.packages import importr
from rpy2 import robjects
from rpy2.robjects.vectors import StrVector

def medicc_ce(D,i):
    if i==0:
        #only need to install once
        utils = importr("utils")
        utils.chooseCRANmirror(ind=1)
        packnames = ('kernlab', 'spatstat')
        utils.install_packages(StrVector(packnames))
    
    #need to import the packages
    kernlab = importr('kernlab')
    spatstat = importr('spatstat')
    
    print(i)
    
    #custom function for the MEDICCquant
    string = """

medicc.clonal.expansion=function(D, n=30){
  ## format distance matrix
  mydist=D
  if (any(is.na(mydist))) {
    warn("Matrix contained NA entries!")
    return(NA)
  }
  
  print(mydist)
  
  if(ncol(mydist) <= 1){
    return(NA)
  }
  
  Kdist = medicc.dist.to.kernmat(mydist)
  KPCA=kpca(Kdist)
  
  print(KPCA)
  
  avgdist = sum(mydist)/2
  nsamples =nrow(Kdist)
  
  if(ncol(rotated(KPCA)) < 2){
    return(NA)
  }
  
  x=rotated(KPCA)[,1]
  y=rotated(KPCA)[,2]
  
  print(x)
  print(y)
  
  if (nrow(Kdist)<=4) {
    f=2
  } else {
    f=1.3/sqrt(1-4/nrow(Kdist))
  }
  PPP=ppp(x,y,window=ripras(x,y,shape="rectangle",f=f))

  ## compute envelope, nsim=19 is 5% alpha level according to http://www.csiro.au/files/files/p10ib.pdf
  if (n>=1) {
    pb=txtProgressBar(min=0,max=n, style=3)
    lstat.replic=sapply(1:n,function(x) {
      setTxtProgressBar(pb,x)
      ENV=suppressWarnings(spatstat::envelope(PPP, Lest, nsim=19, nrank=1,correction="best", global=T, verbose=F))
      obstheo=abs(ENV$obs-ENV$theo) ## distance observed <-> theoretical
      hitheo=abs(ENV$hi-ENV$theo) ## width of confidence bound around theoretical
      maxidx=which.max(obstheo)
      lstat = obstheo[maxidx] / hitheo[maxidx]
      return(lstat)
    })
    lstat = mean(lstat.replic)
    L.stats.all =lstat
  } else {
    L.stats.all=NA
  }
  return(L.stats.all)
}

medicc.dist.to.kernmat = function(D) {
  mydist=D
  return(as.kernelMatrix(exp(-1*mydist/max(mydist))))
}

    """

    #call the MEDICCquant function func
    func = SignatureTranslatedAnonymousPackage(string, "func")
    
    #need to make the python matrix into an R matrix
    D_r = robjects.r.matrix(robjects.IntVector(mat_to_list(D)), nrow=len(D))
    
    return(func.medicc_clonal_expansion(D_r))

#converts a matrix to one long list
def mat_to_list(M):
    out = []
    for c in range(len(M[0])):
        for r in range(len(M)):
            out.append(M[r][c])
    return out