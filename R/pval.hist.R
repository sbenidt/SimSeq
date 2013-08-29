pval.hist <-
function(p)
  
{
  m<-length(p)                            # count the p-values 
  h.cdf<-h.pdf<-p.edf<-rep(NA,length(p))  # initialize h.pdf and p.edf vectors
  ord<-order(p)                           # order the p-values
  p.edf[ord]<-(1:m)/(m+1)                 # define the EDF estimator
  hh<-NULL                                # initialize the histogram heights 
  F0<-p0<-b.edf<-p.brks<-1                # initialize vectors with F0, p0, EDF break-points, and p-value break-points
  
  while(any(p<p0))  # Loop until all p-values are included in the histogram
  {
    slp<-min(1,(2*sum(p[p<p0])+1)/(1+p0*m))  # deterime the average abundance of p-values (under conditional uniform)
    slp.line<-(1-slp*(p0-p)/p0)*F0           # represent this line as a straight line in the EDF  
    ok<-(p<p0)&(p.edf>=slp.line)             # find the p-values that are "OK"
    if (any(ok))                             # Insert a new histogram break-point
    {
      p.cut<-min(p[ok])                   # Find the minimum p-value that is OK
      p.brks<-c(p.cut,p.brks)             # Update the p-value break-point vector
      hh<-c(slp,hh)                       # Update the histogram heights vector 
      p0<-p.cut                           # Update p0
      F0<-min(p.edf[p>=p0])               # Update F0
      b.edf<-c(F0,b.edf)                 # Update the break-points in the EDF histogram estimator
    }
    else # inserts a (p,EDF) point at (0,0)
    {
      p.brks<-c(0,min(p),max(p[p<p0]),p.brks) 
      b.edf<-c(0,min(p.edf[p==min(p)]),min(p.edf[p==max(p[p<p0])]),b.edf)
      p0<-0
    }
  }
  
  last.pt<-min((1:length(p.brks))[p.brks==1])
  p.brks<-p.brks[1:last.pt]
  b.edf<-b.edf[1:last.pt]
  
  p.diffs<-diff(p.brks)
  zero.diffs<-(1:length(p.diffs))[p.diffs==0]
  if (length(zero.diffs)>0)
  {
    p.brks<-p.brks[-(zero.diffs+1)]
    b.edf<-b.edf[-(zero.diffs+1)]   
  }
  
  hh<-exp(log(diff(b.edf))-log(diff(p.brks)))  
  if (any(is.na(hh)))
  {
    print(summary(p))
    print(p.brks)
    print(b.edf)
    print(hh)
    stop("Error: undefined histogram heights.")
    
  } 
  hd<-diff(hh)  # Difference in histogram heights
  p1<-1
  while(any(hd>0))  # Enforce monotonicity
  {
    k<-max((1:length(hd))[hd>0])                        # Find largest histogram bar index k such that bar k is shorter than bar k+1
    p.brks<-c(p.brks[1:k],p.brks[(k+2):length(p.brks)]) # combine histogram bars k and k+1
    b.edf<-c(b.edf[1:k],b.edf[(k+2):length(b.edf)])     # update EDF break-points with monotonocity constraint      
    hh<-exp(log(diff(b.edf))-log(diff(p.brks)))         # compute new histogram heights
    hd<-diff(hh)                                        # compare consecutive histogram heights
  }
  
  h.cdf<-approx(p.brks,b.edf,xout=p)$y                  # Get the histogram-based CDF estimates for each p-value
  p.hist<-exp(log(diff(b.edf))-log(diff(p.brks)))       # Get the height for each histogram bar
  pi0.hat<-min(p.hist)                                  # Get the pi0 estimate from the histogram
  h.ebp<-approx(p.brks,pi0.hat/c(p.hist,p.hist[length(p.hist)]),xout=p)$y # Get the conservative EBP interpolation
  h.fdr<-exp(log(min(pi0.hat))+log(p)-log(h.cdf))                         # Get the histogram-based FDR estimate
  h.ebp[p==0]=0
  h.fdr[p==0]=0
  
  return(list(p=p,              # input p-value vector
              edf=p.edf,        # the p-value edf
              h.cdf=h.cdf,      # the histogram-based cdf estimate
              h.fdr=h.fdr,      # the histogram-based FDR estimate
              h.ebp=h.ebp,      # the histogram-based EBP estimate
              p.brks=p.brks,    # the p-value break-points of the histogram
              p.hist=p.hist,    # the heights of each histogram bar
              edf.brks=b.edf,   # the break-points in the EDF of the histogram estimator
              pi0.hat=pi0.hat)) # the histogram-based estimate of the proportion of tests with a true null hypothesis
}
