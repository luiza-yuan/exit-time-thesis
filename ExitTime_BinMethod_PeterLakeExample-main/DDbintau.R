# Function for bin method to compute tau-corrected D1 and D2
# SRC 2021-02-16

# original program: Step3_Dd+bins+tau+wts+np_Squeal2_2021-02-14.R

DDbins = function(Xvar,bw,ntau,nbin) {

  nx = length(Xvar)
  taus = c(1:ntau)  # choose taus, usually 5 but can go as high as 15
  taumax = ntau

  # Make dx matrix
  X0all = Xvar[1:(nx-taumax)]
  nxs = length(X0all)
  dxmat = matrix(0,nr=nxs,nc=ntau)
  for(i in 1:ntau)  {
    tau = taus[i]
    X0 = Xvar[1:(nx-tau)]
    X1 = Xvar[(tau+1):nx]
    dx = X1-X0
    dxmat[,i] = t(dx[1:nxs])
  }

  # Make bins
  alow <- min(X0all)+bw
  ahigh <- max(X0all)-bw
  na = nbin
  bin.edges <- seq(alow,ahigh,length.out=(na+1)) #edges of the bins
  bin.mid = (bin.edges[1:na] + bin.edges[2:(na+1)])/2 # centers of the bins

  # make bin labels using cut
  Xdtf = as.data.frame(cbind(X0all,dxmat))
  Xdtf$bin <- cut(Xdtf$X0all, breaks = bin.edges, labels = 1:na)

  # moments by bin
  X01 = matrix(0,nr=na,nc=2)  # mean X0 for bins, and bin number
  M1 = matrix(0,nr=na,nc=ntau) # Matrix of M1: na rows, ntau columns
  M2 = M1 # Matrix for M2
  M4 = M1 # Matrix for M4
  # vector for regression correction
  xvec = as.vector(taus)
  # vectors to hold D1 and D1
  D1 = rep(0,na)
  D2 = rep(0,na)
  nbin = rep(0,na)
  # calculate moments, D1, D2 by bin
  for(ib in 1:na) {
    Xb = subset(Xdtf,subset = (Xdtf$bin==ib) )
    X01[ib,1] = ib
    X01[ib,2] = mean(Xb[,1]) # Mean of X0
    nxb = length(Xb[,1])
    nbin[ib] = nxb
    # Calculate moments by bin
    dx1 = Xb[,2:(ntau+1)]
    dx2 = dx1*dx1 # dx is also columns 2 to 1+ntau of Xdtf
    dx4 = dx2*dx2
    M1 = colMeans(dx1)
    M2 = colMeans(dx2)
    M4 = colMeans(dx4)
    # Regressions for D1 and D2
    yvec = M1
    # D1
    inv.erD1 = nxb/(M2 -(M1^2) )
    WD1 = diag(inv.erD1)
    D1[ib] = ( 1/(t(xvec)%*%WD1%*%xvec) )*( t(xvec)%*%WD1%*%yvec )
    D1[ib] = ifelse(nxb <= 1,D1[(ib-1)],D1[ib]) # adjust for bins with only 1 observation
    # D2
    inv.erD2 = nxb/(M4 -(M2^2) )
    WD2 = diag(inv.erD2)
    yvec = M2 - (D1[ib]*xvec)*(D1[ib]*xvec)
    xvec2 = 2*xvec
    D2[ib] = 0.5*( 1/(t(xvec2)%*%WD2%*%xvec2) )*(t(xvec2)%*%WD2%*%yvec)
    D2[ib] = ifelse(nxb <= 1 || D2[ib]<=0,D2[(ib-1)],D2[ib]) # adjust for bins with only 1 observation
  }

  # # KCE: Fix D2[1] = NA (if D2 < 0, previous value is taken, but for first value this results in NA)
  # if (any(is.na(D2) | is.na(D1))){
  #   D2[1] = D2[2]
  # }

  # calculate sigma
  sigma = sqrt(2*D2)

  # Nonparametric smoothing of D1 and D2 using ksmooth
  D1s = ksmooth(x=X01[,2],y=D1,kernel='normal',bandwidth=bw,
                x.points=bin.mid)
  D2s = ksmooth(x=X01[,2],y=D2,kernel='normal',bandwidth=bw,
                x.points=bin.mid)
  sigmas = ksmooth(x=X01[,2],y=sigma,kernel='normal',bandwidth=bw,
                   x.points=bin.mid)

  outlist = list(D1s,D2s,sigmas,bin.mid,D1,D2,sigma)
  # print('output list',quote=F)
  # print('list(D1s,D2s,sigmas,bin.mid,D1,D2,sigma)',quote=F)
  return(outlist)
}
