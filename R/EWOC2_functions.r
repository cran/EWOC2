###############################################################################################################################
#                                                                                                                             #
#                       Internal Functions to be used by EWOC: 2 drugs combination                                            #
#                                                                                                                             #
###############################################################################################################################

# 1. Minimim Dose Increment Function

   Function.Mi.Dose.Increment_simu<-function(Minimum.Dose.Increment, Dose){
     if (Minimum.Dose.Increment==0) {
        rounded.dose<-Dose
     } else {
        int.seq<-seq(from=0, to=1, by=Minimum.Dose.Increment)
        rounded.dose<-int.seq[max(which(int.seq<=Dose))] # round down to a lower level
        if((Dose - rounded.dose) / Minimum.Dose.Increment > 0.5) rounded.dose = rounded.dose + Minimum.Dose.Increment # always round to the nearest level
     }

   # don't allow dose = 0
   if (rounded.dose==0 & Minimum.Dose.Increment!=0) rounded.dose<-Minimum.Dose.Increment

  return(rounded.dose)
}

################################################################################################################################

# 2. Generate JAGS (BUGS) file 
  Function.generate.bugs.file <- function() { 
  
  cat(
"model {

for (i in 1:N) {

Z[i] ~ dbern(p[i])
logit(p[i]) <- logit(rho00) + X[i]*(logit(rho10)-logit(rho00))+Y[i]*(logit(rho01)-logit(rho00))+eta*X[i]*Y[i]
}

rho01 ~ dbeta(a01,b01)
rho10 ~ dbeta(a10,b10)
temp ~ dbeta(a00,b00)
rho00 <- temp*min(rho01,rho10)
eta~dgamma(a,b)

mtdx1<-((logit(theta)-logit(rho00))-(logit(rho01)-logit(rho00))*Y[N-1])/ ((logit(rho10)-logit(rho00))+eta*Y[N-1])

mtdy1<-((logit(theta)-logit(rho00))-(logit(rho10)-logit(rho00))*X[N-1])/ ((logit(rho01)-logit(rho00))+eta*X[N-1])

mtdx2<-((logit(theta)-logit(rho00))-(logit(rho01)-logit(rho00))*Y[N])/ ((logit(rho10)-logit(rho00))+eta*Y[N])

mtdy2<-((logit(theta)-logit(rho00))-(logit(rho10)-logit(rho00))*X[N])/ ((logit(rho01)-logit(rho00))+eta*X[N])

}",
file = "combination_cohort_2.bug" , fill = TRUE, sep = "")

} # end of function

#################################################################################################################################

# 3. function to calculate the dose for next cohort of patients

  Function.nextdose.2d <- function(dose.a,dose.b,resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, a01,b01,a10,b10,a00,b00,a,b, MDI.A, MDI.B, delta1x, delta1y, burn,mm, delta1) {

   # define/create useful variables
   # parameters for priors

   # data
       Z = resp
       X = (dose.a - Min.Dose.A) / (Max.Dose.A - Min.Dose.A)
       Y = (dose.b - Min.Dose.B) / (Max.Dose.B - Min.Dose.B)
       n = length(Z)
       i = n/2 + 1

    # Minimum.Dose.Increment
       Minimum.Dose.Increment.X = MDI.A / (Max.Dose.A - Min.Dose.A)
       Minimum.Dose.Increment.Y = MDI.B / (Max.Dose.B - Min.Dose.B)

    # set up for lower bound
       lbA = - Min.Dose.A / (Max.Dose.A - Min.Dose.A)
       lbB = - Min.Dose.B / (Max.Dose.B - Min.Dose.B)

    # setup mcmc parameters
       chains = 1

    # call jags
       j=jags.model('combination_cohort_2.bug',data=list('Z'=Z,'X'=X,'Y'=Y,'theta'=theta,'a00'=a00,'b00'=b00,'a01'=a01,'b01'=b01,'a10'=a10,'b10'=b10,'a'=a,'b'=b,'N'=n),n.chains=chains,n.adapt=burn, quiet=TRUE)

       s=coda.samples(j,c('rho00','rho01','rho10','eta','mtdx1','mtdy1','mtdx2','mtdy2'),mm)
       ss=as.data.frame(s[[1]])

       # Calculating posterior probability of toxicity
       temp1 = numeric()
       temp2 = numeric()
       for (ll in 1:mm) {
       temp1[ll]<-pdlt(ss$rho00[ll],ss$rho01[ll],ss$rho10[ll],ss$eta[ll],theta, 0, 0)
       if (temp1[ll] > theta + delta1)
       temp2[ll]<-1 else
       temp2[ll]<-0
       }
       postlow <- sum(temp2)/mm

       trcmtdx1<-ss$mtdx1[ss$mtdx1 > lbA]
       xx1<-max(0,quantile(trcmtdx1,alpha))
       xx1<-min(xx1,1)
       if ((xx1 - X[2*i-3]) > delta1x)
       xx1<-X[2*i-3]+delta1x
       xx1 = Function.Mi.Dose.Increment_simu(Minimum.Dose.Increment.X, xx1)

       trcmtdx2<-ss$mtdx2[ss$mtdx2 > lbA]
       xx2<-max(0,quantile(trcmtdx2,alpha))
       xx2<-min(xx2,1)
       if ((xx2 - X[2*i-2]) > delta1x)
       xx2<-X[2*i-2]+delta1x
       xx2 = Function.Mi.Dose.Increment_simu(Minimum.Dose.Increment.X, xx2)

       trcmtdy1<-ss$mtdy1[ss$mtdy1 > lbB]
       yy1<-max(0,quantile(trcmtdy1,alpha))
       yy1<-min(yy1,1)
       if ((yy1 - Y[2*i-3]) > delta1y)
       yy1<-Y[2*i-3]+delta1y
       yy1 = Function.Mi.Dose.Increment_simu(Minimum.Dose.Increment.Y, yy1)

       trcmtdy2<-ss$mtdy2[ss$mtdy2 > lbB]
       yy2<-max(0,quantile(trcmtdy2,alpha))
       yy2<-min(yy2,1)
       yy2 = Function.Mi.Dose.Increment_simu(Minimum.Dose.Increment.Y, yy2)
       if ((yy2 - Y[2*i-2]) > delta1y)
       yy2<-Y[2*i-2]+delta1y
       yy2 = Function.Mi.Dose.Increment_simu(Minimum.Dose.Increment.Y, yy2)

       if(i==2) {
       nextX<-c(X[2*i-3], xx2)
       nextY<-c(yy1,Y[2*i-2])
       } else if (X[2*i-3] == X[2*i-5]) {
         nextX<-c(xx1,X[2*i-2])
         nextY<-c(Y[2*i-3],yy2)
       } else {
         nextX<-c(X[2*i-3],xx2)
         nextY<-c(yy1,Y[2*i-2])
       }

    # convert to original scale
      nextX = nextX * (Max.Dose.A - Min.Dose.A) + Min.Dose.A
      nextY = nextY * (Max.Dose.B - Min.Dose.B) + Min.Dose.B  

    # Return the next dose
      return(list(nextX=nextX, nextY=nextY, vrho00=median(ss$rho00), vrho01=median(ss$rho01), vrho10=median(ss$rho10), veta=median(ss$eta), postlow=postlow))

} # End of ewoc function

################################################################################################################################

# 3. count of numbers after decimal point

  Function.decimalplaces <- function(x) {
  ifelse ( (x %% 1) != 0, nchar(strsplit(sub('0+$', '', as.character(x)), ".", fixed=TRUE)[[1]][[2]]), 0)
  }

################################################################################################################################

# 4. Dose-Response Probability function

  pdlt = function(rho00,rho01,rho10,eta,theta,x,y) {

     # input parameter checking
     if(theta<=0  | theta>=1 )  stop("Error in theta:  theta should be between 0 and 1") 
     if(rho00 <= 0 | rho00 >= 1) stop("rho00 must between 0 and theta")
     if(rho01 <= 0 | rho01 >= 1) stop("rho01 must be 0-1")
     if(rho10 <= 0 | rho10 >= 1) stop("rho10 must be 0-1")
     if(eta <= 0) stop("eta must be positive")
     if(x < 0 | x > 1) stop("doses must be 0-1")
     if(y < 0 | y > 1) stop("doses must be 0-1")

     alpha<-log(rho00/(1-rho00))
     beta<-log(rho10/(1-rho10)) - log(rho00/(1-rho00))
     gamma<-log(rho01/(1-rho01)) - log(rho00/(1-rho00))
     p<-1/(1+exp(-alpha-beta*x-gamma*y-eta*x*y))
     return(p)

  } # End of Dose-Response Probability function

################################################################################################################################

# 5. Simulation function

  Function.simu.2d = function(type, alpha, theta, rho00,rho01,rho10,eta, data1, N, M, a01,b01,a10,b10,a00,b00,a,b, MDI.A, MDI.B, vai, delta1x,delta1y,burn,mm,delta1){

  # determine the initial dose and probability of DLT
    dose.a = c(0,0)
    dose.b = c(0,0)
    
    if(type %in% c("c","C","continuous")) {
    p1 = pdlt(rho00,rho01,rho10,eta,theta,dose.a[1],dose.b[1])
    p2 = pdlt(rho00,rho01,rho10,eta,theta,dose.a[2],dose.b[2])
    } else if(type %in% c("d", "D", "discrete")) {
    p1 = data1$p[data1$x == dose.a[1] & data1$y == dose.b[1]]
    p2 = data1$p[data1$x == dose.a[2] & data1$y == dose.b[2]]	
    }
    
    # simulate response
    resp1 = rbinom(1, 1, p1)
    resp2 = rbinom(1, 1, p2)
    resp = c(resp1, resp2)
    
    # Feasibility bound (alpha)
    alpha.vector = rep(alpha, N)

    # stopping rule
    postlow = numeric()
  # start simulation
    for (i in 1:(N/2)) {
    	
    nextdoses = Function.nextdose.2d(dose.a, dose.b ,resp, theta, alpha.vector[i], Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, a01,b01,a10,b10,a00,b00,a,b, MDI.A, MDI.B, delta1x,delta1y,burn,mm,delta1)
    	
    X = nextdoses$nextX
    Y = nextdoses$nextY
    postlow[i] = nextdoses$postlow	
    
    # determine true probability using either function (continuous) or true value (discrete)
    if(type %in% c("c","C","continuous")) {
    p1 = pdlt(rho00,rho01,rho10,eta,theta,X[1],Y[1])
    p2 = pdlt(rho00,rho01,rho10,eta,theta,X[2],Y[2])
    } else if(type %in% c("d", "D", "discrete")) {
    p1 = data1$p[data1$x == X[1] & data1$y == Y[1]]
    p2 = data1$p[data1$x == X[2] & data1$y == Y[2]]
	}
	
    resp1 = rbinom(1, 1, p1)
    resp2 = rbinom(1, 1, p2)
    
    # update data
    dose.a = c(dose.a,X)
    dose.b = c(dose.b,Y)
    resp = c(resp, resp1, resp2)
    
  # update alpha
    alpha.vector[i+1] = min(0.5, alpha.vector[i] + vai)
    }
    
    # output
    return(list(doseX=dose.a[1:N], doseY=dose.b[1:N], resp=resp[1:N], vrho00=nextdoses$vrho00, vrho01=nextdoses$vrho01, vrho10=nextdoses$vrho10, veta=nextdoses$veta, postlow=postlow))
   } # End of simulation function

################################################################################################################################
# 6. MTD curve finding function
      
     Function.mtd_logistic<-function(rho00,rho01,rho10,eta,theta,x){

     alpha<-log(rho00/(1-rho00))
     beta<-log(rho10/(1-rho10)) - log(rho00/(1-rho00))
     gamma<-log(rho01/(1-rho01)) - log(rho00/(1-rho00))
     a<-(log(theta/(1-theta))-alpha-beta*x)/(gamma+eta*x)
     a}

# End of MTD curve finding function
###############################################################################################################################
# 7. Distance to be minimized

     Function.fdist<-function(x,pt,rho00,rho01,rho10,eta,theta){
     # pt is the point on true MTD curve
     d<-(pt[1]-x)^2+(pt[2]- Function.mtd_logistic(rho00,rho01,rho10,eta,theta,x))^2
     d<-d^0.5
     d}
     
# End of distance function
###############################################################################################################################
# 8. Function to find MTD (discrete dose combination close to the curve)  

  Function.MTDdoses = function(rho00, rho01, rho10, eta, theta, nx, ny) {

   # input parameter checking
   if(theta<=0  | theta>=1 )  stop("Error in theta:  theta should be between 0 and 1") 
   if(rho00 <= 0 | rho00 >= 1) stop("rho00 must between 0 and 1")
   if(rho01 <= 0 | rho01 >= 1) stop("rho01 must be 0-1")
   if(rho10 <= 0 | rho10 >= 1) stop("rho10 must be 0-1")
   if(eta <= 0) stop("eta must be positive")
   if(nx %% 1 != 0 | nx <= 0) stop("number of dose levels for drug A must be positive integer")
   if(ny %% 1 != 0 | ny <= 0) stop("number of dose levels for drug B must be positive integer")

  # special senarios
  # do not recommend MTD if rho00 > theta + 0.1
  if(rho00 > theta + 0.1) {
  doses = data.frame(x=NA, y=NA)
  } else if(Function.mtd_logistic(rho00,rho01,rho10,eta,theta,0) < 0) {
  doses = data.frame(x=0, y=0)
  } else if(Function.mtd_logistic(rho00,rho01,rho10,eta,theta,1) > 1) {
  doses = data.frame(x=1, y=1)
  } else {

  # set number of descrete levels for drug X and Y
  x1 = seq(0,1, length.out=nx)
  y1 = seq(0,1, length.out=ny)

  # find intercept
  xmin = 0
  xmax = 1
  if(Function.mtd_logistic(rho00,rho01,rho10,eta,theta,0) > 1) xmin=uniroot(function(x) (log(theta/(1-theta))-log(rho00/(1-rho00))-(log(rho10/(1-rho10)) - log(rho00/(1-rho00)))*x)/(log(rho01/(1-rho01)) - log(rho00/(1-rho00))+eta*x)-1, c(0,1),extendInt="downX")$root
  if(Function.mtd_logistic(rho00,rho01,rho10,eta,theta,1) < 0) xmax=uniroot(function(x) (log(theta/(1-theta))-log(rho00/(1-rho00))-(log(rho10/(1-rho10)) - log(rho00/(1-rho00)))*x)/(log(rho01/(1-rho01)) - log(rho00/(1-rho00))+eta*x), c(0,1),extendInt="downX")$root

  ######
  # calculate the mimimum distance for each dose combinations
  ld = matrix(nrow=length(x1), ncol=length(y1))

  for (i in 1:length(x1))
  for (j in 1:length(y1)) 
  {
  {
  aal<-optimize(f=Function.fdist,c(xmin,xmax),tol=0.0001,c(x1[i],y1[j]),rho00,rho01,rho10,eta,theta)
  ld[i,j]<-aal$objective
  }
  }

 ######
  # find the best combinations

  # 1. select x fist
  ld1 = ld
  for (i in 1:length(x1)) ld1[i,which(ld1[i,] != min(ld1[i,],na.rm=TRUE))] = NA
  for (j in 1:length(y1)) {
  if (sum(!is.na(ld1[,j])) > 0)
  ld1[which(ld1[,j] != min(ld1[,j],na.rm=TRUE)),j] = NA
  }

  # 2. select y fist
  ld2 = ld
  for (j in 1:length(y1)) ld2[which(ld2[,j] != min(ld2[,j],na.rm=TRUE)),j] = NA
  for (i in 1:length(x1)) {
  if (sum(!is.na(ld2[i,])) > 0)
  ld2[i,which(ld2[i,] != min(ld2[i,],na.rm=TRUE))] = NA
  }

  # intersection ld1 and ld2
  temp = data.frame(ld1 = as.numeric(ld1), ld2 = as.numeric(ld2))
  ld = matrix(apply(temp[, 1:2], 1, mean), nrow=nx, ncol=ny) # intersection
  ld[is.nan(ld)]=  NA

  # determine the dose combinations for MTD
  dosex = matrix(x1,nrow=length(x1), ncol=length(y1))
  dosey = matrix(y1,byrow=TRUE,nrow=length(x1), ncol=length(y1))
  doses = data.frame(x = dosex[which(!is.na(ld))], y=dosey[which(!is.na(ld))])
  #
  #doses is the dose combinations that are close to the true MTD curve
  return(doses)
  }
  }

################################################################################################################################
# 9 Function to compute the posterior probability of DLT at dose combination (x,y) using logistic model
  Function.postdlt = function(NN, theta, nx, ny, triali, dosex, dosey, dlt, priors) {

  # Prior parameters for rho00, rho01, rho10, eta
  a00 <- priors$a00
  b00 <-priors$b00
  a01 <-priors$a01
  b01 <-priors$b01
  a10 <-priors$a10
  b10 <-priors$b10
  a <-priors$a
  b <-priors$b

  # Data for a full trial consists of NN X doses, Y doses, and DLT Z

  X<-dosex[triali,]
  Y<-dosey[triali,]
  Z<-dlt[triali,]

  #Compute the posterior distribution of model parameters give data 
  # mcmc parameters for trial conduct
  chains<-1
  burn<-10000
  mm<-5000

j=jags.model('combination_cohort_2.bug',data=list('Z'=Z,'X'=X,'Y'=Y,'theta'=theta,'a00'=a00,'b00'=b00,'a01'=a01,'b01'=b01,'a10'=a10,'b10'=b10,'a'=a,'b'=b,'N'=NN),n.chains=chains,n.adapt=burn, quiet=TRUE)

  s=coda.samples(j,c('rho00','rho01','rho10','eta','mtdx1','mtdy1','mtdx2','mtdy2'),mm)
  ss=as.data.frame(s[[1]])

  # Calculating Posterior Probability of DLT at MTD (xx,yy)
  temp1 = rep(NA, mm)
  postdlts = data.frame(xx=1:(nx*ny), yy=NA, postdlt=NA, trial=triali)

  k=0
  for (xx in seq(0,1,length.out=nx)) {
  for (yy in seq(0,1,length.out=ny)) {
  k=k+1
  for (ll in 1:mm) temp1[ll]<-pdlt(ss$rho00[ll],ss$rho01[ll],ss$rho10[ll],ss$eta[ll],theta, xx, yy)
  postdlts[k,1:3] = c(xx, yy, mean(temp1 > theta+0.1  | temp1 < theta-0.1))
  }
  }

  # return value
  return(postdlts)
  }
  
################################################################################################################################
# 10. Function to find true parameters for discrete scenario

   Function.parameter = function(data1) {
   	
   	rho00 = data1$p[data1$x==0 & data1$y==0]
   	rho01 = data1$p[data1$x==0 & data1$y==1]
   	rho10 = data1$p[data1$x==1 & data1$y==0]
   	rho11 = data1$p[data1$x==1 & data1$y==1]
   	alpha<-log(rho00/(1-rho00))
    beta<-log(rho10/(1-rho10)) - log(rho00/(1-rho00))
    gamma<-log(rho01/(1-rho01)) - log(rho00/(1-rho00))
    eta = log(rho11/(1-rho11)) -alpha - beta - gamma
    if(eta <=0) eta = 0.01
    return(list(rho00=rho00, rho01=rho01, rho10=rho10, eta=eta))
   }

################################################################################################################################
#                                                                                                                              # 
#                                                                                                                              #
#                                              EWOC R package functions                                                        #
#                                                                                                                              # 
#                                                                                                                              #
################################################################################################################################

# 1. ewoc main function

    ewoc2 <- function(dose.a, dose.b, resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, a01,b01,a10,b10,a00,b00,a,b,delta1x,delta1y,burn,mm,delta1) UseMethod("ewoc2")
    
    ewoc2.default <- function(dose.a, dose.b, resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, a01,b01,a10,b10,a00,b00,a,b, delta1x,delta1y, burn=4000,mm=2000, delta1=0.05) {

    # input parameter checking
      if(length(dose.a) != length(dose.b))  stop("Error in input data: incompatible dimensions for past doses for drug A and drug B") 
      if(length(dose.a) != length(resp))  stop("Error in input data: incompatible dimensions for dose and response") 
      if(length(resp) %% 2 !=0) stop("Error in input data: must have two patients in each cohort")
      if(min(dose.a) < Min.Dose.A)  stop("Error in input data: input dose out of range (min.dose, max.dose)") 
      if(min(dose.b) < Min.Dose.B)  stop("Error in input data: input dose out of range (min.dose, max.dose)") 
      if(max(dose.a) > Max.Dose.A)  stop("Error in input data: input dose out of range (min.dose, max.dose)")
      if(max(dose.b) > Max.Dose.B)  stop("Error in input data: input dose out of range (min.dose, max.dose)")
      if(sum(!as.factor(resp) %in% c("0","1")) > 0)  stop("Error in input data: response should be either 0 or 1")
      if(alpha<=0 | alpha>=1)  stop("Error in alpha: alpha should be between 0 and 1") 
      if(theta<=0  | theta>=1 )  stop("Error in theta:  theta should be between 0 and 1") 
      if(Min.Dose.A<0 | Min.Dose.B<0) stop("Error in min.dose: min.dose should >=0")
      if(Max.Dose.B<0 | Max.Dose.B<0) stop("Error in max.dose: max.dose should >=0")
      if(Max.Dose.A<=Min.Dose.A | Max.Dose.B<=Min.Dose.B) stop("Error in max.dose: max.dose should be bigger than min.dose")
      if(a<0 | b<0 | a01<0 | b01<0 | a10<0 | b10<0 | a00<0 | b00<0) stop("the beta parameters has to be positive")
      if(burn %% 1 != 0 | burn <= 0) stop("number of iterations for adaption must be positive integer")
      if(mm %% 1 != 0 | mm <= 0) stop("number of iterations to monitor must be positive integer")
      if(delta1<=0 | delta1>=1)  stop("Error in delta1: delta1 should be between 0 and 1")

    # delta1x and delta1y
      if(missing(delta1x)) delta1x = 0.2*(Max.Dose.A - Min.Dose.A)   
      if(missing(delta1y)) delta1y = 0.2*(Max.Dose.B - Min.Dose.B)   
      if(delta1x<=0 | delta1x>= (Max.Dose.A-Min.Dose.A))  stop("Error in delta1x: delta1x is out of range")
      if(delta1y<=0 | delta1y>= (Max.Dose.B-Min.Dose.B))  stop("Error in delta1y: delta1y is out of range")
   
    # call function to generate JAGs file
      Function.generate.bugs.file()
       
    # call nextdose function

      dosenext.rounded = Function.nextdose.2d(dose.a, dose.b, resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, a01,b01,a10,b10,a00,b00,a,b, MDI.A=0, MDI.B=0, delta1x/(Max.Dose.A - Min.Dose.A), delta1y/(Max.Dose.B - Min.Dose.B), burn,mm, delta1)
         
    # store the parameters and the results for future uses

      result = list(
      data = data.frame(dose.a = dose.a, dose.b=dose.b, response=resp),
      parameters = list(alpha=alpha, theta=theta, min.dose.a=Min.Dose.A, max.dose.a=Max.Dose.A, min.dose.b=Min.Dose.B, max.dose.b=Max.Dose.B, delta1x=delta1x, delta1y=delta1y, burn=burn, mm=mm, delta1=delta1),
      priors = list(a01=a01, b01=b01, a10=a10, b10=b10, a00=a00, b00=b00, a=a, b=b),
      nextdose.x = dosenext.rounded$nextX,
      nextdose.y = dosenext.rounded$nextY
      )

    # return
      class(result) <- "ewoc2"
      return(result)

      }
################################################################################################################################

# ewoc functions:

    # 1.1 print function for ewoc2

      print.ewoc2 <- function(x, ...){
      	
      # input parameter checking
      if(class(x) != "ewoc2") stop("wrong class, must be ewoc2")	
  
      round.length = max(Function.decimalplaces(x$parameters$min.dose.a), Function.decimalplaces(x$parameters$max.dose.a))
      cohort = nrow(x$data) / 2 + 1
      result = data.frame(Cohort=cohort, Patient=c(2*cohort-1, 2*cohort), Dose.A=NA, Dose.B=NA)
      result$Dose.A = format(x$nextdose.x, digits=round.length+2, nsmall=round.length+2)
      result$Dose.B = format(x$nextdose.y, digits=round.length+2, nsmall=round.length+2)
      print("The recommended doses for next cohort of patients are:")
      print(result)
      }  
    # End of print.ewoc function

#################################################################################################################################
# test run

# library(rjags)
#test = ewoc2(dose.a=c(0,0),dose.b=c(0,0),resp=c(0,0),theta=0.33,alpha=0.25, Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1,a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)
#print(test)

##################################################################################################################################

# 2. Simulation

    ewoc2simu <- function(ntrials,nsamples, type, trho00,trho01,trho10,teta, nx,ny,tp, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, alpha, theta, vai, a01,b01,a10,b10,a00,b00,a,b, delta1x,delta1y,burn,mm,delta1, seed) UseMethod("ewoc2simu")

      ewoc2simu.default <- function(ntrials, nsamples, type, trho00,trho01, trho10,teta, nx,ny,tp,  Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, alpha, theta, vai=0, a01,b01,a10,b10,a00,b00,a,b, delta1x,delta1y, burn=4000,mm=2000, delta1=0.05, seed=0){

    # input parameter checking
      if(ntrials %% 1 != 0 | ntrials <= 0) stop("number of trials must be positive integer")
      if(nsamples > 100)  stop("sample size must between 1 and 100 in phase I clinical trial")
      if(!type %in% c("continuous","discrete","c","d","C","D")) stop("type must be continous or discrete") 
      if(alpha<=0 | alpha>=1)  stop("Error in alpha: alpha should be between 0 and 1") 
      if(theta<=0  | theta>=1 )  stop("Error in theta:  theta should be between 0 and 1") 
      if(Min.Dose.A<0 | Min.Dose.B<0) stop("Error in min.dose: min.dose should >=0")
      if(Max.Dose.B<0 | Max.Dose.B<0) stop("Error in max.dose: max.dose should >=0")
      if(Max.Dose.A<=Min.Dose.A | Max.Dose.B<=Min.Dose.B) stop("Error in max.dose: max.dose should be bigger than min.dose")
      if(a<0 | b<0 | a01<0 | b01<0 | a10<0 | b10<0 | a00<0 | b00<0) stop("the beta parameters has to be positive")
      if(burn %% 1 != 0 | burn <= 0) stop("number of iterations for adaption must be positive integer")
      if(mm %% 1 != 0 | mm <= 0) stop("number of iterations to monitor must be positive integer")
      if(delta1<=0 | delta1>=1)  stop("Error in delta1: delta1 should be between 0 and 1")
   
    # delta1x and delta1y
      if(missing(delta1x)) delta1x = 0.2*(Max.Dose.A - Min.Dose.A)   
      if(missing(delta1y)) delta1y = 0.2*(Max.Dose.B - Min.Dose.B)   
      if(delta1x<=0 | delta1x>= (Max.Dose.A-Min.Dose.A))  stop("Error in delta1x: delta1x is out of range")
      if(delta1y<=0 | delta1y>= (Max.Dose.B-Min.Dose.B))  stop("Error in delta1y: delta1y is out of range")
   
    # continous scenario
      if(type %in% c("continuous","c","C")) {
      if(trho00 <= 0 | trho00 >= theta) stop("rho00 must between 0 and theta")
      if(trho01 <= 0 | trho10 <= 0 | teta <= 0) stop("rho01, rho10 and eta must be positive")
      if(trho01 >= 1 | trho10 >= 1) stop("rho01 and rho10 must be 0-1")
      MDI.A = 0
      MDI.B = 0
      }
   
    # discrete scenario
      if(type %in% c("discrete","d","D")) {
      if(nx %% 1 != 0 | nx <= 0) stop("number of dose levels for drug A must be positive integer")
      if(ny %% 1 != 0 | ny <= 0) stop("number of dose levels for drug B must be positive integer")
      if(min(tp)<=0 | max(tp) > 1) stop("true probabilities must between 0 and 1")
      if(length(tp) != nx*ny) stop("dimention of true probability doesn't match nx and ny")
      data1 = data.frame(x=rep(seq(0,1,length.out=nx),each=ny),y=rep(seq(0,1,length.out=ny),nx),p=tp)
      MDI.A = 1 / (nx-1)
      MDI.B = 1 / (ny-1)
      }
     
    # call function to generate JAGs file
      Function.generate.bugs.file()

    # call simulation function
      set.seed(seed)
      simresults = unlist(lapply(1:ntrials, function(x) Function.simu.2d(type, alpha, theta, trho00,trho01,trho10,teta, data1, N=nsamples, M=ntrials, a01,b01,a10,b10,a00,b00,a,b, MDI.A, MDI.B, vai, delta1x/(Max.Dose.A-Min.Dose.A),delta1y/(Max.Dose.B-Min.Dose.B), burn,mm, delta1)))

    # obtain results
      dose.A = matrix(simresults[substr(names(simresults),1,5) == "doseX"]*(Max.Dose.A-Min.Dose.A) + Min.Dose.A, ncol=nsamples, byrow=TRUE)
      dose.B = matrix(simresults[substr(names(simresults),1,5) == "doseY"]*(Max.Dose.B-Min.Dose.B) + Min.Dose.B, ncol=nsamples, byrow=TRUE)
      resp = matrix(simresults[substr(names(simresults),1,4) == "resp"],  ncol=nsamples, byrow=TRUE)
      postlow = matrix(simresults[substr(names(simresults),1,7) == "postlow"],  ncol=nsamples/2, byrow=TRUE)
      colnames(dose.A) = paste("patient", 1:nsamples)
      colnames(dose.B) = paste("patient", 1:nsamples)
      colnames(resp) = paste("patient", 1:nsamples)
      rownames(dose.A) = paste("trial", 1:ntrials)
      rownames(dose.B) = paste("trial", 1:ntrials)
      rownames(resp) = paste("trial", 1:ntrials)
      rownames(postlow) = paste("trial", 1:ntrials)

      rho00 = simresults[names(simresults)=="vrho00"]
      rho01 = simresults[names(simresults)=="vrho01"]
      rho10 = simresults[names(simresults)=="vrho10"]
      eta   = simresults[names(simresults)=="veta"]

      if(type %in% c("continuous","c","C")) {
      res = list(
      type = "continuous",
      parameters = list(ntrials=ntrials, nsamples=nsamples, trho00=trho00, trho01=trho01, trho10=trho10, teta=teta, alpha=alpha, theta=theta, Min.Dose.A=Min.Dose.A, Max.Dose.A=Max.Dose.A, Min.Dose.B=Min.Dose.B, Max.Dose.B=Max.Dose.B, alpha.increment=vai, delta1x=delta1x, delta1y=delta1y, burn=burn, mm=mm, delta1=delta1),
      priors = list(a01=a01,b01=b01,a10=a10,b10=b10,a00=a00,b00=b00,a=a,b=b),
      Dose.A=dose.A, Dose.B=dose.B, Resp=resp, rho00=rho00, rho01=rho01, rho10=rho10, eta=eta, postlow=postlow)
     
      } else if(type %in% c("discrete","d","D")) {
      
       # get postdlt
	   postdlts = data.frame(xx=NA, yy=NA, postdlt=NA, trial=rep(NA, nx*ny*ntrials))
	   priors = list(a01=a01,b01=b01,a10=a10,b10=b10,a00=a00,b00=b00,a=a,b=b)
	   for (jj in 1:ntrials) {
	   postdlts[(nx*ny*(jj-1)+1):(nx*ny*jj), ]= Function.postdlt(nsamples, theta, nx, ny, jj, dose.A, dose.B, resp, priors)
	   }
			
      res = list(
      type = "discrete",
      parameters = list(ntrials=ntrials, nsamples=nsamples, nx=nx, ny=ny, scenario=data1, alpha=alpha, theta=theta, Min.Dose.A=Min.Dose.A, Max.Dose.A=Max.Dose.A, Min.Dose.B=Min.Dose.B, Max.Dose.B=Max.Dose.B, alpha.increment=vai, delta1x=delta1x, delta1y=delta1y, burn=burn, mm=mm, delta1=delta1),
      priors = list(a01=a01,b01=b01,a10=a10,b10=b10,a00=a00,b00=b00,a=a,b=b),
      Dose.A=dose.A, Dose.B=dose.B, Resp=resp, rho00=rho00, rho01=rho01, rho10=rho10, eta=eta, postlow=postlow, postdlts=postdlts)	
      }
      
      # return
      class(res) <- "ewoc2simu"
      return(res)
    } # End of simulation function

################################################################################################################################
# test run
# b = ewoc2simu(ntrials=5, nsamples=8, type="c", trho00=0.01,trho01=0.2, trho10=0.9,teta=20, Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)
# b = ewoc2simu(ntrials=5, nsamples=8, type="d", nx=6, ny=3, tp=c(0.03,0.05,0.08,0.05,0.08,0.13,0.08,0.13,0.2,0.13,0.2,0.29,0.2,0.29,0.4,0.29,0.4,0.53),  Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)

################################################################################################################################
# ewoc functions for simulation:
   # 2.1 plot (efficiency)

     plot.ewoc2simu <- function(x, type="MTD", conf.reg=0.90, plot.figure="Y", ...){
  
   # input parameter checking
     if(class(x) != "ewoc2simu") stop("wrong class, must be ewoc2simu")
     if(!type %in% c("MTD","bias","percent")) stop("wrong type, must be MTD, bias or percent")
  
   # import the mcmc simulation results
     vrho00 <- x$rho00
     vrho01 <- x$rho01
     vrho10 <- x$rho10
     veta   <- x$eta

     dlt   <- x$Resp
     dosex <- x$Dose.A      
     dosey <- x$Dose.B
     postlow = x$postlow

   # read parameters
     M  = x$parameters$ntrials
     NN = x$parameters$nsamples
     theta = x$parameters$theta

     Min.Dose.A = x$parameters$Min.Dose.A
     Max.Dose.A = x$parameters$Max.Dose.A
     Min.Dose.B = x$parameters$Min.Dose.B
     Max.Dose.B = x$parameters$Max.Dose.B
     round.a = max(Function.decimalplaces(Min.Dose.A), Function.decimalplaces(Max.Dose.A))
     round.b = max(Function.decimalplaces(Min.Dose.B), Function.decimalplaces(Max.Dose.B))

     if(x$type == "continuous") {
     trho00 = x$parameters$trho00
     trho01 = x$parameters$trho01
     trho10 = x$parameters$trho10
     teta   = x$parameters$teta
     } else if(x$type == "discrete") {
     nx = x$parameters$nx
     ny = x$parameters$ny
     data1 = x$parameters$scenario	
     }

     # 1.1 plot TMD
# This program plots the true and estimated MTD of EWOC using drug combinations along with the last trial doses.
     if(type == "MTD") {       
       if(x$type == "continuous") { 
       tmtdx0<-(log(theta/(1-theta)) - log(trho00/(1-trho00))) / (log(trho10/(1-trho10)) - log(trho00/(1-trho00)))

       xtemp1<-seq(0,1,by=0.01)
       xtemp2<-seq(0,1,by=0.01)

       if (tmtdx0 > 1) {
       x2<-seq(0,1,by=0.01)} else {
	 x2<-seq(0,tmtdx0,by=0.01)}

	   ll<-length(x2)
	   y2<-numeric()
	   yest<-numeric()
	   for (i in 1:ll){
	   y2[i]<-Function.mtd_logistic(trho00,trho01,trho10,teta,theta,x2[i])}
	   for (i in 1:101){
	   yest[i]<-Function.mtd_logistic(mean(vrho00),mean(vrho01),mean(vrho10),mean(veta),theta,xtemp1[i])
	   }

	   xtemp2<-x2[y2 >= 0 & y2 <= 1]
         y2<-y2[y2 >= 0 & y2 <= 1]
         xtemp1<-xtemp1[yest >= 0 & yest <= 1]
	   yest<-yest[yest >= 0 & yest <= 1]

	   x3<-seq(0.02,1,by=0.01)
	   y3<-rep(0,99)
	   
	   # Graphing the True and Estimated MTD
	   par(pty="s")

	   plot(xtemp2,y2,type="l",xlim=c(0,1),ylim=c(0,1),xaxt="n", yaxt="n", xlab="Dose (Agent A)", ylab="Dose (Agent B)",main="True and Estimated MTD curve", cex.lab=0.8,cex.axis=0.8)
	   # adjust axes unit
         axis(1, at=axTicks(1), labels=format(axTicks(1)*(Max.Dose.A-Min.Dose.A) + Min.Dose.A,digits=round.a, nsmall=round.a))
         axis(2, at=axTicks(2), labels=format(axTicks(2)*(Max.Dose.B-Min.Dose.B) + Min.Dose.B,digits=round.b, nsmall=round.b))
         points(dosex[,NN],dosey[,NN],pch=1,col="grey",xlab="",ylab="")
	   points(dosex[,NN-1],dosey[,NN-1],pch=1,col="grey",xlab="",ylab="")
	   points(xtemp2,y2,type="l",xlab="",ylab="")
	   points(xtemp1,yest,type="l",lty=2,xlab="",ylab="")
	   legend("topleft",c("True MTD","Estimated MTD"),bty="n",lty=1:2,cex=0.8)

	   xx<-c(dosex[,NN],dosex[,NN-1])
	   yy<-c(dosey[,NN],dosey[,NN-1])
	   dd<-kde2d(xx,yy,n=2*M)
	   zz<-array()
	   for (i in 1:(2*M)){
	   z.x<-max(which(dd$x <= xx[i]))
	   z.y<-max(which(dd$y <= yy[i]))
	   zz[i]<-dd$z[z.x,z.y]
	   }
	   conf.border90<-quantile(zz,probs=1-conf.reg,na.rm=TRUE)
	   contour(dd,levels=conf.border90,labels=conf.reg,lty=1,add=TRUE,labcex=0.7,vfont=c("sans serif","bold"),col="red")
	   
       } else if(x$type == "discrete") {
       par(cex.axis=1.1, cex.lab=1.1, cex.main=1.1, las=1)
       plot(data1$x, data1$y, type="n",xlim=c(0,1),ylim=c(0,1),xlab="Dose (Drug A)",ylab="Dose (Drug B)",xaxt="n",yaxt="n",main=" Dose limiting toxicity scenario")
       abline(h=seq(0,1,length.out=ny),col="lightgray",lty="dotted")
       abline(v=seq(0,1,length.out=nx),col="lightgray",lty="dotted")
       axis(1, seq(0,1,length.out=nx),1:nx)
       axis(2, seq(0,1,length.out=ny),1:ny)
       text(data1$x, data1$y, format(round(data1$p,2),digits=2,nsmall=2))
       #data2 = subset(data1, p==theta)
       data2 = data1[data1$p == theta,]
       text(data2$x, data2$y, format(round(data2$p,2),digits=2,nsmall=2), col=2, font=2)
       }
       
     } else if(type == "bias") {
       if(x$type == "discrete") stop("average bias plot is not applicable for discrete scenario")
       # This program determines the design operating characteristics of EWOC using drug combinations.
       # For a given scenario determined by trho00, trho01, trho10, teta, theta are parameters of the scenario
       # For every estimated MTD from a given trial, and for each point x on the true MTD curve,
       # we cumpute the minimum distance dx from this point to the estimated MTD.
       # We then plot the average distance dx across the m simmulated trials with pointwise 95% confidence intervals. 

      	tmtdx0<-(log(theta/(1-theta)) - log(trho00/(1-trho00))) / (log(trho10/(1-trho10)) - log(trho00/(1-trho00)))
        x1<-seq(0,min(tmtdx0,1),by=0.01)

        kk<-length(x1)
        ymtd<-numeric()
        for(i in 1:kk){
        ymtd[i]<-Function.mtd_logistic(trho00,trho01,trho10,teta,theta,x1[i])}

        x1<-x1[ymtd >= 0 & ymtd <= 1]
        ymtd<-ymtd[ymtd >= 0 & ymtd <= 1]
        kk<-length(x1)

        ld<-matrix(nrow=M,ncol=kk)
        for (i in 1:M){
        for (j in 1 :kk){

        aal<- optimize(f=Function.fdist,c(0,1),tol=0.0001,c(x1[j],ymtd[j]),vrho00[i],vrho01[i],vrho10[i],veta[i],theta)
        ld[i,j]<- aal$objective
        if (ymtd[j] > Function.mtd_logistic(vrho00[i],vrho01[i],vrho10[i],veta[i],theta,aal$minimum)) ld[i,j] <- -ld[i,j]

        }}

        aveld<-numeric()

        for (i in 1 :kk){
        aveld[i]<-mean(ld[,i])
        }

        low<-min(aveld)
        upp<-max(aveld)

        # plot
        plot(x1,aveld,type="l",ylim=c(low,upp),xaxt="n", xlab="Dose (Drug A)", ylab="Average Bias",main="Pointwise average relative minimum distance", cex.lab=0.9,cex.axis=0.8)
        # adjust axes unit
        axis(1, at=axTicks(1), labels=format(axTicks(1)*(Max.Dose.A-Min.Dose.A) + Min.Dose.A,digits=round.a, nsmall=round.a))
        lines(x1,rep(0,kk),lty=2)
     	
     } else if(type == "percent") {
     	if(x$type == "continuous") { 
        tmtdx0<-(log(theta/(1-theta)) - log(trho00/(1-trho00))) / (log(trho10/(1-trho10)) - log(trho00/(1-trho00)))
        x1<-seq(0,min(tmtdx0,1),by=0.01)

        kk<-length(x1)
        ymtd<-numeric()
        for(i in 1:kk){
        ymtd[i]<-Function.mtd_logistic(trho00,trho01,trho10,teta,theta,x1[i])}

        x1<-x1[ymtd >= 0 & ymtd <= 1]
        ymtd<-ymtd[ymtd >= 0 & ymtd <= 1]
        kk<-length(x1)

        ld<-matrix(nrow=M,ncol=kk)
        for (i in 1:M){
        for (j in 1 :kk){

        aal<- optimize(f=Function.fdist,c(0,1),tol=0.0001,c(x1[j],ymtd[j]),vrho00[i],vrho01[i],vrho10[i],veta[i],theta)
        ld[i,j]<- aal$objective

        }}

        # caclulate the percentage
        lpercent1<-numeric()
        lpercent2<-numeric()

        for (i in 1 :kk){
        lpercent1[i]<-(length(ld[,i][ld[,i] <= 0.1*((x1[i]^2+ymtd[i]^2)^0.5)])/M)*100
        lpercent2[i]<-(length(ld[,i][ld[,i] <= 0.2*((x1[i]^2+ymtd[i]^2)^0.5)])/M)*100
        }

        # plot
        plot(x1,lpercent1,type="l",ylim=c(0,100),xaxt="n", xlab="Dose (Agent A)", ylab="Percent",main="Pointwise percent of MTD recommendation", cex.lab=0.9,cex.axis=0.8, col=1)
        # adjust axes unit
        axis(1, at=axTicks(1), labels=format(axTicks(1)*(Max.Dose.A-Min.Dose.A) + Min.Dose.A,digits=round.a, nsmall=round.a))
        lines(x1,lpercent2, lty=2, col=2)
        legend("bottomright",c("p=0.1","p=0.2"),lty=c(1,2),col=1:2, bty="n",cex=0.8)
        
        } else if(x$type == "discrete") {
        
        # parameters (delta)	
		delta = 0.8
		delta5 = 0.35
		
		# get postdlt
		postdlts = x$postdlts
		
		# recommend dose combination based on estiamted MTD curve

		result = data.frame(x=NA,y=NA,trial=NA)
		for (i in 1:M) {

		a = Function.MTDdoses(vrho00[i], vrho01[i], vrho10[i], veta[i], theta, nx, ny)

		#################
		# postdlt control, added 03/06/2015
		#temp1 = subset(postdlts, trial==i)
            temp1 = postdlts[which(postdlts$trial == i),]
		temp1$id = interaction(round(temp1$xx,4), round(temp1$yy,4))
		a$id = interaction(round(a$x,4), round(a$y,4))
		b=merge(a, temp1, by="id", all.x=TRUE)

		# If postdlt > delta5, exclude those dose combinations
		#b=subset(b, postdlt < 1-delta5, select=c(x,y))
            b = b[which(b$postdlt < 1-delta5),c("x","y")]
		a = b
		################

		if(nrow(a)>0) a$trial = i

		# do not recommend MTD if any postlow >= 0.8
		if (max(postlow[i,]) < delta) 
		result = rbind(result, a)
		}

		# remove NAs (including NA in the first row)
		result = result[!is.na(result$x),]

		# summary
		result2 = data.frame(x = 1:(nx*ny), y=NA, percent=NA)
		x1 = seq(0,1, length.out=nx)
		y1 = seq(0,1, length.out=ny)

		k=0
		for (i in 1:nx) {
		for (j in 1:ny) {
		k=k+1
		result2[k,] = c(x1[i], y1[j],sum(result$x==x1[i] & result$y==y1[j])*100/M)
		}
		}

		# plot result
		# true parameters
		parameters = Function.parameter(data1)
		trho00 = parameters$rho00
		trho01 = parameters$rho01
		trho10 = parameters$rho10
		teta = parameters$eta

		#
		data2 = data.frame(x=seq(0, 1, by=0.0001))
		data2$y = Function.mtd_logistic(trho00,trho01,trho10,teta,theta,data2$x)
		data2$y[data2$y<0] = NA
		data2$y[data2$y>1] = NA

        # plot
        if(plot.figure == "Y") {
		par(cex.axis=1.1, cex.lab=1.1, cex.main=1.1, las=1)
		plot(data2$x, data2$y, type="l",xlim=c(0,1),ylim=c(0,1),xlab="Dose (Drug A)",ylab="Dose (Drug B)",xaxt="n",yaxt="n",main="% MTD Recommendation")
		abline(h=seq(0,1,length.out=length(y1)),col="lightgray",lty="dotted")
		abline(v=seq(0,1,length.out=length(x1)),col="lightgray",lty="dotted")
		axis(1, seq(0,1,length.out=length(x1)),1:nx)
		axis(2, seq(0,1,length.out=length(y1)),1:ny)
		text(result2$x, result2$y, format(round(result2$percent,2),digits=2, nsmall=2))
        }
        res = list(result=result, result2=result2)
        invisible(res)
        } # end discrete
    } # end percent
    
  } # End of plot ewoc2simu function

#####################################################################################################################
   # 2.2 summary function for ewoc2simu
       
     print.ewoc2simu <- function(x, ...){
     	
     # input parameter checking
     if(class(x) != "ewoc2simu") stop("wrong class, must be ewoc2simu")
     
     # import the mcmc simulation results
     dlt   <- x$Resp

     # read parameters
     M  = x$parameters$ntrials
     NN = x$parameters$nsamples
     theta = x$parameters$theta
     alpha = x$parameters$alpha
     if(x$type == "discrete") {
     nx = x$parameters$nx
     ny = x$parameters$ny
     data1 = x$parameters$scenario	
     }
     
     # calculate safety characteristics
     propdlt<-numeric()
     for(i in 1:M) propdlt[i] <- sum(dlt[i,]==1)/NN
      
     probdlt2<-sum(propdlt > theta+0.05)/M
     probdlt3<-sum(propdlt > theta+0.1)/M

     # summary
     if(x$type == "continuous") {
     result = data.frame(Value = c(100*sum(propdlt)/M, 100*probdlt2, 100*probdlt3))
     result$Value = format(round(result$Value, 2), digits=2, nsmall=2)
     rownames(result) =c("Average % DLT", "% Trials with DLT rate > theta+0.05", "% Trials with DLT rate > theta+0.1")	
     print("Operating characteristics summarizing trial safety for continuous doses:")
     print(result)
     } else if(x$type == "discrete") {
     
     # define several functions
     pos<-function(a){
	 l<-length(a)
	 y<-numeric()
     for (i in 1:l){
	 if (a[i] >= 0) (y[i] <- a[i])
	 else (y[i]<-0)}
	 y}

     dsquare<-function(a,b){
	 y<-(a-b)^2
	 y}

	 dabs<-function(a,b){
	 y<-abs(a-b)
	 y}

	 dover<-function(a,b,alpha){
     y<-alpha*pos(a-b)+ (1-alpha)*pos(b-a)
	 y}

	 # pt is the true probability of toxicity at a dose combination
	 # pe is estimated probability of selecting dose combination as mtd

	 pt=data1$p
	 pewoc = plot(x, type="percent", plot.figure = "N")$result2$percent / 100
     N<-nx*ny 

     # Accuracy Index
	 Aewocsquare<- 1- (N*(sum(pewoc*dsquare(pt,theta)))/sum(dsquare(pt,theta)))
	 Aewocabs<- 1- (N*(sum(pewoc*dabs(pt,theta)))/sum(dabs(pt,theta)))
	 Aewocover<- 1 - (N*(sum(pewoc*dover(theta,pt,alpha)))/sum(dover(theta,pt,alpha)))
	 
	 # % of selection
	 result = plot(x, type="percent", plot.figure = "N")$result
     data1$id = interaction(round(data1$x,4), round(data1$y,4))
	 result$id= interaction(round(result$x,4),round(result$y,4))
	 bb=merge(result, data1, by="id", all.x=TRUE)
	 bb=bb[order(bb$trial),c(2,3,4,7)]
	 row.names(bb) = NULL
	 names(bb)[1:2] = c("x","y")
	 bb$trial = as.factor(bb$trial)

	 # a simple way
	 k=sum(tapply(bb$p, bb$trial,function(x){all(x %in% data1$p[data1$p >= theta-0.1 & data1$p <= theta+0.1])})) / M *100

     summ = data.frame(Value = c(Aewocsquare, Aewocabs, Aewocover, k, 100*sum(propdlt)/M, 100*probdlt2, 100*probdlt3))
     summ$Value = format(round(summ$Value, 2), digits=2, nsmall=2)
     rownames(summ) =c("Accuracy square discrepancy (sq)", "Accuracy absolute discrepancy (abs)", "Accuracy overdose error (od)", "% Selection", "Average % DLT", "% Trials with DLT rate > theta+0.05", "% Trials with DLT rate > theta+0.1")	
     print("Operating characteristics summarizing trial safety for discrete doses:")
     print(summ)
     }
     
     }  # End of print.ewoc2simu function

##########################################################################################################
# 3. Other functions
#  3.1 Function to plot MTD curve

   mtdcurve = function(rho00,rho01,rho10,eta,theta) {
   	
      # input parameter checking
      if(theta<=0  | theta>=1 )  stop("Error in theta:  theta should be between 0 and 1") 
      if(rho00 <= 0 | rho00 >= 1) stop("rho00 must between 0 and 1")
      if(rho01 <= 0 | rho01 >= 1) stop("rho01 must be 0-1")
      if(rho10 <= 0 | rho10 >= 1) stop("rho10 must be 0-1")
      if(eta <= 0) stop("eta must be positive")

       tmtdx0<-(log(theta/(1-theta)) - log(rho00/(1-rho00))) / (log(rho10/(1-rho10)) - log(rho00/(1-rho00)))
       if (tmtdx0 > 1) {
       x2<-seq(0,1,by=0.01)} else {
	 x2<-seq(0,tmtdx0,by=0.01)}

	 ll<-length(x2)
	 y2<-numeric()
	 for (i in 1:ll){
	 y2[i]<-Function.mtd_logistic(rho00,rho01,rho10,eta,theta,x2[i])}
	 xtemp2<-x2[y2 >= 0 & y2 <= 1]
       y2<-y2[y2 >= 0 & y2 <= 1]
       plot(xtemp2,y2,type="l",xlim=c(0,1),ylim=c(0,1),xlab="Dose (Agent A)", ylab="Dose (Agent B)",main="MTD curve", cex.lab=0.8,cex.axis=0.8)
   
   } # End of function mtdcurve 

##########################################################################################################












