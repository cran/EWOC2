\name{ewoc2simu}
\alias{ewoc2simu}
\alias{ewoc2simu.default}
\title{Generic EWOC2 simulation}
\description{Generic function for simulating EWOC trials for 2 drugs combination}
\usage{
ewoc2simu(ntrials, nsamples, type, trho00, trho01, trho10, teta, nx, ny, tp, 
Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, alpha, theta, vai, a01, 
b01, a10, b10, a00, b00, a, b, delta1x, delta1y, burn, mm, delta1, seed)

\method{ewoc2simu}{default}(ntrials, nsamples, type, trho00, trho01, trho10, teta, nx, ny, tp, 
Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, alpha, theta, vai, a01,
b01, a10, b10, a00, b00, a, b, delta1x, delta1y, burn=4000, mm=2000, delta1=0.05, seed)
}
\arguments{
  \item{ntrials}{a number indicating the number of trials to be simulated}
  \item{nsamples}{a number indicating the number of patients enrolled for each clinical trial}
  \item{type}{a character indicating the type of design, could be 'continous' or 'discrete' or their initials}
  \item{trho00}{a numeric value indicating the true value of the parameter rho00, the probability of DLT when the levels of drugs A and B are both 0}
  \item{trho01}{a numeric value indicating the true value of the parameter rho01, the probability of DLT when the levels of drugs A and B are 0 and 1, respectively}
  \item{trho10}{a numeric value indicating the true value of the parameter rho10, the probability of DLT when the levels of drugs A and B are 1 and 0, respectively}
  \item{teta}{a numeric value indicating the true value of the eta, the interaction parameter}
  \item{nx}{a numeric value indicating the number of dose levels for drug A. It's only necessary if type = 'discrete'}
  \item{ny}{a numeric value indicating the number of dose levels for drug B. It's only necessary if type = 'discrete'}
  \item{tp}{a numerical vector indicating the true probabilities of DLT at each dose combinations, the order is by Drug B first, only necessary if type = 'discrete'}
\item{Min.Dose.A}{a numeric value defining the lower bound of the support of the MTD for drug A}
  \item{Max.Dose.A}{a numeric value defining the upper bound of the support of the MTD for drug A}
  \item{Min.Dose.B}{a numeric value defining the lower bound of the support of the MTD for drug B}
  \item{Max.Dose.B}{a numeric value defining the upper bound of the support of the MTD for drug B}
  \item{alpha}{a numerical value defining the probability that dose selected by EWOC is higher than the MTD.}
  \item{theta}{a numeric value defining the proportion of expectd patients to experience a medically unacceptable, dose-limiting toxicity (DLT) if administered the MTD.}
  \item{vai}{a numeric value indicating variable alpha increment for each new cohort}
  \item{a01}{a numeric value for beta prior distribution associated with parameter rho01}
  \item{b01}{a numeric value for beta prior distribution associated with parameter rho01}
  \item{a10}{a numeric value for beta prior distribution associated with parameter rho10}
  \item{b10}{a numeric value for beta prior distribution associated with parameter rho10}
  \item{a00}{a numeric value for beta prior distribution associated with parameter rho00}
  \item{b00}{a numeric value for beta prior distribution associated with parameter rho00}
  \item{a}{a numeric value for gamma prior distribution associated with parameter eta}
  \item{b}{a numeric value for gamma prior distribution associated with parameter eta}
  \item{delta1x}{Maximum dose escalation at each step for drug A, the default is 0.2*(Max.Dose.A-Min.Dose.A if not assigned)}
  \item{delta1y}{Maximum dose escalation at each step for drug B, the default is 0.2*(Max.Dose.B-Min.Dose.B if not assigned)}
  \item{burn}{Number of iterations for adaption, see n.adapt in jags.model for detail}
  \item{mm}{Number of iterations to monitor, see n.iter in code.samples for detail}
  \item{delta1}{Threshold for toxicity}
  \item{seed}{a numeric value used in random number generation }
}
\value{
\item{type}{same as input parameter type}
\item{parameters}{list of input parameters}
\item{priors}{list of prior parameters}
\item{Dose.A}{a matrix ntrials x nsamples containing the doses of drug A assigned for each patient in a trial and each trial in the simulation}
\item{Dose.B}{a matrix ntrials x nsamples containing the doses of drug B assigned for each patient in a trial and each trial in the simulation}
\item{Resp}{a matrix ntrials x nsamples containing ones and zeros indicating the occurance of DLT (1) and the absence of DLT (0) for each patient in the trial and each trial in the simulation}
\item{rho00}{a numeric vector ntrials x 1 containing the estimated rho00 parameter for each trial in the simulation}
\item{rho01}{a numeric vector ntrials x 1 containing the estimated rho01 parameter for each trial in the simulation}
\item{rho10}{a numeric vector ntrials x 1 containing the estimated rho10 parameter for each trial in the simulation}
\item{eta}{a numeric vector ntrials x 1 containing the estimated eta parameter for each trial in the simulation}
\item{postlow}{a matrix ntrials x nsamples/2 containing posterior probability of DLT at lower doses (both 0 for durg A and B) at each step in a trial and each trial in the simulation}
\item{postdlts}{a matrix (nx x ny x ntrials) x 4 containing posterior probability of DLT at each dose combination sets in each trial in the simulation. This is used to test whether or not a discrete set of MTDs was selected from a continous MTD curve is kept or dropped. It's avaiable only when type = 'discrete'}
}
\references{Tighiouart M, Li Q and Rogatko A. A Bayesian adaptive design for estimating the maximuym tolerated dose curve using drug combinations in cancer phase I clinical trials. Statistics in Medicine. 2017, 36: 280-290.}
\examples{
\dontshow{
# continous 
test1 = ewoc2simu(ntrials=1, nsamples=4, type="c", trho00=0.01,trho01=0.2, trho10=0.9,teta=20, 
Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, a01=1,b01=1,
a10=1,b10=1, a00=1,b00=1,a=0.8,b=0.0384)

print(test1)
plot(test1, type="MTD")
plot(test1, type="bias")
plot(test1, type="percent")

# discrete
tp = c(0.03,0.05,0.08,0.05,0.08,0.13,0.08,0.13,0.2,0.13,0.2,0.29,0.2,0.29,0.4,0.29,0.4,0.53)
test2 = ewoc2simu(ntrials=1, nsamples=4, type="d", nx=6, ny=3, tp=tp, 
Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, 
a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)

print(test2)
plot(test2, type="MTD")
plot(test2, type="percent")
}
\donttest{
# continous 
test1 = ewoc2simu(ntrials=10, nsamples=40, type="c", trho00=0.01,trho01=0.2, trho10=0.9,teta=20, 
Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, a01=1,b01=1,
a10=1,b10=1, a00=1,b00=1,a=0.8,b=0.0384)

print(test1)
plot(test1, type="MTD")
plot(test1, type="bias")
plot(test1, type="percent")

# discrete
tp = c(0.03,0.05,0.08,0.05,0.08,0.13,0.08,0.13,0.2,0.13,0.2,0.29,0.2,0.29,0.4,0.29,0.4,0.53)
test2 = ewoc2simu(ntrials=10, nsamples=40, type="d", nx=6, ny=3, tp=tp, 
Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1, alpha=0.25, theta=0.20, 
a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)

print(test2)
plot(test2, type="MTD")
plot(test2, type="percent")
}
}