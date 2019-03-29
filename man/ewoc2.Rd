\name{ewoc2}
\alias{ewoc2}
\alias{ewoc2.default}
\title{Escalation With Overdose Control for two drugs combination}
\description{Finding the doses of next cohort for a phase I clinical trial based on Escalation with Overdose Control (EWOC) design considering the classic parametrization for binary response and two agents.}
\usage{
ewoc2(dose.a, dose.b, resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, 
a01, b01, a10, b10, a00, b00, a, b, delta1x, delta1y, burn, mm, delta1)

\method{ewoc2}{default}(dose.a, dose.b, resp, theta, alpha, Min.Dose.A, Max.Dose.A, Min.Dose.B, Max.Dose.B, 
a01, b01, a10, b10, a00, b00, a, b, delta1x, delta1y, burn=4000, mm=2000, delta1=0.05)
}
\arguments{
  \item{dose.a}{a numeric vector of allowable doses for drug A}
  \item{dose.b}{a numeric vector of allowable doses for drug B}
  \item{resp}{a numeric vector of allowable responses, 0 or 1}
  \item{theta}{a numeric value defining the proportion of expectd patients to experience a medically unacceptable, dose-limiting toxicity (DLT) if administered the MTD.}
  \item{alpha}{a numerical value defining the probability that dose selected by EWOC is higher than the MTD.}
  \item{Min.Dose.A}{a numeric value defining the lower bound of the support of the MTD for drug A}
  \item{Max.Dose.A}{a numeric value defining the upper bound of the support of the MTD for drug A}
  \item{Min.Dose.B}{a numeric value defining the lower bound of the support of the MTD for drug B}
  \item{Max.Dose.B}{a numeric value defining the upper bound of the support of the MTD for drug B}
  \item{a01}{a numeric value for beta prior distribution associated with parameter rho01}
  \item{b01}{a numeric value for beta prior distribution associated with parameter rho01}
  \item{a10}{a numeric value for beta prior distribution associated with parameter rho10}
  \item{b10}{a numeric value for beta prior distribution associated with parameter rho10}
  \item{a00}{a numeric value for beta prior distribution associated with parameter rho00}
  \item{b00}{a numeric value for beta prior distribution associated with parameter rho00}
  \item{a}{a numeric value for gamma prior distribution associated with parameter eta}
  \item{b}{a numeric value for the gamma prior distribution associated with the parameter eta}
  \item{delta1x}{Maximum dose escalation at each step for drug A, the default is 0.2*(Max.Dose.A-Min.Dose.A if not assigned)}
  \item{delta1y}{Maximum dose escalation at each step for drug B, the default is 0.2*(Max.Dose.B-Min.Dose.B if not assigned)}
  \item{burn}{Number of iterations for adaption, see n.adapt in jags.model for detail}
  \item{mm}{Number of iterations to monitor, see n.iter in code.samples for detail}
  \item{delta1}{Threshold for toxicity}
}
\value{
\item{data}{a data frame containing the current doses and responses set}
\item{parameters}{list of input parameters}
\item{priors}{list of prior parameters}
\item{nextdose.x}{the next recommended doses for drug A}
\item{nextdose.y}{the next recommended doses for drug B}
}
\references{Tighiouart M, Li Q and Rogatko A. A Bayesian adaptive design for estimating the maximuym tolerated dose curve using drug combinations in cancer phase I clinical trials. Statistics in Medicine. 2017, 36: 280-290.}
\examples{
test = ewoc2(dose.a=c(0,0),dose.b=c(0,0),resp=c(0,0),theta=0.33,alpha=0.25, 
Min.Dose.A=0, Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1,a01=1,b01=1,a10=1,b10=1,
a00=1,b00=1,a=0.8,b=0.0384)
print(test)
}
