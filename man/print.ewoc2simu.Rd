\name{print.ewoc2simu}
\alias{print.ewoc2simu}
\title{Summarizing EWOC2 simulation results}
\description{Summarizing EWOC2 simulation results}
\usage{
\method{print}{ewoc2simu}(x, ...)
}
\arguments{
  \item{x}{an object of class "ewoc2simu", usually, a result of a call to ewoc2simu}
  \item{...}{arguments passed to or from methods}
}
\value{a data.frame of 7 x 1 with row represent Accuracy square discrepancy (sq), Accuracy absolute discrepancy (abs), Accuracy overdose (od), percent Selection, Average percent DLT, percent Trials with DLT rate > theta+0.05, percent Trials with LDT rate > theta+0.1}
\references{Tighiouart M, Li Q and Rogatko A. A Bayesian adaptive design for estimating the maximuym tolerated dose curve using drug combinations in cancer phase I clinical trials. Statistics in Medicine. 2017, 36: 280-290.}
\examples{
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
