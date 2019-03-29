\name{plot.ewoc2simu}
\alias{plot.ewoc2simu}
\title{EWOC for 2 drugs combination trial design characteristics}
\description{Function to plot the trial design characteristics from EWOC 2 drugs combination simulation results}
\usage{
\method{plot}{ewoc2simu}(x, type = "MTD", conf.reg=0.9, plot.figure="Y",...)
}
\arguments{
  \item{x}{an object of class "ewoc2simu", usually a result of a call to ewoc2simu}
  \item{type}{a character indicating the type of plots a user requests, could be "MTD", "bias", or "percent". For discrete simulations, "bias" is not available}
  \item{conf.reg}{confidence level that controls the region of the doses from the last trial in the MTD plot}
  \item{plot.figure}{a character indicating whether user wants the plot, 'Y' would be yes, otherwise would be no. It's mainly for internal uses}
  \item{...}{arguments passed to or from methods}
}

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
