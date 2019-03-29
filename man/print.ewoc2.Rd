\name{print.ewoc2}
\alias{print.ewoc2}
\title{Summarizing EWOC2 next doses results}
\description{Summarizing EWOC2 next doses result}
\usage{
\method{print}{ewoc2}(x, ...)
}
\arguments{
  \item{x}{an object of class "ewoc2", usually, a result of a call to ewoc2}
  \item{...}{arguments passed to or from methods}
}
\value{a data.frame  of 2 x 4 with columns for cohort, patients, recommended dose of drug A and recommended dose of drug B for next cohort or 2 patients}
\references{Tighiouart M, Li Q and Rogatko A. A Bayesian adaptive design for estimating the maximuym tolerated dose curve using drug combinations in cancer phase I clinical trials. Statistics in Medicine. 2017, 36: 280-290.}
\examples{
test = ewoc2(dose.a=c(0,0),dose.b=c(0,0),resp=c(0,0),theta=0.33,alpha=0.25, Min.Dose.A=0, 
Max.Dose.A=1, Min.Dose.B=0, Max.Dose.B=1,a01=1,b01=1,a10=1,b10=1,a00=1,b00=1,a=0.8,b=0.0384)
print(test)
}
