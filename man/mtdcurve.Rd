\name{mtdcurve}
\alias{mtdcurve}
\title{Generating MTD curve based on logistic model for two drugs}
\description{Generating MTD curve based on logistic model for two drugs}
\usage{
mtdcurve(rho00, rho01, rho10, eta, theta)
}
\arguments{
  \item{rho00}{a numeric value indicating the true value of the parameter rho00, the probability of DLT when the levels of drugs A and B are both 0}
  \item{rho01}{a numeric value indicating the true value of the parameter rho01, the probability of DLT when the levels of drugs A and B are 0 and 1, respectively}
  \item{rho10}{a numeric value indicating the true value of the parameter rho10, the probability of DLT when the levels of drugs A and B are 1 and 0, respectively}
  \item{eta}{a numeric value indicating the true value of the eta, the interaction parameter}
  \item{theta}{a numerical value defining the proportion of expected patients to experience a medically unacceptable, dose-limiting toxicity (DLT) if administered the MTD}
}
\value{none}

\examples{
mtdcurve(rho00=0.01, rho01=0.2, rho10=0.9, eta=20, theta=0.2)
}
