\name{pdlt}
\alias{pdlt}
\title{Generating probability of DLT based on the EWOC2 model}
\description{Generating probabiity of DLT based on the EWOC 2 drugs combination model}
\usage{
pdlt(rho00, rho01, rho10, eta, theta, x, y)
}
\arguments{
  \item{rho00}{a numeric value indicating the true value of the parameter rho00, the probability of DLT when the levels of drugs A and B are both 0}
  \item{rho01}{a numeric value indicating the true value of the parameter rho01, the probability of DLT when the levels of drugs A and B are 0 and 1, respectively}
  \item{rho10}{a numeric value indicating the true value of the parameter rho10, the probability of DLT when the levels of drugs A and B are 1 and 0, respectively}
  \item{eta}{a numeric value indicating the true value of the eta, the interaction parameter}
  \item{theta}{a numerical value defining the proportion of expected patients to experience a medically unacceptable, dose-limiting toxicity (DLT) if administered the MTD}
  \item{x}{a numeric value of dose level for drug A}
  \item{y}{a numeric value of dose level for drug B}
}
\value{a numeric value indicating the probability of DLT with doses from input based on the logistic model}

\examples{
pdlt(rho00=0.01, rho01=0.2, rho10=0.9, eta=20, theta=0.2, x=0.2, y=0.3)
}
