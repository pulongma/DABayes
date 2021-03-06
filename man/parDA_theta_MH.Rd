% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parDA_theta_MH.R
\name{parDA_theta_MH}
\alias{parDA_theta_MH}
\title{Update parameters from the model in the MCMC algorithm}
\usage{
parDA_theta_MH(outW, theta, prior_theta, ADtheta, loglik_old, chisq_old, N, Lj,
  tol = 1e-10)
}
\arguments{
\item{outW}{a 7 by 1 list containing precomputed quantities associated with W}

\item{theta}{a list of 3 elements contained in parameters in MCMC algorithm}

\item{prior_theta}{quantities (parameters) in prior dist.}

\item{ADtheta}{a 6 by 1 list containing quantities to adjust mean and 
covariance in proposal dist}

\item{loglik_old}{a real number representing likelihood for given parameters}

\item{chisq_old}{a real number representing the quadratic form}

\item{N}{number of temperature ensemble members}

\item{Lj}{an m by 1 vector containing the number of runs for each 
forcing scenario}

\item{tol}{a very small value to avoid numerical instability}
}
\value{
a list of updated parameters, loglikelihood, chisq statistics (quadratic form),
and prior
}
\description{
This function updates parameters for the Bayesian statistical 
model. gamma is updated with a discrete uniform prior; beta is updated with
noninformative prior 1 and random-walk proposal; log of sigma (variance in W)
is updated with normal prior and random-walk proposal; lambda is updated 
with normal prior and random-walk proposal. This function is used in 
parallel computing of the MCMC algorithm.
}
\author{
Pulong Ma <mpulong@gmail.com>
}
\seealso{
DA_theta_MH
}
\keyword{models}

