\name{balanced}
\alias{balanced}

\title{
Compute the M-quantile extension of the mixed linear model for BALANCED DATA ONLY
}
\description{
Compute the M-quantile extension of the mixed linear model for BALANCED DATA ONLY
}
\usage{
balanced(model, k, q, init = list(d = 1, v = 1), store = FALSE, 
         acc = 1e-05, niter = 30, nonzero = TRUE)
}

\arguments{
  \item{model}{
a \code{saemodel} object
}
  \item{k}{
scalar; robustness-tuning constant of the asymmetric Huber psi-function
}
  \item{q}{
scalar; quantile \code{(0,1)}; 
}
  \item{init}{
initial values
}
  \item{store}{
logical; if \code{store=TRUE}, the function returns the iteration-specific estimates
}
  \item{acc}{
scalar; accuarcy measure used in the termination rule of the algorithm
}
  \item{niter}{
scalar; number of iterations
}
  \item{nonzero}{
logical; if \code{nonzero=FALSE}, the method also uses negative variance estimates (sic!)
}
}
\details{
This function computes M-quantile-type estimates of the fixed effects and the random-effect variances of the balanced-data mixed linear model. 

The model must be defined by means of the \code{saemodel} function in the \pkg{rsae} package.
}
\value{
\item{beta}{vector fixed-effects estimates}
\item{s_e}{model-error std. dev. (not variance)}
\item{s_a}{area-level-random-effect std. dev. (not variance)}
\item{niter}{number of iterations used to compute the estimates}

}
\references{
Tzavidis et al. (2011)
}
\author{
Tobias Schoch
}

