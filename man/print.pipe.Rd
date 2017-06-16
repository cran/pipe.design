\name{print.pipe}
\alias{print.pipe}
\alias{print.pipe.sim}
\title{
Print information regarding a trial conducted using the Product of Independent Beta Probabilities Escalation design
}
\description{
Print method for a trial or series of trials conducted using a \code{\link{pipe.design}} model. 
}
\usage{
\method{print}{pipe}(x, \dots)
\method{print}{pipe.sim}(x, pi = x$pi, 
  cut.points = unique(c(0,pmin(1,pmax(0,seq(x$theta-0.15,x$theta+0.15,by=0.1))),1)), 
  digits = 1, print = TRUE, \dots)
}
\arguments{
  \item{x}{
An object of class "pipe" or "pipe.sim" as returned by \code{\link{pipe.design}}
}

 \item{pi}{
A matrix with rows denoting levels of drug A and columns denoting levels of drug B. Each element gives the true probability of the outcome (dose-limiting toxicity) for that dose combination. If omitted then the true probabilities of the outcome will be taken from that used when creating the \code{x} object.
}
 \item{cut.points}{
Cutpoints of toxicity for which the operating characteristics are to be categorised. Default is [0,x$theta-0.1) [x$theta-0.1,x$theta+0.1) [x$theta+0.1,1]
}
 \item{digits}{
	The number of decimal places to print the operating characteristics 
}
 \item{print}{If TRUE then the experimentation and recommendation percentages are printed to the output}
\item{...}{Further arguments passed to or from other methods}
 
}
\details{
If a single trial is conducted, then the \code{\link{print}} function currently produces summary information about the data observed, current best estimate of the MTC and the upper toxicity constraint contour. 
If a simulation study is conducted, then the following operating characteristics are printed:
\describe{
	\item{Experimentation percentages}{Percentage of patients recruited to each true region of toxicity (as specified by \code{cut.points}), across the simulated trials}
	\item{Recommendation percentages}{Percentage of times that recommend doses lie within each true region of toxicity (as specified by \code{cut.points}). As more than one dose combination can be recommended in any trial, the denominator is the total number of recommended phase II doses over the simulations}
	\item{Percentage of times a trial recommends k doses for Phase II}{Provides the percentage of time that each trial recommends 0, 1, ..., k doses for Phase II experimentation. 0 doses indicates that the trial stopped early and recommended no doses.}
}

}


\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk}
}
\references{
 Mander A.P., Sweeting M.J. A product of independent beta probabilities dose escalation design for dual-agent phase I trials. \emph{Statistics in Medicine} (2015) 34(8): 1261--1276.
 }
\seealso{
\code{\link{pipe.design}}, \code{\link{plot.pipe}}
}
