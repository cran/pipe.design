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
\method{print}{pipe.sim}(x, pi = x$pi, cut.points = c(0, 15, 25, 35, 45.1, 100)/100,
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
Cutpoints of toxicity for which the operating characteristics are to be categorised. 
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
	\item{Experimentation percentages}{Percentage of patients recruited to each true region of toxicity (as specified by \code{cutpoints}), across the simulated trials}
	\item{Recommendation percentages}{Percentage of trials that recommend doses within each true region of toxicity (as specified by \code{cutpoint})}
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