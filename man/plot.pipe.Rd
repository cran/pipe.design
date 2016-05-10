\name{plot.pipe}
\alias{plot.pipe}
\alias{plot.pipe.sim}
\title{
Plot escalation or operating characteristics for trials conducted using the Product of Independent Beta Probabilities Escalation design 
}
\description{
Plot method for a trial or series of trials conducted using a \code{\link{pipe.design}} model. 
}
\usage{
\method{plot}{pipe}(x, type = "b", pi = x$pi, theta = x$theta, epsilon.line = TRUE, 
    uppertox.constraint.line = FALSE, add.empirical.data = FALSE, ...)
\method{plot}{pipe.sim}(x, pi = x$pi, theta = x$theta, plot = "both", ...) 
}
\arguments{
  \item{x}{
An object of class "pipe" or "pipe.sim" as returned by \code{\link{pipe.design}}
}
 \item{type}{
 	The type of information to be plotted for a single PIPE trial. Options are "r" (number of outcomes only), "n" (number of patients only), or "b" (both number of outcomes and patients). The default is to plot both the number of outcomes and patients for each dose combination, by using symbols and colour filling, respectively.
 }
 \item{pi}{
A matrix with rows denoting levels of drug A and columns denoting levels of drug B. Each element gives the true probability of the outcome (dose-limiting toxicity) for that dose combination. If omitted then the true probabilities of the outcome will be taken from that used when creating the \code{x} object.
}
\item{theta}{
	The target toxicity probability. If omitted then \code{theta} is taken from that used when creating the \code{x} object.
}
\item{epsilon.line}{
	Should the toxicity safety contour specified by \code{epsilon} in \code{\link{pipe.design}} be shown? Defaults to TRUE.
}
\item{uppertox.constraint.line}{
	Should the toxicity safety contour specified by \code{uppertox.constraint} in \code{\link{pipe.design}} be shown? Defaults to FALSE.
}
\item{add.empirical.data}{
	Should an extra plot be produced with the empirical posterior median toxicity probabilities (from the independent beta distributions), and 95% credible intervals? Defaults to FALSE.
}
 \item{plot}{
	The type of operating characteristics to plot for a series of simulated PIPE trials. Options are "exp" (experimentation percentages), "rec" (recommendation percentages), or "both" (both experimentation and recommendation percentages side-by-side). The default is "both".
}
 \item{...}{Further arguments passed to or from other methods}
}
\details{
If a single trial is conducted, then the \code{\link{plot}} function currently produces a series of matrix plots by cohort showing the cumulative outcome and sample size information for each dose combination. Also shown as a solid black line is the most likely MTC before outcomes from each cohort are observed and, if specified, one or more toxicity safety contours shown in red. The final cohort shows the most likely MTC at the end of the trial together with the true MTC (green dashed line) and the recommended Phase II doses (blue squares), which are always doses below the most likely MTC that have been experimented on (see Mander and Sweeting, Statistics in Medicine 2015 for more details).

If multiple PIPE trials have been simulated then experimentation and recommendation percentages are shown as a heat plot over the dose space together with the true MTC (green dashed line).

}


\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk}
}
\references{
 Mander A.P., Sweeting M.J. A product of independent beta probabilities dose escalation design for dual-agent phase I trials. \emph{Statistics in Medicine} (2015) 34(8): 1261--1276.
 }
\seealso{
\code{\link{pipe.design}}, \code{\link{print.pipe}}
}
