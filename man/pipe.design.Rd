\name{pipe.design}
\alias{pipe.design}

\title{
Dual-Agent Dose Escalation for Phase I Trials using the PIPE Design
}
\description{
Implements the Product of Independent beta Probabilities dose Escalation (PIPE) design for dual-agent Phase I trials as described in Mander A.P. and Sweeting M.J. Statistics in Medicine 2015. doi: 10.1002/sim.6434}
\usage{
pipe.design(N=dim(data)[1]+1, S=1, c, theta, pi=NULL,  prior.med=NULL,	
	prior.ss=NULL, strategy, admis, constraint=NULL, 	epsilon=NULL, 
	mode="sim", data=matrix(nrow=0,ncol=0), a=NULL, b=NULL, 
	alternate=FALSE, 	uppertox.constraint=NULL, stop=NULL)
}

\arguments{
  \item{N}{
	Sample size for the trial. If not specified this is assumed to be one greater than the sample size of the patients recruited thus far (i.e. the number of rows of \code{data} plus one)
}

  \item{S}{
	Number of simulations to conduct. If set to 1 then a single trial is run and invokes different output from \code{\link{print.pipe}} and \code{\link{plot.pipe}} 
}
  \item{c}{
	Cohort size. Must be divisible by \code{N}
}
  \item{theta}{
	The target toxicity probability
}
  \item{pi}{
	A matrix with rows denoting levels of drug A and columns denoting levels of drug B. Each element gives the true probability of the outcome (dose-limiting toxicity) for that dose combination. If omitted then the trial will recommend only for the next dose combination
}
  \item{prior.med}{
	A matrix with rows denoting levels of drug A and columns denoting levels of drug B. Each element gives the prior median probability of the outcome (dose-limiting toxicity) for that dose combination. Should be specified in combination with \code{prior.ss}
}
  \item{prior.ss}{
	A matrix with rows denoting levels of drug A and columns levels of drug B. Each element gives the prior sample size for that dose combination. Should be specified in combination with \code{prior.med}
}
  \item{strategy}{
	A string detailing the dose escalation strategy. Options are "ss" (sample size) and "ss-random" (weighted randomisation by sample size). See details below
}
  \item{admis}{
	A string detailing how to choose admissible doses around the maximum tolerated contour (MTC). Options are "adjacent" and "closest". See details below
}
  \item{constraint}{
	A string detailing any dose-skipping constraint to be applied. Options are "neighbouring" or "no.dose.skip" (see details below). In addition, to prevent diagonal escalation the option "-nodiag" can be appended, i.e. "neighbouring-nodiag" or "no.dose.skip-nodiag". Defaults to NULL
}
  \item{epsilon}{
	Any number between 0 and 1 giving a safety threshold. Doses whose (weighted) posterior probability(dose>MTC)>epsilon will not be experimented on. See details below
}
  \item{mode}{
	A string indicating the running mode of the program. Options are "sim" (the default, used for simulations), "nodlt" (every patient recruited is assumed to be DLT free), or "alldlt" (every patient recruited is assumed to suffer a DLT). 
}
  \item{data}{
	An optional named data frame giving information about dose and toxicity from previously recruited patients. If missing then it is assumed that no data have thus far been collected. Contains the following variables: 
\describe{
     \item{\code{patient}}{ Recruited patient numbers, \code{1,...,n}}
     \item{\code{doseA}}{ Dose levels of recruited patients for drug A}
     \item{\code{doseB}}{ Dose levels of recruited patients for drug B}
     \item{\code{tox}}{ An indicator variable for each patient (1=toxicity, 0=no toxicity)}   
   } 
}
  \item{a}{
	(Optional). A matrix with rows denoting levels of drug A and columns levels of drug B. Each element gives the \code{shape1} parameter (a) from a prior Beta(a,b) distribution. Can be used as an alternative to specifying \code{prior.med} and \code{prior.ss} and must be used in conjunction with \code{b}
}
  \item{b}{
	(Optional). A matrix with rows denoting levels of drug A and columns levels of drug B. Each element gives the \code{shape2} parameter (b) from a prior Beta(a,b) distribution. Can be used as an alternative to specifying \code{prior.med} and \code{prior.ss} and must be used in conjunction with \code{a}
}
  \item{alternate}{
	Logical variable (defaults to FALSE). Should the design always deescalate if above the MTC and escalate if below (subject to admissibility, safety and dose skipping constraints)?
}
\item{uppertox.constraint}{
	Any number between 0 and 1 giving the upper toxicity safety constraint. No dose combination that lies above the most likely contour for this constraint may be dosed. Defaults to NULL
}
\item{stop}{
	(Optional). Value of a stopping threshold. An alternative stopping rule to a fixed sample size specification (through \code{N}) or termination of a trial due to safety constraints (\code{epsilon} or \code{uppertox.constraint}). The trial is stopped if the posterior probability of being > \code{theta} at the lowest dose combination is greater than \code{stop}
}

}
\details{
This function implements a Product of Independent Beta Probabilities (PIPE) Phase I dose-escalation design for a dual-agent trial. An independent Beta prior is placed on the probability of dose limiting toxicity (DLT) for each dose combination. Simple conjugate Bayesian analysis is used to obtain the posterior at each dose combination. This model is considered a `working model` only as it does not account for monotonically increasing toxicity by dose. 

The dose combination assigned to the next cohort of patients is, however, obtained by calculating a maximum tolerated contour (MTC) that bisects the dose space and respects monotonicity. Admissible dose combinations considered for the next cohort are those that are either `adjacent' or `closest' to the MTC (see Mander A.P. and Sweeting M.J., Statistics in Medicine 2015) and can be specified using the \code{admis} argument. If there is more than one dose combination then either the one with the smallest sample size (least experimented on) is chosen, by specifying \code{strategy="ss"}, or admissible doses are chosen at random, weighted by the inverse of the sample size (\code{strategy="ss-random"}).

Admissible doses can be constrained to be within a set of neighbouring doses to the current dose combination (\code{constraint="neighbouring"}) or within the set of doses consisting of up to one dose level above any previously experimented drug A and drug B levels (\code{constraint="no.dose.skip"}). To additionally prevent diagonal escalation (e.g. where both drugs are escalated together), use either (\code{constraint="neighbouring-nodiag"}) or (\code{constraint="no.dose.skip-nodiag"}).

In addition, a safety constraint can be added, to avoid potential overdosing. Dose combinations whose expected probability of being above the MTC is greater than a threshold \code{epsilon} are excluded from the admissible doses for the next cohort. The expected probability is calculated by averaging over all possible MTCs, weighted by their posterior probabilities. An alternative safety constraint can be used by specifying \code{uppertox.constraint=c}, where c is the upper toxicity level from which a contour is derived (e.g. c=0.6 calculates the most likely contour for a 60% chance of DLT). No doses above this contour will be used. 

The programe will either run as a simulation (the default), or by using a deterministic set of outcomes (all DLTs / no DLTs) specified using \code{mode="alldlt"} and \code{mode="nodlt"}, respectively. 

If the user specifies a data-frame using the \code{data} argument, then every simulated trial will begin with the specified doses and DLT outcomes. This option is useful for assess operating characteristics mid-way through the trial, or just to obtain the dose for the next cohort.

A trial of size \code{N} is generally run, except if no doses are deemed safe according to one of the two safety constraints described above, at which point the trial will stop. Users may additionally specify a stopping rule using \code{stop=s}, where the trial is stopped if the posterior probability that the lowest dose combination is above the target toxicity probability is greater than s. Note, this stopping rule relies only on the prior and data collected for the first dose combination and does not borrow strength from other doses. Hence it should be used with caution.

A final option is to force deescalation if above the MTC and escalation if below (subject to meeting other specified constraints) by specifying \code{alternate=TRUE}. This constraint ensures the design follows the properties of `coherence` (Cheung YK. Biometrika 2005). 
}
\value{
\code{pipe.design} returns an object of class "pipe" or "pipe.sim"; the latter occuring when more than one simulation has been conducted.
The function \code{\link{print}} (i.e. \code{\link{print.pipe}} or \code{\link{print.pipe.sim}}) can be used to obtain summary information about the design used, the data observed, the maximum tolerated contour estimate, and the next recommended dose level. 

An object of class "pipe" (when N=1) or "pipe.sim" (when N>1) is a list with the following components:
\describe{
	\item{\code{r.sim}}{A list of matrices showing the number of DLTs for each dose combination in each simulation}
	\item{\code{n.sim}}{A list of matrices showing the number of patients dosed at each combination in each simulation}
	\item{\code{rec.i.sim}}{An S x (N/c) matrix with elements representing the level of drug A for each cohort}
	\item{\code{rec.j.sim}}{An S x (N/c) matrix with elements representing the level of drug B for each cohort}
	\item{\code{exp}}{A matrix showing the experimentation proportions for each dose combination averaged across simulations}
	\item{\code{rec}}{A matrix showing the proportion of times each dose combination is recommended, where the denominator is the total number of recommended phase II doses over the simulations}
	\item{\code{dlts}}{A vector giving the average number of DLTs occuring in each simulated trial}
	\item{\code{mat.list}}{A list of binary matrices showing the most likely MTC before each cohort is recruited. Output only if \code{S=1}}
	\item{\code{uppermat.list}}{A list of binary matrices showing the most likely contour for \code{uppertox.constraint} before each cohort is recruited. Output only if \code{S=1} and \code{uppertox.constraint} is used}
	\item{\code{uppermat2.list}}{A list of binary matrices showing the doses that fail the constraint imposed by \code{epsilon}. Output only if \code{S=1} and \code{epsilon} is used}
	\item{\code{r.list}}{A list of matrices showing the cumulative number of DLTs for each dose combination before each cohort is recruited. Output only if \code{S=1}}
	\item{\code{n.list}}{A list of matrices showing the cumulative number of patients at each dose combination before each cohort is recruited. Output only if \code{S=1}}
	\item{\code{n.rpII}}{A vector listing the number of doses recommended for Phase II from each simulated trial}
	\item{\code{no.not.treated}}{The total number of patients not treated (in the trial) accumulated across the simulated trials due to a trial stopping early}
	\item{\code{pi}}{The true probability of DLT for each dose combination, as specified in the input argument}
	\item{\code{theta}}{The target toxicity probability, as specified in the input argument}
	\item{\code{rpII.list}}{A list giving the recommended Phase II dose combinations for each simulation}
}
}
\references{
 Mander A.P., Sweeting M.J. A product of independent beta probabilities dose escalation design for dual-agent phase I trials. \emph{Statistics in Medicine} (2015) 34(8): 1261--1276.
 
Cheung Y.K. Coherence principles in dose-finding studies. \emph{Biometrika} (2005) 92(4): 863--873.

}

\note{
	The use of the \code{uppertox.constraint} safety constraint is not documented in Mander and Sweeting, Statistics in Medicine 2015, and some preliminary investigations have shown that this constraint can cause rigidity in the escalation process, and therefore should be used with caution.
}

\author{
Michael Sweeting \email{mjs212@medschl.cam.ac.uk} (University of Cambridge, UK)
}

\seealso{
\code{\link{print.pipe}}, \code{\link{print.pipe.sim}}, \code{\link{plot.pipe}}, \code{\link{plot.pipe.sim}}
}
\examples{

## Reproducing Figure 5 from Mander and Sweeting, Statistics in Medicine 2015.
## The true probability of DLT for a 6x6 grid of dose combinations 
## (Scenario 3 from Mander and Sweeting, Statistics in Medicine 2015)
pi <-matrix(c(0.02,0.10,0.20,0.30,0.35,0.45,0.06,0.14,0.24,0.34,0.39,0.49,0.12,0.20,
0.30,0.40,0.45,0.55,0.17,0.25,0.35,0.45,0.50,0.60,0.22,0.30,0.40,0.50,0.60,0.70,0.30,
0.38,0.48,0.58,0.68,0.78),nrow=6,ncol=6)
## Using a weak prior with prior sample size 1/36 for each dose combination and prior
## median of Scenario 1
prior.med<-matrix(c(0.02,0.03,0.06,0.10,0.18,0.23,0.03,0.05,0.09,0.13,0.21,0.27,0.06,
0.09,0.14,0.18,0.26,0.30,0.11,0.14,0.18,0.23,0.30,0.36,0.18,0.21,0.26,0.30,0.39,0.44,
0.23,0.27,0.3,0.36,0.44,0.49),nrow=6,ncol=6)
## Using a neighbouring escalation constraint
## Selecting the closest admissible doses
## Using a safety constraint with epsilon = 0.8
## Running a single trial with sample size of 40, cohort size 2 
## and target toxicity of 0.3
set.seed(700)
example<-pipe.design(N=40,S=1,c=2,theta=0.3,pi=pi,prior.med=prior.med,
prior.ss=matrix(1/36,ncol=6,nrow=6),strategy="ss",constraint="neighbouring",
epsilon=0.8,admis="closest",alternate=FALSE)
example
plot(example)

\dontrun{
## Conducting a simulation study (1000 trials) to investigate 
## operating characteristics for Scenario 3 
## (as reported in Table III in Mander and Sweeting, Statistics in Medicine 2015)
set.seed(262)
scen3<-pipe.design(N=40,S=1000,c=2,theta=0.3,pi=pi,prior.med=prior.med,
prior.ss=matrix(1/36,ncol=6,nrow=6),strategy="ss",constraint="neighbouring",
epsilon=0.8,admis="closest",alternate=FALSE)
scen3
}

}


