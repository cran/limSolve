\name{xsample}
\alias{xsample}
\title{Randomly samples an under- or overdetermined problem with linear equality and inequality conditions}
\description{Bayesian samping of linear problems with linear inequality conditions. 
 Uses either the "hit and run" (or directions algorithm) or a mirroring technique for sampling 
 
 The monte carlo method produces a sample solution for 
 \deqn{Ex=f}
 \deqn{Ax\sim=B}
 \deqn{Gx>=h}
 where Ex=F have to be met exactly, and Ax$\sim$=B approximately.}
\usage{xsample(A=NULL, B=NULL, E=NULL, F=NULL, 
G=NULL, H=NULL, sdB=1, iter=3000, type = "mirror", 
jmp = 0.1, tol=sqrt(.Machine$double.eps), 
x0=NULL, fulloutput = FALSE)}
\arguments{
  \item{A }{numeric matrix containing the coefficients of the (approximate) equality constraints, \eqn{Ax~=B}}
  \item{B }{numeric vector containing the right-hand side of the (approximate) equality constraints}
  \item{E }{numeric matrix containing the coefficients of the (exact) equality constraints, \eqn{Ex=F}}
  \item{F }{numeric vector containing the right-hand side of the (exact) equality constraints}
  \item{G }{numeric matrix containing the coefficients of the inequality constraints, \eqn{Gx>=H}}
  \item{H }{numeric vector containing the right-hand side of the inequality constraints}
  \item{sdB }{vector with standard deviation on B, used for weighing; default = equal weights}
  \item{iter }{integer determining the number of iterations}
  \item{type }{type of algorithm: one of: "mirror", (mirroring algorithm), "rda" (random directions algorithm) or "cda" (coordinates directions algorithm)}
  \item{jmp }{jump length of the transformed variables q: \eqn{x=p+Zq} (only if \code{type}=="mirror")}
  \item{tol }{tolerance for equality and inequality constraints; numbers whose absolute value is smaller than \code{tol} are set to zero}
  \item{x0 }{initial (particular) solution}
  \item{fulloutput }{if TRUE, also outputs the transformed variables q}
}
\value{
  a list containing:
  \item{X }{matrix whose rows contain the sampled values of x}
  \item{acceptedratio }{ratio of acceptance (i.e. the ratio of the accepted runs / total iterations)}
  \item{p }{only outputted if fulloutput is TRUE: probability vector for all samples (e.g. one value for each row of \code{X})}
  \item{Q }{only outputted if fulloutput is TRUE: the transformed samples Q}
}
\author{Karel Van den Meersche<k.vdMeersche@nioo.knaw.nl>
  \cr Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Sample the underdetermined Mink diet problem
E <- rbind(Minkdiet$Prey,rep(1,7))
F <- c(Minkdiet$Mink,1)

pairs(xsample(E=E,F=F,G=diag(7),H=rep(0,7),iter=1000)$X,
     main="Minkdiet 1000 solutions")

# sample the overdetermined Chemtax problem
Nx     <-nrow(Chemtax$Ratio)
  
# equations that have to be met exactly Ex=f: 
# sum of all fraction must be equal to 1.
E <- rep(1,Nx)
F <- 1
  
# inequalities, Gx>=h:
# all fractions must be positive numbers
G <- diag(nrow=Nx)
H <- rep(0,Nx)
  
# equations that must be reproduced as close as possible, Ax ~ b
# the field data; the input ratio matrix and field data are rescaled
A     <- t(Chemtax$Ratio/rowSums(Chemtax$Ratio))
B     <- Chemtax$Field/sum(Chemtax$Field)
  
# Sample
pairs(xsample(A=A,B=B,E=E,F=F,G=G,H=H,iter=1000)$X,
     main="Chemtax 1000 solutions")}
\details{The algorithm proceeds in two steps.
	\item (1) the equality constraints \eqn{Ex=F} are eliminated, and the system \eqn{Ex=f}, \eqn{Gx>=h} is rewritten as \eqn{G(p+Zq)>= h}, 
  i.e. containing only inequality constraints.
	\item (2) the probability distribution of \eqn{q} is sampled numerically using a random walk (based on the Metropolis algorithm). 
  There are three algorithms for selecting new samples: \code{rda}, \code{cda} (two hit-and-run algorithms) and a novel \code{mirror} algoithm.
	\item In the \code{rda} algorithm first a random direction is selected, and the new sample obtained by uniformly sampling the line 
  connecting the old sample and the intersection with the planes defined by the inequality constraints.
  \item the \code{cda} algorithm is similar, except that the direction is chosen along one of the coordinate axes.  
  \item the \code{mirror} algorithm is unpublished; it uses the inequality constraints as "reflecting planes" along which jumps are reflected.
\cr In contrast to \code{cda} and \code{rda}, this algorithm also works with unbounded problems (i.e. for which some of the unknowns can attain Inf).

For more information, see \file{xsamplevignette.pdf} in the packages \file{docs} subdirectory. 
  
  \item Raftery and Lewis (1996) suggest a minimum of 3000 iterations to reach the extremes 
  \item If provided, then \code{p} should be a valid particular solution (i.e. Ex=b and G>=h), else the algorithm will fail.
  \item If the particular solution (\code{x0}) is not provided, then the parsimonious solution is sought, see \code{\link{ldei}}. 
  For underdetermined systems, this may not be the most efficient way to start the algorithm}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

