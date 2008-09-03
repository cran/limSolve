\name{xsample}
\alias{xsample}
\title{Randomly samples an under- or overdetermined problem with linear equality and inequality conditions}
\description{Bayesian samping of linear problems with linear inequality conditions. 
 Uses either the "hit and run" (or directions algorithm) or a mirroring technique for sampling 
 
 The monte carlo method produces a sample solution for 
 \deqn{Ex=f}
 \deqn{Ax\simeq B}
 \deqn{Gx>=h}
 where Ex=F have to be met exactly, and \eqn{Ax\simeq B} approximately.}
\usage{xsample(A=NULL, B=NULL, E=NULL, F=NULL, 
G=NULL, H=NULL, sdB=1, iter=3000, outputlength=iter,
type = "mirror", jmp = NULL, tol=sqrt(.Machine$double.eps),
x0=NULL, fulloutput = FALSE, test=TRUE)}
\arguments{
  \item{A }{numeric matrix containing the coefficients of the
    (approximate) equality constraints, \eqn{Ax\simeq B}}
  \item{B }{numeric vector containing the right-hand side of the (approximate) equality constraints}
  \item{E }{numeric matrix containing the coefficients of the (exact) equality constraints, \eqn{Ex=F}}
  \item{F }{numeric vector containing the right-hand side of the (exact) equality constraints}
  \item{G }{numeric matrix containing the coefficients of the inequality constraints, \eqn{Gx>=H}}
  \item{H }{numeric vector containing the right-hand side of the inequality constraints}
  \item{sdB }{vector with standard deviation on B, used for weighing; default = equal weights}
  \item{iter }{integer determining the number of iterations}
  \item{outputlength} { number of iterations kept in the output; at most equal to \code{iter}}
  \item{type }{type of algorithm: one of: "mirror", (mirroring algorithm), "rda" (random directions algorithm) or "cda" (coordinates directions algorithm)}
  \item{jmp }{jump length of the transformed variables q: \eqn{x=x0+Zq} (only if \code{type}=="mirror"); if jmp is NULL, a reasonable value is determined by xsample,
  depending on the size of the NULL space.}
  \item{tol }{tolerance for equality and inequality constraints; numbers whose absolute value is smaller than \code{tol} are set to zero}
  \item{x0 }{initial (particular) solution}
  \item{fulloutput }{if TRUE, also outputs the transformed variables q}
  \item{test }{if TRUE, xsample will test for hidden equalities (see details). This may be necessary for large problems, but slows down execution a bit}

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

pairs(xsample(E=E,F=F,G=diag(7),H=rep(0,7),iter=5000,output=1000,type="cda")$X,
     main="Minkdiet 1000 solutions, - cda")

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
  
# equations that must be reproduced as closely as possible, Ax ~ b
# the field data; the input ratio matrix and field data are rescaled
A     <- t(Chemtax$Ratio/rowSums(Chemtax$Ratio))
B     <- Chemtax$Field/sum(Chemtax$Field)
  
# Sample
xs <- xsample(A=A,B=B,E=E,F=F,G=G,H=H,iter=2000,outputlen=1000)$X
pairs(xs, main="Chemtax 1000 solutions - 2000 iterations")

     }
\details{The algorithm proceeds in two steps.
	\item (1) the equality constraints \eqn{Ex=F} are eliminated, and the system \eqn{Ex=f}, \eqn{Gx>=h} is rewritten as \eqn{G(p+Zq)>= h}, 
  i.e. containing only inequality constraints and where Z is a basis for the null space of E.
	\item (2) the probability distribution of \eqn{q} is sampled numerically using a random walk (based on the Metropolis algorithm). 
  There are three algorithms for selecting new samples: \code{rda}, \code{cda} (two hit-and-run algorithms) and a novel \code{mirror} algoithm.
	\item In the \code{rda} algorithm first a random direction is selected, and the new sample obtained by uniformly sampling the line 
  connecting the old sample and the intersection with the planes defined by the inequality constraints.
  \item the \code{cda} algorithm is similar, except that the direction is chosen along one of the coordinate axes.  
  \item the \code{mirror} algorithm is unpublished; it uses the inequality constraints as "reflecting planes" along which jumps are reflected.
\cr In contrast to \code{cda} and \code{rda}, this algorithm also works with unbounded problems (i.e. for which some of the unknowns can attain Inf).

For more information, see the package vignette \code{vignette(xsample)} or the file xsample.pdf in the packages \file{docs} subdirectory.
  
Raftery and Lewis (1996) suggest a minimum of 3000 iterations to reach the extremes

If provided, then \code{x0} should be a valid particular solution (i.e. \eqn{E*x0=b} and \eqn{G*x0>=h}), else the algorithm will fail.

For larger problems, a central solution may be necessary as a starting point for the \code{rda} and \code{cda} algorithms.
A good starting value is provided by the "central" value when running the function \code{\link{xranges}} with option \code{central} equal to TRUE.

If the particular solution (\code{x0}) is not provided, then the parsimonious solution is sought, see \code{\link{ldei}}.
\cr For underdetermined systems, this may not be the most efficient way to start the algorithm. The parsimonious solution is usually
  located near the edges, and the \code{rda} and \code{cda} algorithms may not get out of this corner.
 The \code{mirror} algorithm is insensitive to that. Here it may be even better to start in a corner (as this position will always never be reached by random sampling).

The algorithm will fail if there are hidden equalities. For instance, two inequalities may together impose an equality on an unknown,
or, inequalities may impose equalities on a linear combination of two or more unknowns.

In this case, the basis of the null space Z will be deficient. Therefore, \code{xsample} starts by checking if such hidden equalities exist.
If it is suspected that this is NOT the case, set \code{test} to FALSE. This will speed up execution slightly.

It is our experience that for small problems either the \code{rda} and \code{cda} algorithms are often more efficient.
For really large problems, the \code{mirror} algorithm is usually much more efficient; select a jump length (\code{jmp}) that ensures good random coverage, while
still keeping the number of reflections reasonable. If unsure about the size of jmp, the default will do.

See \code{\link{E_coli}} for an example where a relatively large problem is sampled.
  }
\keyword{algebra}
\keyword{array}
\keyword{optimize}

