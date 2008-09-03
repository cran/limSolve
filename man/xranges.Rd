\name{xranges}
\alias{xranges}
\title{Calculates ranges of inverse unknowns}
\description{Given the linear constraints
 \deqn{Ex=f}
 \deqn{Gx>=h}
 finds the minimum and maximum values of all elements of vector \eqn{x}
 \cr by successively minimising and maximising each \code{x}, using linear programming
 \cr uses linear programming function \code{lp} from package lpSolve (see note)}
\usage{xranges(E=NULL, F=NULL, G=NULL, H=NULL, ispos=FALSE, tol=1e-8,
       central = FALSE, full=FALSE)}
\arguments{
  \item{E }{numeric matrix containing the coefficients of the equalities \eqn{Ex=F}}
  \item{F }{numeric vector containing the right-hand side of the equalities}
  \item{G }{numeric matrix containing the coefficients of the inequalities \eqn{Gx>=H}}
  \item{H }{numeric vector containing the right-hand side of the inequalities}
  \item{ispos }{if TRUE, it is imposed that unknowns are positive quantities}
  \item{tol }{tolerance for equality and inequality constraints}
  \item{central }{if TRUE, the mean value of all range solutions is also outputted}
  \item{full }{if TRUE, all range solutions are also outputted}
  }
\value{
  a matrix with at least two columns:

  column 1 and 2: the minimum and maximum value of each \code{x}

  if \code{central} is TRUE: column 3 = the central value

  if \code{full} is TRUE: next columns contain all valid range solutions
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Estimate the ranges in the Diet Composition of Mink
A <- rbind(Minkdiet$Prey,rep(1,7))
B <- c(Minkdiet$Mink,1)
(DC<-xranges(A,B,ispos=TRUE))
(DC<-xranges(A,B,G=diag(7),H=rep(0,7)))
}
\details{
The ranges are estimated by successively minimising and maximising each unknown, and using
linear programming.

by default linear programming assumes that all unknowns are positive. If this is the case,
set \code{ispos} equal to TRUE; the function will execute faster then. If \code{ispos}
is FALSE, then a system double the size of the original system must be solved for.

\code{xranges} outputs only the minimum and maximum value of each flow unless:

\code{full} is TRUE. In this case, all the results of the successive minimisation and maximisation will be outputted,
i.e. for each linear programming application, not just the value of the unknown being optimised but also the corresponding values of the other unknowns will be outputted.

If \code{central} is TRUE, then the mean of all the results of the linear programming will be outputted.
This may be a good starting value for \code{\link{xsample}}

Note: the columns corresponding to the \code{central} value and the \code{full} results are
valide solutions of the equations \eqn{Ex=F} and \eqn{Gx>=H}. This is not the case for the
first two columns (with the minimal and maximal values).
}

\seealso{\code{\link{varranges}}, for range estimation of variables,
\cr \code{\link{xsample}}, to randomly sample the lsei problem
\cr \code{\link[lpSolve]{lp}}: linear programming from package lpSolve}
\note{we have experienced that \code{lpSolve} sometimes fails (and \R stops) if it is frequently repeated, and for small problems.}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

