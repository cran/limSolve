\name{xranges}
\alias{xranges}
\title{Calculates ranges of inverse unknowns}
\description{Given the linear constraints
 \deqn{Ex=f}
 \deqn{Gx>=h}
 finds the minimum and maximum values of all elements of vector \eqn{x}
 \cr by successively minimising and maximising each \code{x}, using linear programming
 \cr Uses \code{lpSolve} - see note}
\usage{xranges(E=NULL, F=NULL, G=NULL, H=NULL)}
\arguments{
  \item{E }{numeric matrix containing the coefficients of the equaliies \eqn{Ex=F}}
  \item{F }{numeric vector containing the right-hand side of the equalities}
  \item{G }{numeric matrix containing the coefficients of the inequalities \eqn{Gx>=H}}
  \item{H }{numeric vector containing the right-hand side of the inequalities}
}
\value{
  a 2-column matrix with the minimum and maximum value of each \code{x}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Estimate the ranges in the Diet Composition of Mink
A <- rbind(Minkdiet$Prey,rep(1,7))
B <- c(Minkdiet$Mink,1)
(DC<-xranges(A,B))}
\seealso{\code{\link{varranges}}, 
\cr \code{\link[lpSolve]{lp}}: linear programming from package lpSolve}
\note{we have experienced that \code{lpSolve} sometimes fails (and \R stops) if it is frequently repeated, and for small problems.}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

