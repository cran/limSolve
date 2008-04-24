\name{varranges}
\alias{varranges}
\title{Calculates ranges of inverse variables}
\description{Given the linear constraints

 \deqn{Ex=f}
 \deqn{Gx>=h}
 and a set of "variables" described by the linear equations 
 \deqn{Var = EqA.x+EqB}
 finds the minimum and maximum values of the variables
  by successively minimising and maximising each variable equation
 \cr uses linear programming function lpSolve (see note)}
\usage{varranges(E=NULL, F=NULL, G=NULL, H=NULL, EqA, EqB=NULL)}
\arguments{
  \item{E }{numeric matrix containing the coefficients of the equalities \eqn{Ex=F}}
  \item{F }{numeric vector containing the right-hand side of the equalities}
  \item{G }{numeric matrix containing the coefficients of the inequalities \eqn{Gx>=H}}
  \item{H }{numeric vector containing the right-hand side of the inequalities}
  \item{EqA }{numeric matrix containing the coefficients that define the variable equations}
  \item{EqB }{numeric vector containing the right-hand side of the variable equations}
}
\value{
  a 2-column matrix with the minimum and maximum value of each equation (variable)
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Ranges in the contribution of food 3+4+5 in the diet of Mink

A <- rbind(Minkdiet$Prey,rep(1,7))
B <- c(Minkdiet$Mink,1)
EqA   <- c(0,0,1,1,1,0,0)            # sum of food 3,4,5
(isoA<-varranges(A,B,EqA=EqA))        # ranges of part of food 3+4+5}
\seealso{\code{\link{xranges}}, 
\cr \code{\link[lpSolve]{lp}}: linear programming from package lpSolve}
\note{we have experienced that \code{lpSolve} sometimes fails (and \R stops) if it is frequently repeated, and for small problems.}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

