\name{linp}
\alias{linp}
\title{Linear Programming}
\description{Solves a linear programming problem,
 \deqn{\min(\sum {Cost_i.x_i})} subject to        
 \deqn{x_i>=0}
 \deqn{Ex=f}
 \deqn{Gx>=h}

 
 Note: this function is simply a wrapper around \code{\link[lpSolve]{lp}} from package lpSolve
 \cr This R-code may fail and terminate \R for very small problems that are repeated frequently...
 \cr Sometimes multiple solutions exist for the same problem.}
\usage{linp(E=NULL, F=NULL, G=NULL, H=NULL, Cost, verbose=TRUE, ...)}
\arguments{
  \item{E }{numeric matrix containing the coefficients of the equality constraints \eqn{Ex=F}; if the columns of E have a names attribute, they will be used to label the output}
  \item{F }{numeric vector containing the right-hand side of the equality constraints}
  \item{G }{numeric matrix containing the coefficients of the inequality constraints \eqn{Gx>=H}; if the columns of G have a names attribute, and the columns of E do not, they will be used to label the output}
  \item{H }{numeric vector containing the right-hand side of the inequality constraints}
  \item{Cost }{numeric vector containing the coefficients of the cost function; if the Cost a names attribute, and neither the columns of E nor G have a name, they will be used to label the output}
  \item{verbose }{logical to print error messages}
  \item{... }{extra arguments passed to R-function \code{lp}}
}
\value{
  a list containing:
  \item{X }{vector containing the solution of the linear programming problem.}
  \item{residualNorm }{scalar, the sum of absolute values of residuals of equalities and violated inequalities. Should be very small or zero for a feasible linear programming problem}
  \item{solutionNorm }{scalar, the value of the minimised \code{Cost} function, i.e. the value of \eqn{\sum {Cost_i.x_i}}}
  \item{IsError }{logical, \code{TRUE} if an error occurred}
  \item{type }{the string "linp", such that how the solution was obtained can be traced}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{#--------------------------------------------
# Linear programming problem 1, not feasible
#--------------------------------------------

# maximise x1 + 3*x2 
# subject to
#-x1 -x2    < -3
#-x1 + x2   <-1
# x1 + 2*x2 < 2
G    <- matrix(nrow=3,data=c(-1,-1,1, -1,1,2))
H    <- c(3,-1,2)
Cost <- c(-1,-3)
(L<-linp(E=NULL,F=NULL,Cost=Cost,G=G,H=H))
L$residualNorm

#--------------------------------------------
# Linear programming problem 2, feasible
#--------------------------------------------

# minimise x12 + 8*x13 + 9*x14 + 2*x23 + 7*x24 + 3*x34  
# subject to:
#-x12             + x23 + x24       = 0
#     - x13       - x23       + x34 = 0
# x12 + x13 + x14                   > 1
#             x14       + x24 + x34 < 1
A    <- matrix(nrow=2,byrow=TRUE,data=c(-1,0,0,1,1,0,
                                        0,-1,0,-1,0,1))
B    <- c(0,0)
G    <- matrix(nrow=2,byrow=TRUE,data=c(1,1,1,0,0,0,
                                        0,0,-1,0,-1,-1))
H    <- c(1,-1)
Cost <- c(1,8,9,2,7,3)
(L<-linp(E=A,F=B,Cost=Cost,G=G,H=H))
L$residualNorm}
\seealso{\code{\link{ldei}},  \code{\link{lsei}}, 
\cr \code{\link[lpSolve]{lp}} the original function from package lpSolve
\cr \code{\link{Blending}}}
\keyword{algebra}
\keyword{array}
\keyword{optimize}

