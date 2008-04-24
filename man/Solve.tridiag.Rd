\name{Solve.tridiag}
\alias{Solve.tridiag}
\title{Solution of tridiagonal system of linear equations}
\description{Solves the linear system of equations \deqn{Ax=B}
 \cr where A has to be square and \emph{tridiagonal}, i.e with nonzero elements only on, one band above, and one band below the diagonal.}
\usage{Solve.tridiag(diam1, dia, diap1, rhs=rep(0,times=length(dia)))}
\arguments{
  \item{diam1 }{a vector with (nonzero) elements below the diagonal}
  \item{dia }{a vector with (nonzero) elements on the diagonal}
  \item{diap1 }{a vector with (nonzero) elements above the diagonal}
  \item{rhs }{numeric vector containing right hand side}
}
\value{
  vector with the solution, x, of the tridiagonal system of equations Ax=B
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# create tridagonal system: bands on diagonal, above and below
nn   <- 50                          # nr rows and columns of A
aa   <- runif(nn)
bb   <- runif(nn)
cc   <- runif(nn)
# full matrix
A                       <- matrix(nrow=nn,ncol=nn,data=0)
diag(A)                 <- bb
A[cbind(1:(nn-1),2:nn)] <- cc[-nn]
A[cbind(2:nn,1:(nn-1))] <- aa[-1]
B<-runif(nn)

# solve as full matrix
solve(A,B)                           

# same,  now using tridiagonal algorithm
Solve.tridiag(aa[-1],bb,cc[-nn],B)}
\details{If the length of the vector *dia* is N, then the lengths of diam1 and diap1 should be equal to N-1}
\seealso{\code{\link{Solve.banded}}, the function to solve a banded system of linear equations.
\cr \code{\link{Solve}} the generalised inverse solution, 
\cr \code{\link{solve}} the R default}
\keyword{array}

