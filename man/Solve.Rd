\name{Solve}
\alias{Solve}
\title{
  Generalised inverse solution of Ax=B
}
\description{
  Generalised inverse solution of \deqn{Ax=B}
  
  Solve uses the Moore-Penrose generalized inverse of matrix A
  (function \code{\link[MASS]{ginv}} from package MASS).
 
  Note: \code{\link{solve}}, the \R default requires a square, positive
  definite A. Solve does not have this restriction.
}

\usage{
Solve (A, B=diag(nrow=nrow(A)), tol=sqrt(.Machine$double.eps))
}

\arguments{
  \item{A }{numeric matrix containing the coefficients of the equations
    \eqn{Ax=B}.
  }
  \item{B }{numeric matrix containing the right-hand sides of the equations;
    the default is the unity matrix, in which case the function will return
    the Moore-Penrose generalized inverse of matrix A.
  }
  \item{tol }{tolerance for selecting singular values.
  }
}

\value{
  a vector with the generalised inverse solution.
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}

\examples{
A <- matrix(nrow=4,ncol=3,data=c(1:8,6,8,10,12)) # col3 = col1+col2
B <- 0:3
X <- Solve(A,B)            # generalised inverse solution
A \%*\% X - B                # should be zero (except for roundoff)
(gA <- Solve(A))           # generalised inverse of A
}
\seealso{
  \code{\link[MASS]{ginv}} to estimate the Moore-Penrose generalized inverse
  of a matrix, in package MASS,
\cr \code{\link{solve}} the R default
}
\keyword{array}

