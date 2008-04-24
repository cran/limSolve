\name{resolution}
\alias{resolution}
\title{Row and column resolution of a matrix}
\description{Given an input matrix or its singular value decomposition,
 \cr calculates the resolution of the equations (rows) and of the unknowns (columns)
 of the matrix}
\usage{resolution(s, tol=sqrt(.Machine$double.eps))}
\arguments{
  \item{s }{either matrix or its singular value decomposition}
  \item{tol }{tolerance for singular values}
}
\value{
  a list containing:
  \item{row }{resolution of the rows  (equations)}
  \item{col }{resolution of the columns (variables)}
  \item{nsolvable }{number of solvable unknowns - rank of the matrix}
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>
          \cr Dick van Oevelen<d.vanoevelen@nioo.knaw.nl>}
\examples{resolution(matrix(nr=3,runif(9)))  #3rows,3columns
resolution(matrix(nr=3,runif(12))) #3rows,4columns
resolution(matrix(nr=3,runif(6)))  #3rows,2columns
resolution(cbind(c(1,2,3),c(2,3,4),c(3,5,7))) # r3=r1+r2,c3=c1+c2}
\references{Several references, for instance:
\cr Menke, W., 1989. Geophysical Data Analysis: Discrete Inverse Theory. Revised edition.
International Geophysics Series. Academic Press, London.}
\seealso{\code{\link{svd}}, the singluar value decomposition}
\keyword{array}

