\name{Solve.banded}
\alias{Solve.banded}
\title{Solution of a banded system of linear equations}
\description{Solves the linear system of equations \deqn{Ax=B} by Gaussion elimination  
 \cr where A has to be square, and \emph{banded}, i.e. with the only nonzero elements in bands near the diagonal. 
 
 The matrix A is either inputted as a full square matrix or as the non-zero bands
 \cr uses subroutines dgbfa and dgbsl (FORTRAN) from the Linpack routines}
\usage{Solve.banded(abd, nup, nlow,rhs=rep(0,times=ncol(abd)), 
full=(nrow(abd)==ncol(abd)))}
\arguments{
  \item{abd }{either a matrix containing the (nonzero) bands, rotated row-wise (anti-clockwise) only or a full square matrix}
  \item{nup }{number of nonzero bands above the diagonal; ignored if full matrix is inputted}
  \item{nlow }{number of nonzero bands below the diagonal; ignored if full matrix is inputted}
  \item{rhs }{numeric vector containing right hand side}
  \item{full }{if true: full matrix is passed, if false:banded}
}
\value{
  vector with solution, x, of banded system of equations Ax=B
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# Generate a banded matrix of random numbers
nup  <- 2                         # nr nonzero bands above diagonal
ndwn <- 3                         # nr nonzero bands below diagonal
nn   <- 10                        # nr rows and columns of A
A <- matrix(nrow=nn,ncol=nn,data=runif(1:(nn*nn)))
A [row(A)<col(A)-nup | row(A)>col(A)+ndwn]<-0 
diag(A) <- 1                      # 1 on diagonal is easily recognised 
# right hand side
B<-runif(nrow(A))                
# solve it, using the default solver and banded (inputting full matrix)
Full <- solve(A,B)
Band1<- Solve.banded(A,nup,ndwn,B)

# create banded form of matrix A
Aext <- rbind(matrix(nc=ncol(A),nr=nup,0),
              A,matrix(nc=ncol(A),nr=ndwn,0))
abd  <- matrix(nrow=nup+ndwn+1,ncol=nn,
  data=Aext[(col(Aext))<=row(Aext)&col(Aext)>=row(Aext)-ndwn-nup])
# print both to screen
A
abd
# solve problem with banded version
Band2 <- Solve.banded(abd,nup,ndwn,B)
# compare 3 methods of solution
rbind(Full,Band1,Band2)}
\references{J.J. Dongarra, J.R. Bunch, C.B. Moler, G.W. Stewart, LINPACK Users' Guide, SIAM, 1979.}
\details{If the input matrix \code{abd} is square, it is assumed that the full, square A is inputted, unless \code{full} is false.
	\cr is \code{abd} is not square, then the number of columns denote the number of unknowns, while the number of rows equals the nonzero bands, i.e. nup+nlow+1}
\seealso{\code{\link{Solve.tridiag}} To solve a tridiagonal system of linear equations.
\cr \code{\link{Solve}} the generalised inverse solution, 
\cr \code{\link{solve}} the R default}
\keyword{array}

