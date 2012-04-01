\name{ldei}
\alias{ldei}
\title{
  Weighted Least Distance Programming with equality and inequality constraints.
}

\description{
  Solves the following \bold{underdetermined} inverse problem:
  \deqn{\min(\sum {x_i}^2)} subject to
  \deqn{Ex=f}
  \deqn{Gx>=h}

  uses least distance programming subroutine ldp (FORTRAN) from Linpack

  The model has to be UNDERdetermined, i.e. the number of independent
  equations < number of unknowns.
}

\usage{
ldei(E, F, G = NULL, H = NULL,
     tol = sqrt(.Machine$double.eps), verbose = TRUE)
}

\arguments{
  \item{E }{numeric matrix containing the coefficients of the equality
    constraints \eqn{Ex=F}; if the columns of \code{E} have a names attribute,
    they will be used to label the output.
  }
  \item{F }{numeric vector containing the right-hand side of the
    equality constraints.
  }
  \item{G }{numeric matrix containing the coefficients of the inequality
    constraints \eqn{Gx>=H}; if the columns of \code{G} have a names attribute
    and the columns of \code{E} do not, they will be used to label the output.
  }
  \item{H }{numeric vector containing the right-hand side of the inequality
    constraints.
  }
  \item{tol }{tolerance (for singular value decomposition, equality and
    inequality constraints).
  }
  \item{verbose }{logical to print \code{ldei} error messages.
  }
}
\value{
  a list containing:
  \item{X }{vector containing the solution of the least distance with
    equalities and inequalities problem.
  }
  \item{unconstrained.solution }{vector containing the unconstrained solution
    of the least distance problem, i.e. ignoring \eqn{Gx>=h}.
  }
  \item{residualNorm }{scalar, the sum of absolute values of residuals of
    equalities and violated inequalities; should be zero or very small if
    the problem is feasible.
  }
  \item{solutionNorm }{scalar, the value of the quadratic function at the
    solution, i.e. the value of \eqn{\sum {w_i*x_i}^2}.
  }
  \item{IsError }{logical, \code{TRUE} if an error occurred.
  }
  \item{type }{the string "ldei", such that how the solution was obtained
    can be traced.
  }
}

\author{
  Karline Soetaert <karline.soetaert@nioz.nl>.
}

\examples{
#-------------------------------------------------------------------------------
# A simple problem
#-------------------------------------------------------------------------------
# minimise x1^2 + x2^2 + x3^2 + x4^2 + x5^2 + x6^2
# subject to:
#-x1            + x4 + x5      = 0
#    - x2       - x4      + x6 = 0
# x1 + x2 + x3                 > 1
#           x3       + x5 + x6 < 1
# xi > 0

E <- matrix(nrow = 2, byrow = TRUE, data = c(-1, 0, 0, 1, 1, 0,
                                              0,-1, 0, -1, 0, 1))
F <- c(0, 0)
G <- matrix(nrow = 2, byrow = TRUE, data = c(1, 1, 1, 0, 0, 0,
                                             0, 0, -1, 0, -1, -1))
H    <- c(1, -1)
ldei(E, F, G, H)

#-------------------------------------------------------------------------------
# parsimonious (simplest) solution of the mink diet problem
#-------------------------------------------------------------------------------
E <- rbind(Minkdiet$Prey, rep(1, 7))
F <- c(Minkdiet$Mink, 1)

parsimonious <- ldei(E, F, G = diag(7), H = rep(0, 7))
data.frame(food = colnames(Minkdiet$Prey),
           fraction = parsimonious$X)
dotchart(x = as.vector(parsimonious$X),
         labels = colnames(Minkdiet$Prey),
         main = "Diet composition of Mink extimated using ldei",
         xlab = "fraction")
         
}

\references{
  Lawson C.L.and Hanson R.J. 1974. Solving Least Squares Problems, Prentice-Hall

  Lawson C.L.and Hanson R.J. 1995. Solving Least Squares Problems.
   SIAM classics in applied mathematics, Philadelphia.  (reprint of book)
}

\seealso{
  \code{\link{Minkdiet}}, for a description of the Mink diet example.

  \code{\link{lsei}}, \code{\link{linp}}

  \code{\link{ldp}}
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

