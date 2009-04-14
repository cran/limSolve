\name{varranges}
\alias{varranges}
\title{
  Calculates ranges of inverse variables
}
\description{
  Given the linear constraints

  \deqn{Ex=f}
  \deqn{Gx>=h}
  and a set of "variables" described by the linear equations
  \deqn{Var = EqA.x+EqB}
  finds the minimum and maximum values of the variables
  by successively minimising and maximising each variable equation

  uses linear programming function \code{\link[lpSolve]{lp}} from package
  \code{lpSolve}.
}
\usage{
varranges(E=NULL, F=NULL, G=NULL, H=NULL, EqA, EqB=NULL, 
  ispos=FALSE, tol=1e-8)
}

\arguments{
  \item{E }{numeric matrix containing the coefficients of the equalities
    \eqn{Ex=F}.
  }
  \item{F }{numeric vector containing the right-hand side of the
    equalities.
  }
  \item{G }{numeric matrix containing the coefficients of the inequalities
    \eqn{Gx>=H}.
  }
  \item{H }{numeric vector containing the right-hand side of the
    inequalities.
  }
  \item{EqA }{numeric matrix containing the coefficients that define the
    variable equations.
  }
  \item{EqB }{numeric vector containing the right-hand side of the variable
    equations.
  }
  \item{ispos }{if \code{TRUE}, it is imposed that unknowns are positive
    quantities.
  }
  \item{tol }{tolerance for equality and inequality constraints.
  }
}
\value{
  a 2-column matrix with the minimum and maximum value of each equation
  (variable)
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{# Ranges in the contribution of food 3+4+5 in the diet of Mink

A <- rbind(Minkdiet$Prey,rep(1,7))
B <- c(Minkdiet$Mink,1)
EqA   <- c(0,0,1,1,1,0,0)                 # sum of food 3,4,5
(isoA<-varranges(A,B,EqA=EqA,ispos=TRUE)) # ranges of part of food 3+4+5
varranges(A,B,EqA=EqA,G=diag(7),H=rep(0,7)) 
}
\seealso{
  \code{\link{xranges}},

  \code{\link[lpSolve]{lp}}: linear programming function from package lpSolve.
}

\keyword{algebra}
\keyword{array}
\keyword{optimize}

