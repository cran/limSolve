
################################################################################
##                                                                            ##
## LINEAR INVERSE MODELLING    -   LIM                                        ## 
##    Karline Soetaert and Karel Van den Meersche                             ##
##                                                                            ##
##   routines that relate to the following problem:                           ##
##   1. min or max (f(x)), where f(x) = ||Ax-b||^2 or f(x) = sum(ai*xi)       ##  
##      subject to equality constraints Ex=f and inequality constraints Gx>=h ##
##   2. find the root of a nonlinear function                                 ##
##                                                                            ##
## part 1: solving Linear Inverse Models                                      ##
## -----------------------------                                              ##
## ldei        : Least distance with equalities and inequalities              ##
## ldp         : Least Distance Programming (subject to inequalities)         ##
## linp        : linear programming, uses lp from package lpSolve             ##
## lsei        : Least Squares with Equalities and Inequalities               ##
## resolution  : calculates the resolution of equations and variables         ##
## xranges     : Calculates ranges of inverse unknowns                        ##
## varranges   : Calculates ranges of inverse equations (variables)           ##
##                                                                            ##
## part 2: utilities for solving systems of linear equations                  ##
## -----------------------------                                              ##
## Solve         : Generalised inverse solution                               ##
## Solve.tridiag : Solves tridiagonal system of linear equations              ##
## Solve.banded  : Solves banded system of linear equations                   ##
##                                                                            ##
## part 3: utility for sampling underdeterminged system of linear equations   ##
## -----------------------------                                              ##
## xample  : coordinates directions or mirror algorithm for sampling LIM      ##
################################################################################


##############################################################################
##                                                                          ##
## LINEAR INVERSE MODELLING                                                 ## 
##    SOLVING LIM                                                           ##
##                                                                          ##
##############################################################################


################################################################################
## ldei        : Solves underdetermined problem, Least Distance Programming   ##
################################################################################

ldei <- function(E,      # numeric matrix containing the coefficients of the equality constraints Ex=F
                 F,       # numeric vector containing the right-hand side of the equality constraints
                 G=NULL,  # numeric matrix containing the coefficients of the inequality constraints G*X>=H
                 H=NULL,  # numeric vector containing the right-hand side of the inequality constraints
                 tol=sqrt(.Machine$double.eps),  # Tolerance (for singular value decomposition, equality and inequality constraints)
                 verbose=TRUE)                   # Logical to print ldei error message

{
  ## input consistency
  if (! is.matrix(E) & ! is.null(E)) E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G)) G <- t(as.matrix(G))

  ##------------------------------------------------------------------------
  ## 0. Setup problem
  ##------------------------------------------------------------------------
  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Neq    <- nrow(E)   # number of equations
  Nx     <- ncol(E)   # number of unknowns
  Nin    <- nrow(G)   # number of inequalities

  ## copies
  EE <- E
  GG <- G


  ## Consistency of input
  if (!is.null(G)) 
    {
      if (ncol(G)   != Nx)  stop("cannot solve least distance problem - E and G not compatible")
      if (length(H) != Nin) stop("cannot solve least distance problem - G and H not compatible")
    }

  if (length(F) != Neq) stop("cannot solve least distance problem - E and F not compatible")


  ##------------------------------------------------------------------------
  ## 1. Decompose matrix by singular value decomposition
  ##------------------------------------------------------------------------
  S          <- svd(E,nrow(E),ncol(E))

  ## number of 'solvable unknowns' - the rank of the problem: 
  solvable   <- sum(S$d > tol * S$d[1])     # inequality gives TRUE (1) or FALSE (0)
  if (solvable > Nx) stop("cannot solve problem by LDEI - overdetermined system")

  ## number of 'unsolvable' unknowns
  unsolvable <- Nx-solvable

  ##------------------------------------------------------------------------
  ## 2. Backsubstitution
  ##------------------------------------------------------------------------
  Xb         <- S$v[,1:solvable] %*% (t(S$u[,1:solvable])/S$d[1:solvable])%*%F

  ## Check if solution is correct
  CC         <- E %*%Xb - F
  if (any(abs(CC)>tol)) stop("cannot solve problem by LDEI - equations incompatible")

  ##------------------------------------------------------------------------
  ## 3. Constrained solution
  ##------------------------------------------------------------------------
  ## Check if inequalities are already satisfied

  ifelse (!is.null(G), CC        <- G%*%Xb-H, CC<-0)

  if (all(CC > -tol))     # done
    {
      X    <- Xb
      IsError <-FALSE
    } else   {

      ## Number of unknowns to be solved by the least distance programming 
      LDPNx   <-unsolvable

      ## Construct an orthogonal basis of V2 = (I-V*VT)*Random
      rnd     <- matrix(runif(Nx*unsolvable),nrow=Nx,ncol=unsolvable)
      V2      <- diag(1,Nx) - S$v[,1:solvable]%*%t(S$v[,1:solvable])
      V2      <- V2 %*% rnd

      ## Orthogonal basis, constructed by Singular Value Decomposition on V2. 
      Vort    <- svd(V2,nu=0,nv=unsolvable)  ## only right-hand side needed
      ortho   <- V2 %*% Vort$v[,1:unsolvable]
      for (i in 1:Nx) ortho[i,]<- ortho[i,]/Vort$d[1:unsolvable]

      ## Least distance programming matrices LDPG = G*ortho, LDPH = H - G*Xb 
      LDPG    <- G %*% ortho
      LDPH    <- H - G %*% Xb

      ## call the DLL with the least distance programming routine ldp
      IsError <- FALSE
      NW      <- (LDPNx+1)*(Nin+2) +2*Nin
      sol  <-.Fortran("ldp",G=LDPG,H=LDPH,
                      NUnknowns=as.integer(LDPNx),NConstraints=as.integer(Nin),
                      NW=as.integer(NW),X=as.vector(rep(0,LDPNx)),XNorm=0.,
                      W=as.double(rep(0.,NW)),xIndex=as.integer(rep(0,Nin)),
                      Mode=as.integer(0),
                      verbose=as.logical(verbose),IsError=as.logical(IsError))

      IsError<-sol$IsError 

      ## The solution, corrected for weights
      X    <- ortho[,1:unsolvable] %*% as.matrix(sol$X) + Xb
      X[which(abs(X)<tol)] <- 0         ## zero very tiny values
      
    } ## end CC > - tol

  ## Residual of the inequalities
  residin  <- 0
  if (!is.null(GG)) {
    ineq     <- GG %*% X - H
    residin  <- -sum(ineq[ineq<0])
  }
  ## Total residual norm
  residual <- sum(abs(EE %*% X - F))+residin

  ## The solution norm
  solution <- sum ((X)^2)
  X    <- as.vector(X)
  xnames <- colnames(E)
  if (is.null(xnames)) xnames <- colnames(G)
  names (X) <- xnames

  return(list(X=X,                            # vector containing the solution of the least distance problem.
              unconstrained.solution=as.vector(Xb),  # vector containing the unconstrained solution of the least distance problem
              residualNorm=residual,          # scalar, the sum of residuals of equalities and violated inequalities
              solutionNorm=solution,          # scalar, the value of the quadratic function at the solution
              IsError=IsError,                # if an error occurred
              type="ldei"))

} ########## END OF ldei ########## 

################################################################################
## nnls        : Solves nonnegative least squares problem                     ##
################################################################################

nnls <- function(A,  # numeric matrix containing the coefficients of the equality constraints A*X~=B
                 B,  # numeric vector containing the right-hand side of the equality constraints
                 tol=sqrt(.Machine$double.eps),   # Tolerance (for singular value decomposition, equality and inequality constraints)
                 verbose=TRUE)                   # Logical to print ldp error message

{

  ##------------------------------------------------------------------------
  ## 0. Setup problem
  ##------------------------------------------------------------------------
  ## input consistency
  if (! is.matrix(A) & ! is.null(A)) A <- t(as.matrix(A))

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Nx     <- ncol(A)   # number of unknowns
  Neq    <- nrow(A)   # number of inequalities
  if (length(B) != Neq) stop("cannot solve nnls problem - A and B not compatible")

  sol  <-.Fortran("xnnls",A=A,MDA=as.integer(Neq),M=as.integer(Neq),
                  N=as.integer(Nx),B=B,X=as.vector(rep(0,Nx)),RNorm=0.,
                  W=as.double(rep(0.,Nx)),ZZ=as.double(rep(0.,Neq)),
                  Index=as.integer(rep(0,Nx)), Mode=as.integer(0))
  IsError <- FALSE
  Mode <- sol$Mode

  if (Mode != 1) 
    {IsError <- TRUE
     if (verbose) {
       if (Mode ==2) warning("nnls: The dimensions of the problem are bad")
       if (Mode ==3) warning("nnls: iteration count exceeded: more than 3*N iterations")
     }
   }
  
  ## The solution
  X    <- sol$X

  ## Residual of the inequalities
  residual  <- -sum(X[X<0])

  ## The solution norm
  solution <- sum(abs(A %*% X - B))

  xnames <- colnames(A)
  names (X) <- xnames

  return(list(X=X,                            # vector containing the solution of the least distance problem.
              residualNorm=residual,          # scalar, the sum of violated inequalities (i.e where x<0)
              solutionNorm=solution,          # scalar, the value of the quadratic function at the solution
              IsError=IsError,                # if an error occurred
              type="nnls"))

} ########## END OF nnls ########## 

################################################################################
## ldp         : Solves Least Distance Programming                            ##
################################################################################

ldp <- function(G,  # numeric matrix containing the coefficients of the inequality constraints G*X>=H
                H,  # numeric vector containing the right-hand side of the inequality constraints
                tol=sqrt(.Machine$double.eps),   # Tolerance (for singular value decomposition, equality and inequality constraints)
                verbose=TRUE)                   # Logical to print ldp error message

{

  ##------------------------------------------------------------------------
  ## 0. Setup problem
  ##------------------------------------------------------------------------
  ## input consistency
  if (! is.matrix(G) & ! is.null(G)) G <- t(as.matrix(G))

  if (is.null(tol)) tol <- sqrt(.Machine$double.eps)

  ## Problem dimension
  Nx     <- ncol(G)   # number of unknowns
  Nin    <- nrow(G)   # number of inequalities
  if (length(H) != Nin) stop("cannot solve least distance problem - G and H not compatible")

  IsError <- FALSE
  NW      <- (Nx+1)*(Nin+2) +2*Nin
  sol  <-.Fortran("ldp",G=G,H=H,
                  NUnknowns=as.integer(Nx),NConstraints=as.integer(Nin),
                  NW=as.integer(NW),X=as.vector(rep(0,Nx)),XNorm=0.,
                  W=as.double(rep(0.,NW)),xIndex=as.integer(rep(0,Nin)),
                  Mode=as.integer(0),
                  verbose=as.logical(verbose),IsError=as.logical(IsError))
  IsError<-sol$IsError 
  
  ## The solution
  X    <- sol$X

  ## Residual of the inequalities
  residual  <- 0
  if (!is.null(G)) {
    ineq     <- G %*% X - H
    residual  <- -sum(ineq[ineq<0])
  }
  ## The solution norm
  solution <- sum ((X)^2)

  xnames <- colnames(G)
  names (X) <- xnames

  return(list(X=X,                            # vector containing the solution of the least distance problem.
              residualNorm=residual,          # scalar, the sum of residuals of violated inequalities
              solutionNorm=solution,          # scalar, the value of the quadratic function at the solution
              IsError=IsError,                # if an error occurred
              type="ldp"))

} ########## END OF ldp ########## 

################################################################################
## linp        : linear programming, uses lp from package lpSolve             ##
################################################################################

linp <- function(E=NULL, # numeric matrix containing the coefficients of the equality constraints Ex=F
                 F=NULL, # numeric vector containing the right-hand side of the equality constraints        
                 G=NULL, # numeric matrix containing the coefficients of the inequality constraints G*X>=H
                 H=NULL, # numeric vector containing the right-hand side of the inequality constraints
                 Cost,   # numeric vector containing the coefficients of the cost function
                 ispos=TRUE,    #
                 int.vec=NULL,  # Numeric vector giving the indices of variables that are required to be integer. The length of this vector will therefore be the number of integer variables.
                 verbose=TRUE,  # Logical to print error message
                 ...)           # extra arguments passed to R-function lp
  
  ##------------------------------------------------------------------------
  ## Solves a linear programming problem, 
  ## Minimise               Cost 
  ## subject to             E*x=f 
  ##                        G*x>=h
  ##
  ## Note: uses lp from package lpSolve 
  ## This R-code sometimes fails and terminates R
  ##      for very small problems that are repeated frequently...
  ##------------------------------------------------------------------------

{
  ## input consistency
  if (! is.matrix(E) & ! is.null(E)) E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G)) G <- t(as.matrix(G))


  ## problem dimension
  Neq  <- nrow(E)   # Number equalities
  Nx   <- ncol(E)   # Number unknowns
  Nin  <- nrow(G)   # Number inequalities

  if (is.null(Nx )) Nx  <- ncol(G)
  if (is.null(Nin)) Nin <- 0
  if (is.null(Neq)) Neq <- 0

  NumEq <- Nin+Neq  # total number of equations

  ## consistency of input
  if (!is.null(G)) 
    {
      if (ncol(G)   != Nx)   stop("cannot solve linear programming problem - E and G not compatible")
      if (length(H) != Nin)  stop("cannot solve linear programming problem - G and H not compatible")
    }

  if (!is.null(E)) 
    {
      if (length(F) != Neq)  stop("cannot solve linear programming problem - E and F not compatible")
    }

  if (length(Cost)!= Nx) stop("cannot solve linear programming problem - Cost not compatible")
  IsError <- FALSE

  ## con: constraints ; rhs: right hand side

  ## the equalities:
  con   <- E
  rhs   <- F
  dir   <- rep("==",Neq)

  ## inequalities:
  if (Nin > 0)
    {
      con   <- rbind(con,G)
      rhs   <- c(rhs,H)
      dir   <- c(dir,rep(">=",Nin))
    }

  if (!ispos)
    {
      con  <- cbind(con, -con)
      Cost <- c(Cost, -Cost)
      if (! is.null(int.vec)) int.vec<-c(int.vec,int.vec+Nx)
    }

  ## the solution
  sol    <- lp("min",Cost,con,dir,rhs,int.vec=int.vec,...)
  mode   <- sol$status
  ## from the original code lp_lib.h
  if (mode == -5) {print("unknown error");IsError<-TRUE}
  if (mode == -4) {print("data ignored");IsError<-TRUE}
  if (mode == -3) {print("no bfp");IsError<-TRUE}
  if (mode == -2) {print("no memory");IsError<-TRUE}
  if (mode == -1) {print("problem not run");IsError<-TRUE}
  if (mode == 2) {print("problem infeasible");IsError<-TRUE}
  if (mode == 3) {print("problem unbounded");IsError<-TRUE}
  if (mode == 4) {print("problem degenerate");IsError<-TRUE}
  if (mode == 5) {print("problem failed");IsError<-TRUE}
  if (mode == 6) {print("user aborted");IsError<-TRUE}
  if (mode == 7) {print("timed out");IsError<-TRUE}
  if (mode == 8) {print("running");IsError<-TRUE}
  if (mode == 9) {print("presolved");IsError<-TRUE}

  X           <- sol$solution
  if (!ispos) X <- X[1:Nx]-X[(Nx+1):(2*Nx)]
  solutionNorm<- sol$objval

  ## Total residual norm
  residual <- 0
  if (!is.null(E))  residual <- sum(abs(E %*% X - F))+residual 
  if (!is.null(G))  
    {
      ineq     <- G %*% X - H
      residual <- residual -sum(ineq[ineq<0])
    }
  xnames <- colnames(E)
  if (is.null(xnames)) xnames <- colnames(G)
  if (is.null(xnames)) xnames <- names(Cost)
  names (X) <- xnames

  return(list(X=X,                        # vector containing the solution of the linear programming problem.
              residualNorm=residual,      # scalar, the sum of residuals of equalities and violated inequalities
              solutionNorm=solutionNorm,  # scalar, the value of the minimised cost function
              IsError=IsError,            # if an error occurred
              type="linp"))

} ########## END OF linp ##########                 


################################################################################
## lsei        : Least Squares with Equalities and Inequalities               ##
################################################################################

lsei <- function(A=NULL,                     # numeric matrix containing the coefficients of the quadratic function to be minimised, ||Ax-B||
                 B=NULL,                     # numeric vector containing the right-hand side of the quadratic function to be minimised
                 E=NULL,                     # numeric matrix containing the coefficients of the equality constraints, Ex=F
                 F=NULL,                     # numeric vector containing the right-hand side of the equality constraints
                 G=NULL,                     # numeric matrix containing the coefficients of the inequality constraints, Gx>=H
                 H=NULL,                     # numeric vector containing the right-hand side of the inequality constraints
                 Wx=NULL,                    # x-Weighting coefficients 
                 Wa=NULL,                    # weights of the quadratic function (A, b matrix) to be minimised
                 type = 1,                        # integer code determining algorithm to use 1=lsei , 2=solve.QP from R-package quadprog
                 tol=sqrt(.Machine$double.eps),   # Tolerance (for singular value decomposition, equality and inequality constraints)
                 tolrank=NULL,   # Tolerance for equations
                 fulloutput = FALSE, 
                 verbose=TRUE)                    # Logical to print error messages

  ##------------------------------------------------------------------------
  ## Solves an lsei inverse problem, 
  ## lsei = Least Squares with Equality and Inequality constraints 
  ## Minimise             ||A*X-B|| 
  ## subject to             E*X=F 
  ##                        G*X>=H
  ## uses either the LINPACK package lsei or solve.QP from package quadprog
  ##------------------------------------------------------------------------

{

  ##------------------------------------------------------------------------
  ## Setup problem
  ##------------------------------------------------------------------------
  ## input consistency

  if (is.vector(E) & length(F)==1) E <- t(E) else if (! is.matrix(E) & ! is.null(E)) E <- as.matrix(E)
  if (is.vector(A) & length(B)==1) A <- t(A) else if (! is.matrix(A) & ! is.null(A)) A <- as.matrix(A)
  if (is.vector(G) & length(H)==1) G <- t(G) else if (! is.matrix(G) & ! is.null(G)) G <- as.matrix(G)
  if (! is.matrix(F) & ! is.null(F)) F <- as.matrix(F)
  if (! is.matrix(B) & ! is.null(B)) B <- as.matrix(B)
  if (! is.matrix(H) & ! is.null(H)) H <- as.matrix(H)
  if (is.null(A) && is.null (E)) {
    if(is.null(G))   stop("cannot solve least squares problem - A, E AND G are NULL")
    A <- matrix(nrow=1,ncol=ncol(G),0)
    B <- 0                       }  else if (is.null(A)) {
      A <- matrix(nrow=1,data=E[1,]);B<-F[1]}

  ## Problem dimension
  Neq  <- nrow(E)
  Napp <- nrow(A)
  Nx   <- ncol(A)
  Nin  <- nrow(G)
  if (is.null (Nx))   Nx  <- ncol(E)
  if (is.null (Nx))   Nx  <- ncol(G)
  if (is.null (Neq))  Neq  <- 1
  if (is.null (Nin))  Nin  <- 1
  ## If equalities/inequalities absent: "empty" constraints...
  if (is.null(E)) E <- matrix(nrow=1,ncol=Nx,0)
  if (is.null(F)) F <- 0
  if (is.null(G)) G <- matrix(nrow=1,ncol=Nx,0)
  if (is.null(H)) H <- 0

  if (ncol(G)   != Nx)   stop("cannot solve least squares problem - A and G not compatible")
  if (ncol(E)   != Nx)   stop("cannot solve least squares problem - A and E not compatible")
  if (length(B) != Napp) stop("cannot solve least squares problem - A and B not compatible")
  if (length(F) != Neq)  stop("cannot solve least squares problem - E and F not compatible")
  if (length(H) != Nin)  stop("cannot solve least squares problem - G and H not compatible")

  if (! is.null(Wa))
    {
      if (length(Wa) != Napp) stop ("cannot solve least squares problem - Wa should have length = number of rows of A") 

      A <- A*Wa
      B <- B*Wa
    }

  Tol <- tol 
  if (is.null(Tol)) Tol <- sqrt(.Machine$double.eps)

  ##------------------------------------------------------------------------
  ## Solution
  ##------------------------------------------------------------------------

  IsError <- FALSE

  if (type == 1)     # use LINPACKs lsei
    {
      ineq <- Nin+Nx
      mIP  <- ineq+2*Nx+2

      ## extra options?
      lpr <- 1
      if (fulloutput) lpr <- lpr+3
      if (! is.null(tolrank)) lpr <- lpr + 6
      if (! is.null(Wx))   {lw  <- length (Wx) ; lpr <- lpr + 2 + lw}
      
      ProgOpt <- rep(1.0,lpr)
      if (lpr>1)
        {
          ipr <- 1
          if (fulloutput) {ProgOpt[ipr:(ipr+2)] <- c(ipr+3,1,1); ipr <- ipr+3}
          if (! is.null(tolrank))
            {
              if (length(tolrank) == 1) tolrank <- rep(tolrank,len=2)
              ProgOpt[ipr:(ipr+5)] <- c(ipr+6,4,tolrank[1],ipr+6,5,tolrank[2]) 
              ipr <- ipr+6
            }
          if (! is.null(Wx)) 
            {
              lw  <- length (Wx) 
              if (lw ==1) {ProgOpt[ipr:(ipr+2)]<- c(ipr+3,2,1)} else
              {if (lw != Nx) stop("cannot solve least squares problem - number of weighs should be =1 or =number of unknowns") 
               lw <- lw + ipr+1
               ProgOpt[ipr:lw] <- c(lw+1,3,Wx)
             }
            }
        } 
      mdW <- Neq + Napp + ineq
      if (fulloutput) mdW <- max(mdW, Nx)
      mWS <- 2*(Neq+Nx)+max(Napp+ineq,Nx)+(ineq+2)*(Nx+7)

      sol <-.Fortran("lsei",NUnknowns=Nx,NEquations=Neq,NConstraints=Nin,NApproximate=Napp, 
                     A=A,B=B,E=E,F=F,G=G,H=H,X=as.vector(rep(0,Nx)),
                     mIP=as.integer(mIP),mdW=as.integer(mdW),mWS=as.integer(mWS),
                     IP=as.integer(rep(0,mIP)),
                     W=as.double(matrix(nrow=mdW,ncol=Nx+1,0.)),WS=as.double(rep(0.,mWS)),
                     lpr=as.integer(lpr),ProgOpt=as.double(ProgOpt),
                     verbose=as.logical(verbose),IsError=as.logical(IsError))

      if (fulloutput) 
        {covar<-matrix(nrow=mdW,ncol=Nx+1,data=sol$W)[1:Nx,1:Nx]
         RankEq <- sol$IP[1]
         RankApp <- sol$IP[2]
       }
    } else if (type == 2) 


      {                 ## use R's solve.QP, package quadprog
        if (! is.null(Wx)) stop ("cannot solve least squares problem - weights not implemented for type 2") 
        if (! is.null(Wa)) stop ("cannot solve least squares problem - weights not implemented for type 2") 

        dvec  <- t(A) %*% B
        Dmat  <- t(A) %*% A
        diag(Dmat) <- diag(Dmat)+1e-8
        Amat  <- t(rbind(E,G))
        bvec  <- c(F,H)

        sol   <- solve.QP(Dmat ,dvec, Amat , bvec, meq=Neq)
        sol$IsError <- FALSE

        if (sol$value>Tol) 
          {
            sol $IsError<- TRUE
            print("problem incompatible")
          }
        sol$X <- sol$solution

      } else  stop("cannot solve least squares problem - type unknown")



  X <- sol$X
  X[which(abs(X)<Tol)] <- 0         # zero very tiny values

  ## Residual of the inequalities
  residual <- 0
  if (Nin > 0) 
    { 
      ineq     <- G %*% X - H
      residual  <- residual -sum(ineq[ineq<0])
    }

  ## Total residual norm
  if (Neq> 0) 
    {
      residual <- residual + sum(abs(E %*% X - F))  
    }
  if (residual>Tol) 
    {
      sol$IsError<- TRUE
    }

  ## The solution norm
  
  solution <- 0
  if (Napp > 0)
    {
      solution <- sum ((A %*% X - B)^2)
    }

  xnames <- colnames(A)
  if (is.null(xnames)) xnames <- colnames(E)
  if (is.null(xnames)) xnames <- colnames(G)
  names (X) <- xnames

  res <- list(X=X,                            # vector containing the solution of the least squares problem.
              residualNorm=residual,          # scalar, the sum of residuals of equalities and violated inequalities
              solutionNorm=solution,          # scalar, the value of the minimised quadratic function at the solution
              IsError=sol$IsError,            # if an error occurred
              type="lsei")

  if (fulloutput && type == 1)
    {
      res$covar<-covar
      res$RankEq <- sol$IP[1]
      res$RankApp <- sol$IP[2]
    }

  return(res)

} ########## END OF lsei ##########                 



################################################################################
## resolution  : calculates the resolution of equations and variables         ##
################################################################################

resolution <- function (s,              # either matrix or its singular value decomposition
                        tol=sqrt(.Machine$double.eps))  # tolerance for singular values

  ##------------------------------------------------------------------------
  ## Given the singular value decomposition s, or the input matrix  
  ## calculates the resolution of the equations (rows) and of the variables (columns)
  ## of the matrix
  ##------------------------------------------------------------------------
{

  if (is.numeric(s)) s <- svd(s)
  solvable <- sum (s$d>tol*s$d[1])       # number of sufficiently large singular values

  ## The resolution of the equations
  resolutioneq <- diag(s$u[,1:solvable]%*%t(s$u[,1:solvable]))

  ## The resolution of the variables
  resolutionvar <- diag(s$v[,1:solvable]%*%t(s$v[,1:solvable]))

  return(list(row=resolutioneq,      # resolution of the rows  (equations)
              col=resolutionvar,     # resolution of the columns (variables)
              nsolvable=solvable))   # number of solvable unknowns - rank

}


xranges  <-  function (E = NULL, F = NULL, G = NULL, H = NULL, 
                       ispos=FALSE, tol = 1e-8, central = FALSE, full=FALSE)
  ##------------------------------------------------------------------------
  ## Given the linear constraints
  ##                        E*X=F 
  ##                        G*X>=H
  ##
  ## finds the minimum and maximum values of all elements of vector X 
  ## by successively minimising and maximising each x, using linear programming
  ## uses lpSolve - may fail (if frequently repeated)                       
  ## unknowns can possibly be negative unless ispos=TRUE
  ## if all are positive, then it is solved much faster.
  ##------------------------------------------------------------------------

{
  ## input consistency
  if (!is.matrix(E) & !is.null(E))
    E <- t(as.matrix(E))
  if (!is.matrix(G) & !is.null(G))
    G <- t(as.matrix(G))
  ## Dimensions of the problem
  Neq <- nrow(E)
  Nx  <- ncol(E)
  if (is.null(Nx)) Nx <- ncol(G)
  Nineq <- nrow(G)
  if (is.null(Nineq))
    Nineq <- 0
  if (is.null(Neq))
    Neq <- 0
  Range <- matrix(ncol = 2, nrow = Nx, NA)

  ## con: constraints ; rhs: right hand side
  ## First the equalities 
  con <- E
  rhs <- F
  dir <- rep("==", Neq)
  ## then the inequalities
  if (Nineq > 0) {
    con <- rbind(con, G)
    rhs <- c(rhs, H)
    dir <- c(dir, rep(">=", Nineq))
  }
  
  AllX   <- NULL
  Summed <- rep(0,Nx)
  nsum   <- 0

  if (ispos) {

    for (i in 1:Nx) {
      obj <- rep(0, Nx)
      obj[i] <- 1
      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0) Range[i, 1] <- lmin$objval else
      if(lmin$status == 3) Range[i, 1] <- -1e30 else
                           Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0) Range[i, 2] <- lmax$objval else
      if(lmax$status == 3) Range[i, 2] <- 1e30 else
                           Range[i, 2] <- NA
      if (central)
        {
          if (! any (is.na(lmin$solution)) && lmin$status==0) {Summed<- Summed + lmin$solution;nsum<-nsum+1}
          if (! any (is.na(lmax$solution)) && lmax$status==0) {Summed<- Summed + lmax$solution;nsum<-nsum+1}
        }
      if (full) {
        if( lmin$status==0)     AllX <- cbind(AllX,lmin$solution)
        if( lmax$status==0)     AllX <- cbind(AllX,lmax$solution)
      }
    }
  } else{
    ## First test if problem is solvable...
    Sol <- lsei(E=E,F=F,G=G,H=H)
    if (Sol$residualNorm > tol)
      { 
        warning (paste("cannot proceed: problem not solvable at requested tolerance",tol))
        return(Range)
      } 

    ## double the number of unknowns: x -> x1 -x2, where x1>0 and x2>0
    con <- cbind(con,-1*con)

    for (i in 1:Nx) {
      obj <- rep(0, 2*Nx)
      obj[i]    <- 1
      obj[Nx+i] <- -1

      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0) Range[i, 1] <- lmin$objval else
      if(lmin$status == 3) Range[i, 1] <- -1e30 else
                           Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0) Range[i, 2] <- lmax$objval else
      if(lmax$status == 3) Range[i, 2] <- 1e30 else
                           Range[i, 2] <- NA
      if (central)
        {
          if (! any (is.na(lmin$solution)) && lmin$status==0) {Summed<- Summed + lmin$solution;nsum<-nsum+1}
          if (! any (is.na(lmax$solution)) && lmax$status==0) {Summed<- Summed + lmax$solution;nsum<-nsum+1}
        }
      if (full) {
        if( lmin$status==0)     AllX <- cbind(AllX,lmin$solution)
        if( lmax$status==0)     AllX <- cbind(AllX,lmax$solution)
      }
    }
    if (central) Summed <- Summed[1:Nx]- Summed[(Nx+1):(2*Nx)]
    if (full)    AllX   <- AllX[1:Nx,]-AllX[(Nx+1):(2*Nx),]
  }
  colnames(Range) <- c("min", "max")
  xnames <- colnames(E)
  if (is.null(xnames))
    xnames <- colnames(G)
  rownames(Range) <- xnames
  if (central) Range<-cbind(Range, central = Summed/nsum)
  if (full) Range<-cbind(Range, AllX)
  return(Range)
}


################################################################################
## varranges    : Calculates ranges of inverse equations (variables)          ##
################################################################################

varranges <- function(E=NULL, # numeric matrix containing the coefficients of the equaliies Ex=F
                      F=NULL, # numeric vector containing the right-hand side of the equalities
                      G=NULL, # numeric matrix containing the coefficients of the inequalities G*X>=H
                      H=NULL, # numeric vector containing the right-hand side of the inequalities
                      EqA,    # numeric matrix containing the coefficients that define the variables
                      EqB=NULL, # numeric vector containing the right-hand side of the variable equation
                      ispos=FALSE,
                      tol=1e-8) 

  ##------------------------------------------------------------------------
  ## Given the linear constraints
  ##                        E*X=F 
  ##                        G*X>=H
  ## and a set of variables described by the linear equations Var = EqA*X+EqB
  ##
  ## finds the minimum and maximum values of the variables
  ## by successively minimising and maximising each linear combination, 
  ## using linear programming
  ## uses lpSolve - may fail (if frequently repeated)                       
  ##------------------------------------------------------------------------

{
  ## input consistency
  if (! is.matrix(E) & ! is.null(E)) E <- t(as.matrix(E))
  if (! is.matrix(G) & ! is.null(G)) G <- t(as.matrix(G))
  if (! is.matrix(EqA) & ! is.null(EqA)) EqA <- t(as.matrix(EqA))

  ## Dimensions of the problem
  Neq    <- nrow(E)    # number of equations
  Nx     <- ncol(E)    # number of unknowns
  Nineq  <- nrow(G)    # number of inequalities

  if (is.null(Nineq))
    Nineq <- 0
  if (is.null(Neq))
    Neq <- 0

  NVar   <- nrow(EqA)  # number of equations to minimise/maximise
  ## con: constraints ; rhs: right hand side
  ## First the equalities 

  con   <- E
  rhs   <- F
  dir   <- rep("==",Neq)
  if (Nineq > 0)
    {
      con   <- rbind(con,G)
      rhs   <- c(rhs,H)
      dir   <- c(dir,rep(">=",Nineq))
    }
  Range <- matrix(ncol=2,nrow=NVar,NA)

  if (ispos) {

    obj   <- vector(length = Nx)

    for (i in 1:NVar)
      {
        obj        <- EqA[i,]
        lmin       <- lp("min",obj,con,dir,rhs)
        if (lmin$status == 0) Range[i,1] <- lmin$objval else
        if (lmin$status == 3) Range[i,1] <- -1e30       else
                              Range[i,1] <- NA
        lmax       <- lp("max",obj,con,dir,rhs)
        if (lmax$status == 0) Range[i,2] <- lmax$objval else
        if (lmax$status == 3) Range[i,2] <- 1e30        else
                              Range[i,2]  <- NA
      }
  } else{
    ## First test if problem is solvable...
    Sol <- lsei(E=E,F=F,G=G,H=H)
    if (Sol$residualNorm > tol)
      { 
        warning (paste("cannot proceed: problem not solvable at requested tolerance",tol))
        return(Range)
      } 
    ## double the number of unknowns: x -> x1 -x2, x1>0 and x2>0
    con <- cbind(con,-1*con)
    EqA <- cbind(EqA,-1*EqA)

    for (i in 1:NVar) {
      obj <- EqA[i,]
      lmin <- lp("min", obj, con, dir, rhs)
      if(lmin$status == 0) Range[i, 1] <- lmin$objval else
      if(lmin$status == 3) Range[i, 1] <- -1e30 else
                           Range[i, 1] <- NA
      lmax <- lp("max", obj, con, dir, rhs)
      if(lmax$status == 0) Range[i, 2] <- lmax$objval else
      if(lmax$status == 3) Range[i, 2] <- 1e30 else
                           Range[i, 2] <- NA
    }
  }

  if (!is.null(EqB))
    {Range[,1]<-Range[,1]-EqB
     Range[,2]<-Range[,2]-EqB
   }
  colnames(Range) <- c("min","max")
  rownames(Range) <- rownames(EqA)
  return(Range)    # a 2-column matrix with the minimum and maximum value of each equation (variable)

} ########## END OF varranges ########## 

##############################################################################
##                                                                          ##
## LINEAR INVERSE MODELLING                                                 ## 
##    VARIOUS UTILITIES                                                     ##
##                                                                          ##
##############################################################################

################################################################################
## Solve       : Generalised inverse solution                                 ##
################################################################################

Solve <- function(A,                     # numeric matrix containing the coefficients of the equation Ax=B
                  B=diag(nrow=nrow(A)),  # numeric matrix containing the right-hand sides of the equation; default=unity matrix
                  tol=sqrt(.Machine$double.eps)) # tolerance for selecting singular values

  ##------------------------------------------------------------------------
  ## Generalised inverse solution of Ax=B
  ## note: solve is the R default, but requires square, positive definite A
  ##------------------------------------------------------------------------
  {
    M <-ginv(A,tol)
    if (is.null(M)) return(NULL)
    B <- matrix(nrow=nrow(A),data=B)
    X <- M %*% B
    if (ncol(B) == 1)
      {
        xnames <- colnames(A)
        X <- as.vector(X)
        names (X) <- xnames
      }
    return(X)   ## the generalised inverse solution, x, of Ax=B; if B is not specified: returns the generalised inverse of A
  }


################################################################################
## tridiag     : Solves tridiagonal system of linear equations                ##
################################################################################

Solve.tridiag  <- function(diam1,   ## (nonzero) elements below the diagonal
                           dia,     ## (nonzero) elements on the diagonal
                           diap1,   ## (nonzero) elements above the diagonal
                           rhs=rep(0,times=length(dia))) # numeric vector containing right hand side

  ##------------------------------------------------------------------------
  ## Solves a tridiagonal system of equations
  ##------------------------------------------------------------------------

  {
    Nmx <- length(dia) 
    if (length(diam1) != Nmx-1) stop("cannot solve tridiagonal problem - diam1 and dia not compatible")
    if (length(diap1) != Nmx-1) stop("cannot solve tridiagonal problem - diap1 and dia not compatible")
    if (length(rhs  ) != Nmx)   stop("cannot solve tridiagonal problem - rhs and dia not compatible")
    value=rep(0.,Nmx)
    sol <-.Fortran("tridia",Nmx=Nmx, 
                   au=as.double(c(0,diam1)),bu=as.double(dia),cu=as.double(c(diap1,0)),
                   du=as.double(rhs),value=as.double(value))

    return(sol$value)    # vector with solution of tridiagonal system

  }

################################################################################
## Solve.banded      : Solves banded system of linear equations               ##
################################################################################

Solve.banded <- function(abd,   # (nonzero) bands, rotated row-wise or full square matrix
                         nup,   # number of nonzero bands above the diagonal
                         nlow,  # number of nonzero bands below the diagonal
                         rhs=rep(0,times=ncol(abd)), # numeric vector containing right hand side
                         full=(nrow(abd)==ncol(abd))) # if true: full matrix is passed, if false:banded
  ##------------------------------------------------------------------------
  ## Solves a banded system of equations
  ##------------------------------------------------------------------------

  {
    Nmx  <- ncol(abd)
    nr   <- nrow(abd)
    if (full)     # full matrix was specified
      {
        A <- abd
        if (nrow(abd) != ncol(abd))  stop("cannot solve banded problem - nrows and ncols of abd are not the same, while the input matrix is said to be full")
        Aext <- rbind(matrix(nc=ncol(A),nr=nup,0),A,matrix(nc=ncol(A),nr=nlow,0))
        abd  <- matrix(nr=nup+nlow+1,nc=Nmx,data=
                       Aext[(col(Aext))<=row(Aext)&col(Aext)>=row(Aext)-nlow-nup])
        nr   <- nrow(abd)
      }
    Nabd <- 2*nlow+nup+1
    if (nr != nlow+nup+1)  stop("cannot solve banded problem - abd not compatible with nup and nlow")
    if (length(rhs  ) != Nmx)  stop("cannot solve banded problem - rhs and abd not compatible")

    ABD <- rbind(matrix(nr=nlow,nc=Nmx,0),abd)
    IsError <- FALSE
    sol <-.Fortran("banded",abd=ABD,beta=as.double(rhs),NumABD=as.integer(Nabd), 
                   NumSvar=Nmx,BandDown=as.integer(nlow),BandUp=as.integer(nup),
                   indx=as.integer(rep(0,Nmx)),info=as.integer(0),IsError=as.logical(IsError))
    return(sol$beta)    # vector with solution of tridiagonal system                 

  }

################################################################################
## xsample  : Samples linear problem with equality and inequality constraints ##
################################################################################

xsample <- function(A=NULL,             #Ax~=B
                    B=NULL,
                    E=NULL,              #Ex=F
                    F=NULL,    
                    G=NULL,             #Gx>H; 
                    H=NULL,             #Gx>H; 
                    sdB=1,              #standard deviations on B (weighting)
                    iter=3000,          #number of iterations
                    outputlength = iter, # number of rows of output matrix
                    burninlength = NULL,
                    type="mirror", # one of mirror, cda, da ; cda and da need to have a closed space (inequality constraints)!!
                    jmp=NULL,             #jump length of the transformed variables q: x=x0+Zq (only if type=mirror)
                    tol=sqrt(.Machine$double.eps), # accuracy of Z,g,h: smaller numbers are set to zero to avoid rounding errors
                    x0=NULL,            #particular solution
                    fulloutput=FALSE,   # provide diagnostic output such as the transformed variables q, and sample probabilities
                    test=TRUE)          # if "test" test for hidden equalities
  {

#########################################
### 1. definition of internal functions ##
#########################################
    
    mirror <- function(q1,g,h,k=length(q),jmp)
      ## function ensuring that a jump from q1 to q2
      ## fulfills all inequality constraints formulated in g and h
      ## gq=h can be seen as equations for planes that are considered mirrors. when a jump crosses one or more of these
      ## mirrors, the vector describing the jump is deviated according to rules of mirroring.
      ## the resulting new vector q will always be within the subspace of R^n for which all inequalities are met.
      ## also are the requirements for a MCMC met: the probability in the subspace is constant,
      ## the probability out of the subspace is 0.
      ## q1 has to fulfill constraints by default!
      ## Karel Van den Meersche 20070921
      {
        ##if (any((g%*%q1)<h)) stop("starting point of mirroring is not in feasible space")
        q2 <- rnorm(k,q1,jmp)
        if (!is.null(g))
          {
            residual <- g%*%q2-h
            q10 <- q1

            while (any(residual<0))                 #mirror
              {
                epsilon <- q2-q10                       #vector from q1 to q2: our considered light-ray that will be mirrored at the boundaries of the space
                w <- which(residual<0)                        #which mirrors are hit?
                alfa <- ((h-g%*%q10)/g%*%epsilon)[w]    #alfa: at which point does the light-ray hit the mirrors? g*(q1+alfa*epsilon)-h=0
                whichminalfa <- which.min(alfa)
                j <- w[whichminalfa]                    #which smallest element of alfa: which mirror is hit first?
                d <- -residual[j]/sum(g[j,]^2)     #add to q2 a vector d*Z[j,] which is oriented perpendicular to the plane Z[j,]%*%x+p; the result is in the plane.
                q2 <- q2+2*d*g[j,]                      #mirrored point
                residual <- g%*%q2-h
                q10 <- q10+alfa[whichminalfa]*epsilon   #point of reflection
              }
          }
        q2
      }

    ## hit-and-run algorithms
    ## modeled after Smith, R.L. Efficient Monte Carlo Procedures for
    ## Generating Points Uniformly Distributed over Bounded Regions. Operations Research 32, pp. 1296-1308,1984.

    cda <- function(q,g,h,k=length(q),...)            # coordinates direction algorithm
                                        # samples a new point in the feasible range of the solution space along a random coordinate
      {
        i <- sample(1:k,1)                  #
        h1 <- h-as.matrix(g[,-i])%*%q[-i]              # g[,i]q[i]>h1
        maxqi <- min((h1/g[,i])[g[,i]<0])
        minqi <- max((h1/g[,i])[g[,i]>0])
        q[i]  <-  runif(1,minqi,maxqi)
        return(q)
      }

    rda <- function(q,g,h,k=length(q),...)            #random direction algorithm
                                        #samples a new point in the feasible range of the solution space in a random direction
      {
        ##if (any((g%*%q)<h)) stop("starting point is not in feasible space")
        e <- rnorm(k)
        d <- e/norm(e)                        #d: random direction vector; q2 = q + alfa*d
        
        alfa <- ((h-g%*%q)/g%*%d)             #
        if (any(alfa>0)) alfa.u <- min(alfa[alfa>0]) else alfa.u <- 0
        if (any(alfa<0)) alfa.l <- max(alfa[alfa<0]) else alfa.l <- 0

        q.u <- q+alfa.u*d
        q.l <- q+alfa.l*d

        if (any(g%*%q.u<h)) alfa.u <- 0
        if (any(g%*%q.l<h)) alfa.l <- 0
        q <- q+runif(1,alfa.l,alfa.u)*d
        return(q)
      }

    norm <- function(x) sqrt(x%*%x)

    automatedjump <- function(a,b,g,h,g.scale=5,a.scale=5)
      {
        if (is.null(g)) s1 <- rep(NA,k)
        else
          {
            q_ranges <- xranges(E=NULL,F=NULL,g,h)
            s1 <- abs(q_ranges[,1]-q_ranges[,2])/g.scale
          }
        if (is.null(A)) s2 <- rep(NA,k)
        else
          {
            estVar <- solve(t(a)%*%diag(sdB^-2)%*%a) # estimated variance on the parameters, simplified from Brun et al 2001
            estSd  <- sqrt(diag(estVar))
            s2 <- estSd/a.scale
          }
        s <- pmin(s1,s2,na.rm=T)
        s[s>tol^-2] <- NA
        if (any (is.na(s)))
          {
           if (all(is.na(s)))
             {
               warning(" problem is unbounded - all jump lengths are set to 1")
               s[] <- 1
             } else {
               warning(" problem is unbounded - some jump lengths are set arbitrarily")
               s[is.na(s)] <- mean(s,na.rm=T)*100
             }
         }
        return(s)
      }

    
#############################
### 2. the xsample function ##
#############################

    ## conversions vectors to matrices and checks
    if (is.data.frame(A)) A <- as.matrix(A)
    if (is.data.frame(A)) A <- as.matrix(A)
    if (is.data.frame(A)) A <- as.matrix(A)
    if (is.vector(A)) A <- t(A)
    if (is.vector(E)) E <- t(E)
    if (is.vector(G)) G <- t(G)
    
    ## find a particular solution x0
    if (is.null(x0))
      {
        l <- lsei(A=A,B=B,E=E,F=F,G=G,H=H)#,Wa=1/sdB)
        if (l$residualNorm>1e-6)
          stop("no particular solution found;incompatible constraints")
        else
          x0 <- l$X
      }
    n <- length(x0)   
    
    ## KS: test for equalities, hidden in inequalities...
    if (test && !is.null(G))
      {
        xv <- varranges(E,F,G,H,EqA=G)
        ii <- which (xv[,1]-xv[,2]==0)
        if (length(ii)>0) { # if they exist: add regular equalities !
          E  <- rbind(E,G[ii,])
          F  <- c(F,xv[ii,1])

          G  <- G[-ii,]
          H  <- H[-ii]
          if (length(H)==0) G <- H <- NULL
        }
        xr <- xranges(E,F,G,H)
        ii <- which (xr[,1]-xr[,2]==0)
        if (length(ii)>0)  # if they exist: add regular equalities !
          {
            dia <- diag(nrow=nrow(xr))
            E  <- rbind(E,dia[ii,])
            F  <- c(F,xr[ii,1])
          }
      }

    ## Z is an orthogonal matrix for which E%*%Z=0; it can serve as basis for the null space of E.
    ## all solutions for the equalities have the form x = x0 + Zq with q a random vector. 
    ## the vector q is varied in a random walk, using a MCMC with acceptance rate = 1. The inequality constraints Gx>H
    ## can be rewritten as Gx0-H + (GZ)q >0
    
    if (!is.null(E))
      {
        Z <- Null(t(E)); Z[abs(Z)<tol] <- 0  #x=x0+Zq ; EZ=0
      } else { Z <- diag(n) }
    if(length(Z)==0)
      {
        warning("the problem has a single solution; this solution is returned as function value")
        return(x0)
      }
    k <- ncol(Z)

    if (!is.null(G))
      {
        g <- G%*%Z
        h <- H-G%*%x0                                            #gq-h>=0
        g[abs(g)<tol] <- 0
        h[abs(h)<tol] <- 0

      } else { g <- G; h <- H }
    

    if(!is.null(A))
      {
        a <- A%*%Z
        b <- B-A%*%x0                          #aq-b~=0
        v <- svd(a,nv=k)$v                     #transformation q <- t(v)q for better convergence
        a <- a%*%v                             #transformation a <- av
        if (!is.null(G)) g <- g%*%v            #transformation g <- gv
        Z <- Z%*%v                             #transformation Z <- Zv
        sdB <- rep(sdB,length=length(b))
        prob <- function(q) prod(dnorm(b,a%*%q,sdB))
        test <- function(q2) (prob(q2)/prob(q1))>runif(1) #metropolis criterion
      } else {
                                        # transforming g and h
                                        #        if (!is.null(G))
                                        #          {
                                        #            q_ranges <- xranges(G=g,H=h,full=T)
                                        #                                        # ##########################
                                        #                                        # # if there's not enough elements in q_ranges, we make the system bounded
                                        #                                        # ###########################
                                        #            if (ncol(q_ranges)<k+2& ncol(q_ranges)>3)
                                        #              {
                                        #                warning(" problem is unbounded - all jump lengths are set to 1")
                                        #                g2 <- rbind(g,diag(k),-diag(k))
                                        #                q_maxdist <- min(tol^-1,max(dist(q_ranges[,-(1:2)])) * 1000)
                                        #                h2 <- c(h,rep(-q_maxdist,2*k))
                                        #                q_ranges <- xranges(G=g2,H=h2,full=T)
                                        #              }
        
                                        #            q_extreme <- t(q_ranges[,-(1:2)])         # these are the values of q at the edges of the feasible space
                                        #            r <- princomp(q_extreme)$loadings   # transformation matrix r; q'=t(r)q
                                        #            g <- g%*%r
                                        #            Z <- Z%*%r

                                        #            ## v <- svd(g,nv=k)$v
                                        #            ## g <- g%*%v
                                        #            ## Z <- Z%*%v
                                        #          }
        prob <- function(q) 1
        test <- function(q2) TRUE
      }

    outputlength <- min (outputlength,iter)
    ou <- ceiling(iter/outputlength)

    q1 <- rep(0,k)
    x <- matrix(ncol=n,nrow=outputlength,dimnames=list(NULL,colnames(A)))
    x[1,] <- x0
    naccepted <- 1
    p <- vector(length=outputlength) # probability distribution
    p[1] <- prob(q1)

    if (fulloutput)
      {
        q <- matrix(ncol=k,nrow=outputlength)
        q[1,] <- q1
      }
    
    if (is.null(jmp)) jmp <- automatedjump(a,b,g,h) # automatedjump(g,h)
    if (type=="mirror") newq <- mirror
    if (type=="rda") newq <- rda
    if (type=="cda") newq <- cda

    ## ##################
    ## the random walk ##
    ## ##################
    if (!is.null(burninlength))
      for(i in 1:burninlength)
        {
          q2 <- newq(q1,g,h,k,jmp)
          if (test(q2)) q1 <- q2
        }
    for (i in 2:outputlength)
      {
        for (ii in 1:ou)
          {
            q2 <- newq(q1,g,h,k,jmp)
            if (test(q2)) { q1 <- q2
                            naccepted <- naccepted+1 }
          }
        x[i,] <- x0+Z%*%q1
        p[i] <- prob(q1)
        if (fulloutput)  q[i,] <- q1
      }
    ## ##################
    ## end random walk ##
    ## ##################

    xnames <- colnames(A)
    if (is.null(xnames)) xnames <- colnames(E)
    if (is.null(xnames)) xnames <- colnames(G)
    colnames (x) <- xnames

    xsample <- list(X=x,acceptedratio=naccepted/iter,p=p,jmp=jmp)
    if (fulloutput) xsample <- list(X=x,acceptedratio=naccepted/iter,Q=q,p=p,jmp=jmp)

    return(xsample)
  }



################################################################################
### varsample  : Samples variable equations                              ##
################################################################################

varsample <- function(X,        # matrix of valid x-values, e.g. generated with x-sample
                      EqA,      # numeric matrix containing the coefficients that define the variables
                      EqB=NULL) # numeric vector containing the right-hand side of the variable equation
{
  if (! is.matrix(EqA) & ! is.null(EqA)) EqA <- t(as.matrix(EqA))
  if (! is.matrix(X) & ! is.null(X)) X <- t(as.matrix(X))  

  if (is.null(EqB)) EqB <- rep(0,nrow(EqA))
  Var <- NULL
  if (ncol(X) != ncol(EqA)) stop("matrix X and EqA not compatible")
  for (i in 1:nrow(X))
    Var <- rbind(Var,t(EqA%*%X[i,]-EqB))
  return(Var)
}
