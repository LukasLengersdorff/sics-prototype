### Constraints ###

#'
#'
#' @export
compute_U = function(M, orthonormalize = TRUE) {
  ## Add matrix U to matrix M such that the extended matrix rbind(M,U) is invertible
  E = tryCatch({matlib::echelon(diag(ncol(M)) - t(M) %*% solve(M %*% t(M)) %*% M)},
               error = function(e) {
                 if (det(M %*% t(M)) < 1e-10) stop("The specified constraints are not linearly independent. M = [A; B; C] is rank-deficient")
               })
  allzeros = apply(E, 1, function(x) all(x==0))
  U = E[!allzeros, ,drop=FALSE]
  if (orthonormalize) {
    if (nrow(U) == 1) {U = U/sqrt(sum(U^2))}
    if (nrow(U) > 1) {U = t(matlib::GramSchmidt(t(U)))}
  }
  U
}

#' Define a set of linearly independent constraints
#'
#' Define a set of equality, inequality, and/or interval constraints for a parameter of dimension \eqn{n}, using Matrix equation syntax.
#'
#' The `Constraint` class is used to define sets of linearly independent equality, inequality, and interval
#' constraints. `Constraint` objects are created by this function, and by all flavors of helper functions
#' (e.g., `EqualityConstraint`, `LinearConstraint`, etc.). Constraints can be combined/intersected with the
#' `+` operator, and this is the suggested way to build up constraints.
#'
#' `Constraint()` acts as the constructor for the `Constraint` class. Directly using it
#'  is the most flexible, but also the least "user-friendly" way to specify constraints
#'  (e.g. the exact matrices and vectors need to be figured out beforehand).
#'
#' @param A,a A \eqn{k \times n} matrix `A` and a \eqn{n} vector `a`, defining the equality
#' constraint \eqn{A x = a}
#'
#' @param B,b A \eqn{k \times n} matrix `B` and a \eqn{n} vector `b`, defining the inequality
#' constraint \eqn{B x > b}
#'
#' @param C,c0,c1 A \eqn{k \times n} matrix `C` and \eqn{n} vectors `c0` and `c1`, defining the interval
#' constraint \eqn{C x \in (c_0,c_1)}
#'
#' @param n An integer greater than 0. The dimension of the parameter. Only needed when `A`,`B`, and `C` are undefined (thus, only for empty constraints)
#'
#' @param a_default,b_default Real numbers. If `a` or `b` are left undefined, they default to vectors with `n` repetitions of `a_default`/`b_default`.
#'
#' @returns An object of class `Constraint`.
#'
#' @examples
#' # NULL
#'
#' @export
#'
Constraint = function(A=NULL, B=NULL, C=NULL,
                      a=NULL,
                      b=NULL,
                      c0=NULL, c1=NULL,
                      n = NULL,
                      a_default = 0,
                      b_default = 0) {

  ## Promote vectors to matrix rows
  if (is.vector(A)) {A = matrix(A,nrow=1)}
  if (is.vector(B)) {B = matrix(B,nrow=1)}
  if (is.vector(C)) {C = matrix(C,nrow=1)}

  ## Derive dimension n from matrices, if needed
  if (is.null(n)) {
    if (!is.null(A)) {n = ncol(A)} else if (!is.null(B)) {n = ncol(B)} else if (!is.null(C)) {n = ncol(C)} else {stop("If A, B, and C are not specified, n needs to be specified.")}
  }
  ## Deal with null matrices, all-zero matrices, or empty matrices A, B, C
  if (is.null(A) || all(A==0) || length(A)==0) {A = matrix(0, nrow=0, ncol = n)}
  if (is.null(B) || all(B==0) || length(B)==0) {B = matrix(0, nrow=0, ncol = n)}
  if (is.null(C) || all(C==0) || length(C)==0) {C = matrix(0, nrow=0, ncol = n)}

  #-- Check inputs ----
  if (nrow(A)==0) {
    mA = 0
    a = array(dim = 0)
  } else {
    mA = nrow(A)
    if (is.null(a)) a = rep(a_default, times = mA)
  }

  if (nrow(B)==0) {
    mB = 0
    b = array(dim = 0)
  } else {
    mB = nrow(B)
    if (is.null(b)) b = rep(b_default, times = mB)
  }

  if (nrow(C)==0) {
    mC = 0
    c0 = array(dim = 0)
    c1 = array(dim = 0)
  } else {
    mC = nrow(C)
    if (is.null(c0) || is.null(c1)) {
      stop("If C is specified, both c0 and c1 need to be specified as well")
    } else {
      if (length(c0) != mC) stop("c0 and c1 need to have length nrow(C)")
      if (length(c0) != length(c1)) stop("c0 and c1 are not of same length")
      if (any(c1 <= c0)) stop("c1 must be strictly greater than c0")
    }
  }

  mU = n - (mA + mB + mC)

  ## Calculate stuff

  # Here would be space to check for a "prepare" tag and skip calculations & set M, MA... = NULL if prepare=FALSE.
  # Could be interesting to allow representation of linearly dependent constraints for additional utility.

  if (mU == n) { # Special handling of "no constraints at all"
    U = diag(n)
    empty = TRUE
  } else {
    U = compute_U(rbind(A, B, C))
    empty = FALSE
  }
  M = solve(rbind(A,B,C,U))

  colnames(M) = NULL

  if (mA == 0) {iA = 0} else {iA = 1:mA}
  if (mB == 0) {iB = 0} else {iB = (mA+1):(mA+mB)}
  if (mC == 0) {iC = 0} else {iC = (mA+mB+1):(mA+mB+mC)}
  if (mU == 0) {iU = 0} else {iU = (mA+mB+mC+1):n}

  structure(list(M = M,
                 MA = M[,iA,drop=FALSE], MB = M[,iB,drop=FALSE], MC = M[,iC,drop=FALSE], MU = M[,iU,drop=FALSE],
                 A = A, B = B, C = C, U = U,
                 a = a, b = b, c0 = c0, c1 = c1,
                 n = n, mA = mA, mB = mB, mC = mC, mU = mU),
            class = "Constraint", empty = empty)
}

#'
#'
#' @export
compute_logdet_jacobian = function(con) {
  MO = cbind(con$MB, con$MC, con$MU)
  0.5*log(det(crossprod(MO)))
}


#' Define an empty constraint
#'
#' Define an empty constraint of dimension \eqn{n}, i.e., a constraint that lets the n-dimensional parameter vary freely.
#'
#' This is useful for sampling from unconstrained encompassing models, and as a basis for building up more complex constraints
#'
#' @param n An integer. The dimension of the parameter
#' @returns An object of class `Constraint`
#' @examples
#' # NULL
#'
#' @export
EmptyConstraint = function(n) {
  Constraint(n = n)
}

#' Row-bind padded with zeros
#'
#' A version of `rbind` that pads the shorter of the two vectors or matrices with zeroes.
#'
#' @param x,y Two vectors or matrices
#' @returns A matrix
#' @examples
#' # NULL
#'
#' @export
rbind_padded = function(x,y) {
  nx = ncol(x); ny = ncol(y)
  if (any(nx == ny, is.null(nx), is.null(ny))) return(rbind(x,y))
  n = max(nx,ny)
  x2 = matrix(0,nrow=nrow(x),ncol=n)
  x2[,1:nx] = x
  y2 = matrix(0,nrow=nrow(y),ncol=n)
  y2[,1:ny] = y
  #message("Matrix was padded with zeros.\n")
  return(rbind(x2,y2))
}

#'
#'
#' @export
`+.Constraint` = function(con1,con2) {
  Constraint(A = rbind_padded(con1$A, con2$A),
             B = rbind_padded(con1$B, con2$B),
             C = rbind_padded(con1$C, con2$C),
             a = c(con1$a, con2$a),
             b = c(con1$b, con2$b),
             c0 = c(con1$c0, con2$c0),
             c1 = c(con1$c1, con2$c1))
}

#'
#'
#' @export
print.Constraint = function(con) {
  cat        ("-----------------------------------------------\n")
  cat(sprintf("--- Constraint on %i-dimensional parameter x ---\n", con$n))
  cat        ("-----------------------------------------------\n\n")
  if (attr(con, "empty")) {
    cat("Empty constraint. All components of x may vary \nfreely.\n\n")
    cat("-----------------------------------------------")
    return(invisible(NULL))
  }
  if (con$mA > 0) {
    cat(sprintf("%i equality constraint%s:\n\n", con$mA, ifelse(con$mA==1,"","s")))
    #cat("--- Equality constraints ---\n\n")
    .printhelper(con$A, con$a, "=")
    cat("\n")
  }
  if (con$mB > 0) {
    cat(sprintf("%i inequality constraint%s:\n\n", con$mB, ifelse(con$mB==1,"","s")))
    #cat("--- Inequality constraints ---\n\n")
    .printhelper(con$B, con$b, ">")
    cat("\n")
  }
  if (con$mC > 0) {
    cat(sprintf("%i interval constraint%s:\n\n", con$mC, ifelse(con$mC==1,"","s")))
    #cat("--- Interval constraints ---\n\n")
    .printhelper2(con$C, con$c0, con$c1)
    cat("\n")
  }
  cat("-----------------------------------------------")
}

### print helpers ###
.printhelper = function(m,v,s,zerostring=".", rounding_digits = 3) {
  m = round(m, rounding_digits)
  v = round(v, rounding_digits)
  maxl_m = max(nchar(m))
  form_m = sprintf("%%%is",maxl_m)
  maxl_v = max(nchar(v))
  form_v = sprintf("|%%%is|",maxl_v)
  nr = nrow(m)
  nc = ncol(m)
  mid = floor(nr/2+1)
  for (i in 1:nr) {
    cat(" |")
    for (j in 1:nc) {
      cat(sprintf(form_m, ifelse(m[i,j]==0, zerostring, m[i,j])))
      if (j < nc) cat(" ")
    }
    if (i == mid) {
      cat(sprintf("| * x %s ",s))
    } else {
      cat("|       ")
    }
    cat(sprintf(form_v, v[i]))
    if (i == mid) {
      #cat(".\n")
      cat(" \n")
    } else {
      cat(" \n")
    }
  }
}

#'
.printhelper2 = function(m,v,v2, zerostring=".", rounding_digits = 3) {
  m = round(m, rounding_digits)
  v = round(v, rounding_digits)
  v2 = round(v2, rounding_digits)
  maxl_m = max(nchar(m))
  form_m = sprintf("%%%is",maxl_m)
  maxl_v = max(nchar(v))
  maxl_v2 = max(nchar(v2))
  form_v = sprintf("(%%%is, %%%is)",maxl_v, maxl_v2)
  nr = nrow(m)
  nc = ncol(m)
  mid = floor(nr/2+1)
  for (i in 1:nr) {
    cat(" |")
    for (j in 1:nc) {
      cat(sprintf(form_m, ifelse(m[i,j]==0, zerostring, m[i,j])))
      if (j < nc) cat(" ")
    }
    if (i == mid) {
      cat(sprintf("| * x %s ","in"))
    } else {
      cat("|        ")
    }
    cat(sprintf(form_v, v[i], v2[i]))
    if (i == mid) {
      #cat(".\n")
      cat(" \n")
    } else {
      cat(" \n")
    }
  }
}




##### Utility constraint builders #####

equalityMatrix = function(indices, n = max(indices)){
  M = matrix(0, ncol = n, nrow = length(indices)-1)
  for (i in 1:(length(indices)-1)) {
    M[i, indices[i]] = 1
    M[i, indices[i+1]] = -1
  }
  M
}

indexMatrix = function(indices, n = max(indices)){
  M = matrix(0, ncol = n, nrow = length(indices))
  for (i in 1:length(indices)) M[i,indices[i]] = 1
  M
}

#' ...
#'
#' @export
EqualityConstraint = function(indices, n = max(indices)){
  if (length(indices) <= 1) stop("indices must be a vector of length > 1")
  # Could also default to an empty constraint for length 1 indices, but I think
  # this is the safer approach
  A = equalityMatrix(indices, n)
  Constraint(A = A, a = numeric(length(indices)-1))
}

#'
#'
#' @export
OrdinalConstraint = function(ordered_indices, n = max(ordered_indices)){
  if (length(ordered_indices) <= 1) stop("ordered_indices must be a vector of length > 1")
  # Could also default to an empty constraint for length 1 indices, but I think
  # this is the safer approach
  B = equalityMatrix(ordered_indices, n)
  Constraint(B = B, b = numeric(length(ordered_indices)-1))
}

#'
#'
#' @export
ConstantConstraint = function(indices, constant = 0, n = max(indices)){
  A = indexMatrix(indices,n)
  Constraint(A = A, a = rep(constant, length(indices)))
}

#'
#'
#' @export
MeanConstraint = function(n, m = 0){
  A = matrix(1/n, nrow=1, ncol=n)
  Constraint(A = A, a = m)
}



.expression2vector = function(string, vars, append_const = FALSE) {
  n = length(vars)
  v = numeric(n) %>% as.list
  names(v) = vars
  vec = numeric(n)
  const = with(v, eval(str2expression(string)))
  for (i in seq_len(n)) {
    v2 = v; v2[[vars[i]]] = 1
    vec[i] = with(v2, eval(str2expression(string)))
  }
  vec = vec - const
  if (append_const) vec[n+1] = const
  return(vec)
}

#'
#'
#'... Cannot handle interval constraints
#' @export
LinearConstraint = function(equation_string, vars = NULL, n = NULL, default_var = "x") {
  # Default handling
  if (is.null(vars)) {
    # If vars == NULL, but n supplied, set vars = x0, x1, ..., xn. Else, derive n from string
    if (is.null(n)) {
      pattern = paste0(default_var, "[0-9]+")
      n = stringr::str_extract_all(equation_string, pattern, simplify = TRUE) %>%
        stringr::str_extract_all("[0-9]+", simplify = TRUE) %>%
        as.numeric %>%
        max
    }
    vars = paste0(default_var, seq_len(n))
  }

  # Simplify and flatten equations
  eqs = unlist(strsplit(equation_string, "[,;]"))
  num_eq_signs = sum(sapply(eqs, function(i)
    length(unlist(
      gregexpr("[=<>]", i)
    ))))
  if (num_eq_signs > length(eqs)) {
    .eqs = character(num_eq_signs)
    k = 1
    for (cur_eq in eqs) {
      eq_sign_idx = unlist(gregexpr("[=<>]", cur_eq))
      eq_signs = sapply(eq_sign_idx, function(i)
        substr(cur_eq, i, i))
      exprs = strsplit(cur_eq, "[=<>]")[[1]]
      tidy_exprs = character(length(eq_signs))
      for (i in 1:length(eq_signs)) {
        .eqs[k] = paste(exprs[i], eq_signs[i], exprs[i + 1])
        k = k + 1
      }
    }
    eqs = .eqs
  }
    eqs = stringr::str_replace_all(eqs,"[[:whitespace:]]", "")

    # Extract left hand side, sign, and right hand side
    .charmatrix = stringr::str_match(eqs, "(.+)([=<>])(.+)")[,-1,drop=FALSE]
    .lhs = .charmatrix[,1]
    .sgn = .charmatrix[,2]
    .rhs = .charmatrix[,3]
    # Flip "<" to ">"
    lhs = if_else(.sgn == "<", .rhs, .lhs)
    rhs = if_else(.sgn == "<", .lhs, .rhs)
    sgn = if_else(.sgn == "<", ">", .sgn)
    # All to lefthand side
    expressions = sprintf("%s - (%s)", lhs, rhs)
    M = sapply(expressions, .expression2vector, vars = vars, append_const = TRUE)
    M = t(M)
    const = M[,ncol(M)] * (-1)
    M = M[,-ncol(M),drop=FALSE]

    A = M[sgn == "=",]
    B = M[sgn == ">",]
    a = const[sgn == "="]
    b = const[sgn == ">"]

    Constraint(A=A,B=B,a=a,b=b)
  }


#'
#'
#' @export
IntervalConstraint = function(equation_string, c0 = 0, c1 = 1,
                              vars = NULL, n = NULL, default_var = "x") {
  if (length(equation_string)>1) stop("IntervalConstraint can only handle one constraint per call.")
  # Default handling
  if (is.null(vars)) {
    # If vars == NULL, but n supplied, set vars = x0, x1, ..., xn. Else, derive n from string
    if (is.null(n)) {
      pattern = paste0(default_var, "[0-9]+")
      n = stringr::str_extract_all(equation_string, pattern, simplify = TRUE) %>%
        stringr::str_extract_all("[0-9]+", simplify = TRUE) %>%
        as.numeric %>%
        max
    }
    vars = paste0(default_var, seq_len(n))
  }

  .C = .expression2vector(string = equation_string, vars = vars, append_const = TRUE)
  C = .C[-length(.C)]
  const = .C[length(.C)]

  Constraint(C=C, c0=c0-const, c1=c1-const)

}


#'
#'
#' @export
transform_par = function(con, theta, unbounded = FALSE, g = list(B = log, C = qlogis)) {
  suppressWarnings({
  wb = NULL
  wc = NULL
  wu = NULL
  if (con$mB > 0) wb = con$B %*% theta - con$b
  if (con$mC > 0) wc = ((con$C %*% theta) - con$c0)/(con$c1 - con$c0)
  if (con$mU > 0) wu = con$U %*% theta
  })

  if (con$mA > 0) {if(any(c(con$A %*% theta) != con$a)) warning("Equality constraint not satisfied!")}
  if (con$mB > 0) {if(any(wb <= 0)) warning("Inequality constraint not satisfied! ")}
  if (con$mC > 0) {if(any(wc < 0 | wc > 1)) warning("Interval constraint not satisfied! ")}

  if (unbounded) { suppressWarnings({
    if (con$mB > 0) wb = g$B(wb)
    if (con$mC > 0) wc = g$C(wc)
  })}

  c(wb, wc, wu)
}

#'
#'
#' @export
transform_par_inv = function(con, w, unbounded = FALSE, g_inv = list(B = exp, C = plogis)) {
  wb = w[seq_len(con$mB)]
  wc = w[(con$mB) + seq_len(con$mC)]
  wu = w[(con$mB+con$mC) + seq_len(con$mU)]

  if (unbounded) {
    if (length(wb) > 0) wb = g_inv$B(wb)
    if (length(wc) > 0) wc = g_inv$C(wc)
  }

  y = c(con$a,
        wb + con$b,
        (con$c1-con$c0)*wc + con$c0,
        wu)
  c(con$M %*% y)
}


#'
#'
#' @export
check_constraint = function(con, theta, detailed_output = FALSE) {
    wb = NULL
    wc = NULL
    wa = NULL
    if (con$mA > 0) wa = c(c(con$A %*% theta) == con$a)
    if (con$mB > 0) wb = c(con$B %*% theta > con$b)
    if (con$mC > 0) wc = c(con$c0 < con$C %*% theta & con$C %*% theta < con$c1)

  if(!detailed_output) return(all(c(wa,wb,wc)))
  list(A = wa, B = wb, C = wc)
}
