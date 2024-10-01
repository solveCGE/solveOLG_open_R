# production function and marginal products
fY     = function(K, L, TFP) TFP * K^alpha * L^(1 - alpha)
MPK    = function(K, L, TFP) TFP * alpha * (K / L)^(alpha - 1)
MPKK   = function(K, L, TFP) TFP * alpha * (alpha - 1) * K^(alpha - 2) * L^(1 - alpha)
MPL    = function(K, L, TFP) TFP * (1 - alpha) * (K / L)^alpha

firmSS = function(tt) {
  pK[tt]     <<- (r[tt] + delta * (1 - tauprof[tt])) / (1 - tauprof[tt])
  K[tt]      <<- (pK[tt] / alpha / TFP[tt])^(1 / (alpha - 1)) * LD[tt]
  Inv[tt]    <<- delta * K[tt]
  qTob[tt]   <<- (1 + r[tt])
}

firmTE = function(vars, tt) {

  # contemporaneous
  K1    = vars[1]
  qTob1 = vars[2]
  Inv1  = vars[3]

  # leads
  K2    = vars[4]
  qTob2 = vars[5]
  # Inv2  = vars[6]

  resid    = zeroscol(nendo)
  partials = zerosmat(nendo, nvar)

  # residuals
  resid[1] = K2 - (1 - delta) * K1 - Inv1
  resid[2] = (1 + r[tt]) * (1 + (1 - tauprof[tt]) * psi * (Inv1 / K1 - delta)) - qTob2
  resid[3] = (1 - tauprof[tt]) * MPK(K1, LD[tt], TFP[tt]) - (1 - tauprof[tt]) * psi / 2 * (delta^2 - (Inv1 / K1)^2) + tauprof[tt] * delta + qTob2 / (1 + r[tt]) * (1 - delta) - qTob1

  # partials
  neq = 1

  # eq 1:
  partials[neq, 1] = -(1 - delta)
  # partials[neq,2] = 0
  partials[neq, 3] = -1
  partials[neq, 4] = 1
  # partials[neq,5] = 0
  neq = neq + 1

  # eq 2:
  partials[neq, 1] = -(1 + r[tt]) * (1 - tauprof[tt]) * psi * Inv1 / (K1)^2
  # partials[neq,2] = 0
  partials[neq, 3] = (1 + r[tt]) * (1 - tauprof[tt]) * psi / K1
  # partials[neq,4] = 0
  partials[neq, 5] = -1
  neq = neq + 1

  # eq 3:
  partials[neq, 1] = (1 - tauprof[tt]) * MPKK(K1, LD[tt], TFP[tt]) - (1 - tauprof[tt]) * (psi * (Inv1 / K1)^2) / K1
  partials[neq, 2] = -1
  partials[neq, 3] = (1 - tauprof[tt]) * (psi * (Inv1 / K1)^2) / Inv1
  # partials[neq,4] = 0
  partials[neq, 5] = 1 / (1 + r[tt]) * (1 - delta)
  neq = neq + 1

  return(cbind(partials, resid))

}

firmprobR = function(maxiter = 20L, tol = 1e-10, verbose = F) {

  yy  = rbind(K, qTob, Inv)  ## use initial steady state as guess
  yyy = matrix(yy, ncol = 1); yyy = matrix(yyy[(nstate + 1):(ntot - ncont), 1], ncol = 1)

  # create Jacobian and Resids
  Jac    = Matrix(0.0, nrow = ntot2, ncol = ntot2, sparse = T)
  Res    = Matrix(0.0, nrow = ntot2, ncol = 1, sparse = F)

  err    = 10.0
  iter   = 1L

  while (err > tol) {

    # Period 1
    tt = 1

    ins     = c(K[1], yyy[1:(nvar - nstate)])
    outTemp = firmTE(ins, tt)

    Res[1:nendo]                  = outTemp[, nvar + 1]
    Jac[1:nendo, 1:(nvar - nstate)]  = outTemp[, (nstate + 1):nvar]

    # Periods 2:(tend-2)
    for (tt in 2:(tend - 2)) {

      ins     = yyy[(nendo * (tt - 1) - nstate + 1):(nendo * (tt - 1) - nstate + nvar)]
      outTemp = firmTE(ins, tt)

      Res[(nendo * (tt - 1) + 1):(nendo * tt)]                                                     = outTemp[, nvar + 1]
      Jac[(nendo * (tt - 1) + 1):(nendo * tt), (nendo * (tt - 1) - nstate + 1):(nendo * (tt - 1) - nstate + nvar)]  = outTemp[, 1:nvar]

    }

    # Period tend-1
    tt = tend - 1

    ins     = c(yyy[(nendo * (tt - 1) - nstate + 1):ntot2], qTob[tend]) # OK

    outTemp = firmTE(ins, tt)

    Res[(nendo * (tt - 1) + 1):(nendo * tt)]                                                              = outTemp[, nvar + 1]
    Jac[(nendo * (tt - 1) + 1):(nendo * tt), (nendo * (tt - 1) - nstate + 1):(nendo * (tt - 1) - nstate + nendo * 2 - ncont)]  = outTemp[, 1:(nendo * 2 - ncont)]

    # Newton update of guess
    yyy = yyy - solve(Jac, Res)

    err = max(abs(Res))
    if (verbose == T) cat(paste0("Iteration: ", iter, " \t\tError: ", err, "\n"))

    if (err <= tol) {
      if (verbose == T) cat("Convergence!\t\t")
      break
    }

    iter = iter + 1

    if (iter > maxiter) {

      cat("No convergence! Maximum iterations reached!\t\t")
      break
    }

  }

  yyy = rbind(K[1], yyy, qTob[tend], Inv[tend])
  yy  = matrix(yyy, nrow = nendo, ncol = tend)

  K    <<- yy[1, ]
  qTob <<- yy[2, ]
  Inv  <<- yy[3, ]

}

### testing ###
# ins0    = c(K0,qTob0,Inv0,K0,qTob0)
# testout = firmTE(ins0,1)
###############
