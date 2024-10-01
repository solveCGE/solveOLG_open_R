# HH.R.. solve household problem for given prices wz[,z], abz[,z], taxes, etc.

HH_root = function(lambdain, sage = fag, z = 1) {

  # EULER EQUATION: solve forward in age
  lambdaz[sage, z]    <<- lambdain
  if (sage < nag) {
    for (a in sage:(nag - 1)) {
      lambdaz[a + 1, z] <<- lambdaz[a, z] / ((1 / (1 + rho)) * gamz[a, z] * (1 + rz[a, z]))
    }
  }

  # CONSUMPTION
  pcz[sage:nag, z]      <<- 1 + tauCz[sage:nag, z]
  Consz[sage:nag, z]    <<- (pcz[sage:nag, z] * lambdaz[sage:nag, z])^(-sigma)

  # HOURS SUPPLY
  ellz[sage:nag, z]     <<- ((wz[sage:nag, z] * (1 - tauWz[sage:nag, z]) * thetaz[sage:nag, z] / pcz[sage:nag, z] * (Consz[sage:nag, z]^(-1 / sigma))) / parlv0[sage:nag])^sigL
  dis_totz[sage:nag, z] <<- (sigL / (1 + sigL)) * parlv0[sage:nag] * ellz[sage:nag, z]^((1 + sigL) / sigL) - parlv1[sage:nag]

  # CONSUMPTION AND SAVINGS
  yz[sage:nag, z]       <<- notretz[sage:nag, z] * (wz[sage:nag, z] * (1 - tauWz[sage:nag, z]) * ellz[sage:nag, z] * thetaz[sage:nag, z]) + (1 - notretz[sage:nag, z]) * (1 - tauWz[sage:nag, z]) * pz[sage:nag, z] - taulz[sage:nag, z]

  # ASSETS: solve forward in age
  Az[1, z]         <<- 0

  if (sage < nag) {
    for (i in sage:(nag - 1)) {
      Az[i + 1, z]   <<- (1 + rz[i, z]) * (Az[i, z] + yz[i, z] + ivz[i, z] + abz[i, z] - pcz[i, z] * Consz[i, z]) # if sage > 1 take previous age entry in Az as starting value! (i.e. has to be given globally not passed in function)
    }
  }
  Savz[sage:nag, z]  <<- Az[sage:nag, z] + yz[sage:nag, z] + ivz[sage:nag, z] + abz[sage:nag, z] - pcz[sage:nag, z] * Consz[sage:nag, z]

  return(Savz[nag, z])
}

HH = function(sage = fag, z = 1, maxiter = 30, stol = 1e-10, atol = 0.1) {

  err            = Inf
  iter           = 0
  trys           = 0
  stepsize       = 1e-6 # for numerical gradient

  lambdatrys     = c(1.0, 0.5, 1.5, 0.25, 1.25, 0.1, 1.0)
  maxtrys        = length(lambdatrys)
  while_continue = TRUE

  while (while_continue) {

    while_continue = FALSE
    lambdazsave    = lambdaz[sage, z]

    while (((err > stol) || (abs(Savz[nag, z]) > atol)) && (trys < maxtrys)) {

      trys = trys + 1
      iterpertry = 0
      lambdaz1 = lambdazsave * lambdatrys[trys]

      breakwhile = FALSE
      while ((err > stol) && (iterpertry < maxiter) && (breakwhile == FALSE)) {
        if (iterpertry == 0) { # Newton step for first iteration
          f2 = HH_root(lambdaz1 + stepsize, sage, z); iter = iter + 1

          if (!is.finite(f2)) {breakwhile = TRUE; break}
          f1 = HH_root(lambdaz1, sage, z); iter = iter + 1

          if (!is.finite(f1)) {breakwhile = TRUE; break}
          lambdaz2 = lambdaz1 - f1 * stepsize / (f2 - f1)
          if (!is.finite(lambdaz2) || (lambdaz2 < 0)) {breakwhile = TRUE; break}
        } else { # Secant method
          f1 = HH_root(lambdaz1, sage, z); iter = iter + 1

          if (!is.finite(f1)) {breakwhile = TRUE; break}
          lambdaz2 = lambdaz1 - f1 * (lambdaz1 - lambdaz0) / (f1 - f0)
          if (!is.finite(lambdaz2) || (lambdaz2 < 0)) {breakwhile = TRUE; break}
        }
        err = abs(lambdaz2 - lambdaz1)
        lambdaz0 = lambdaz1
        lambdaz1 = lambdaz2
        f0       = f1
        iterpertry = iterpertry + 1
      }
    }
  }

  if (abs(Savz[nag, z]) > atol) {
    HH_nonconvz[nag, z] <<- 1 # counter
  }

}

HHall = function(starttime = 1, calibinit = F, scaleA = 1) {
  for (z in starttime:ncoh) {
    if (z <= nag - fag + starttime - 1) {
      if (calibinit == T) {
        Az[, z] <<- Av0
      }
      Az[nag - (z - starttime), z] <<- Az[nag - (z - starttime), z] * scaleA
      HH(nag - (z - starttime), z)
    } else {
      HH(z = z)
    }
  }
}
