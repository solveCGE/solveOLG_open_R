# computes demographic transition (and updates intervivo transfers accordingly)
compdemo = function() {

  # compute demography transition
  for (tt in 2:tend) {
    Nv[1, tt]  <<- NB[tt]
    for (i in 2:nag) {
      Nv[i, tt] <<- Nv[i - 1, tt - 1] * gamv[i - 1, tt - 1]
    }
  }

  Nz   <<- per2coh(Nv)
  N    <<- aggcoh2per(Nz)
  Nc[] <<- colSums(Nv[1:(fag - 1), ])

  # Compute neutral intervivo-transfers by rescaling received transfers
  for (tt in 1:tend) {
    ivgiven          = -sum(Nv[, tt] * ivv[, tt] * (ivv[, tt] < 0))
    ivreceived       = sum(Nv[, tt] * ivv[, tt] * (ivv[, tt] > 0))

    ivv[ivv[, tt] > 0, tt] <<- ivv[ivv[, tt] > 0, tt] * (ivgiven / ivreceived)

    if (abs(sum(ivv[, tt] * Nv[, tt])) > 1e-10) cat("ERROR IN RECOMPDEMO: Unbalanced intervivo transfers!")

  }

  ivz <<- per2coh(ivv)

}

# computes the further life expectancy
lifeexpect = function(gamv) {

  nag = length(gamv)

  lifeexpectageN = zeroscol(nag)
  for (a in (nag - 1):1) {
    lifeexpectageN[a] = (lifeexpectageN[a + 1] + 1) * gamv[a] + gamv[a + 1] * (1 - gamv[a])
  }

  return(lifeexpectageN)
}

# main routine that solves the transition path of the full model
solveOLG = function(starttime = 1, maxiter = 200, tol = 1e-4, damping_budget = 1.0, damping_assets = 1.0, damping_ab = 1.0) {

  cat("\nRunning Tatonnement Algorithm for Transition:\n\n")

  tstart_loop       = proc.time()

  scaleA            = 1 # initialize
  scaleab           = onesrow(tend) # initialize

  # ===== demography ======#
  compdemo() # recomputes demographic transition

  for (iter in 1:maxiter) {
    tstart_iter      = proc.time()

    # FIRM PROBLEM (for given labor supply)
    ### terminal conditions ###
    firmSS(tend)
    firmprobR(verbose = F) # updates: K, qTob and Inv for given LD

    Y               <<- fY(K, LD, TFP)
    w               <<- MPL(K, LD, TFP) / (1 + tauF)
    wz              <<- per2coh(w, nag)
    V               <<- qTob * K
    J               <<- psi / 2 * K * (Inv / K - delta)^2
    pK              <<- MPK(K, LD, TFP)
    TaxF            <<- tauprof * (Y - (1 + tauF) * w * LD - J - delta * K)
    Div             <<- Y - (1 + tauF) * w * LD - Inv - J - TaxF

    # ===== solve the households' problem for given prices and tax rates ======#
    HHall(starttime = starttime, calibinit = (iter == 1), scaleA = scaleA)

    # ===== aggregation ======#
    Cons      <<- aggcoh2per(Consz * Nz)
    LS        <<- aggcoh2per(notretz * ellz * thetaz * Nz)
    A         <<- aggcoh2per(Az * Nz)
    ab        <<- aggcoh2per(abz * Nz)
    iv        <<- aggcoh2per(ivz * Nz) # should be 0 by construction
    Nw        <<- aggcoh2per(notretz * Nz)
    Nr        <<- aggcoh2per((1 - notretz) * Nz)

    # government budget
    P           <<- aggcoh2per((1 - notretz) * pz * Nz)
    tauW        <<- aggcoh2per(tauWz * notretz * ellz * thetaz * Nz) / LS
    TaxP        = aggcoh2per((1 - notretz) * tauWz * pz * Nz)
    Taxl        = aggcoh2per(taulz * Nz)
    Rev         <<- TaxF + (tauF * LD + tauW * LS) * w + Taxl + tauC * Cons + TaxP
    CG[]        <<- colSums(cGv * Nv)
    Exp         <<- CG + P

    # follow given debt-path
    PB[starttime:(tend - 1)]  <<- DG[starttime:(tend - 1)] - DG[(starttime + 1):tend] / (1 + r[starttime:(tend - 1)])
    PB[tend]                <<- r[tend] * DG[tend] / (1 + r[tend])

    DF[starttime:tend] <<- A[starttime:tend] - V[starttime:tend] - DG[starttime:tend]

    TB[starttime:(tend - 1)] <<- DF[(starttime + 1):tend] / (1 + r[starttime:(tend - 1)]) - DF[starttime:(tend - 1)]
    TB[tend]     <<- -DF[tend] * r[tend] / (1 + r[tend])

    # ===== excess demands ======#
    edy       <<- Inv + J + Cons + CG + TB - Y
    edg       <<- Rev - Exp - PB
    edl       <<- LD - LS
    eda       <<- DG + DF + V - A
    ediv      <<- -iv
    edab      <<- aggcoh2per((1 - gamz) * Savz * Nz) - ab
    edw       <<- 1 * edy + w * edl + ediv + edab + edg + eda - c(eda[2:tend], eda[tend]) / (1 + r) # Walras' Law

    # check Walras' Law: this always has to hold (even out of equilibrium)! If not there is something wrong with accounting in the model
    # if (max(abs(edw[starttime:(tend-1)]))> 1e-10) stop("Error: Walras Law does not hold!")

    tstop_iter      = proc.time()

    # ===== checking error and breaking loop ======#
    err             = sum(abs(edy[starttime:tend])) + sum(abs(edg[starttime:tend])) + sum(abs(edl[starttime:tend])) + sum(abs(eda[starttime:tend])) + sum(abs(ediv[starttime:tend])) + sum(abs(edab[starttime:tend]))
    err2            = log(err / tol)

    cat(paste0("Iteration:  ", formatC(iter, width = 3), "   scaleA: ", format_dec(scaleA, 6), "   scaleab: ", format_dec(mean(scaleab), 6), "   non-conv.HH: ", format_dec(sum(HH_nonconvz), 0), "   Time: ",  format_dec(tstop_iter[3] - tstart_iter[3], 3), " sec   log of err/tol: ", format_dec(err2, 8), "\n"))

    if (err2 < 0.0) {
      cat(paste0(rep(" ", 107), collapse = ""), "Convergence!\n\n")
      break
    }
    if (iter == maxiter) {
      cat(paste0(rep(" ", 107), collapse = ""), "No Convergence!\n\n")
      break
    }

    HH_nonconvz[] <<- 0 # reset convergence counter

    # ======= updating for next iteration =======#
    # budget rules
    budget_surplus  = edg * damping_budget

    if (budget_bal == 1) {
      tauWv     <<- tauWv - kronecker(budget_surplus / (w * LS), onescol(nag))
      tauWz     <<- per2coh(tauWv)
    }
    if (budget_bal == 2) {
      tauF       <<- tauF - budget_surplus / (w * LD)
    }
    if (budget_bal == 3) {
      tauC       <<- tauC - budget_surplus / Cons
      tauCv      <<- kronecker(tauC, onescol(nag))
      tauCz      <<- per2coh(tauCv)
    }
    if (budget_bal == 4) {
      taul            <<- taul - budget_surplus / (N - Nc)
      taulv[fag:nag, ] <<- kronecker(taul, onescol(nag - fag + 1))
      taulz           <<- per2coh(taulv)
    }
    if (budget_bal == 5) {
      tauprof    <<- tauprof - budget_surplus / (Y - (1 + tauF) * w * LD - J - delta * K)
    }
    if (budget_bal == 6) {
      cGv        <<- cGv + kronecker(budget_surplus / N, onescol(nag))
      CG         <<- colSums(cGv * Nv)
    }

    scaleab         = 1 + (aggcoh2per((1 - gamz) * Savz * Nz) / ab - 1) * damping_ab
    abv             <<- abv * kronecker(matrix(scaleab, 1, tend), onescol(nag))
    abz             <<- per2coh(abv)
    LD              <<- LS
    scaleA          = 1 + ((DG[starttime] + DF[starttime] + V[starttime]) / A[starttime] - 1) * damping_assets

  }

  # convert cohort-view variables back to period-view variables
  # (those where only cohort-view variables were altered in solveOLG)
  Av       <<- coh2per(Az)
  Consv    <<- coh2per(Consz)
  lambdav  <<- coh2per(lambdaz)
  Savv     <<- coh2per(Savz)
  dis_totv <<- coh2per(dis_totz)
  ellv     <<- coh2per(ellz)
  rv       <<- coh2per(rz)
  wv       <<- coh2per(wz)
  pcv      <<- coh2per(pcz)
  yv       <<- coh2per(yz)

  tstop_loop       = proc.time()

  cat(paste0("Computation time:\t", format_dec(tstop_loop[3] - tstart_loop[3], 4), " sec\n"))
  cat(paste0("CHECK SOLUTION:\t\t", max(abs(edy) + abs(edl) + abs(edg) + abs(eda) + abs(ediv) + abs(edab))), "\n")

}
