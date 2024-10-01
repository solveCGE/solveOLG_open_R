cat("Start calibration: \n")

# parameters
delta                    = 0.05                 # depreciation rate
psi                      = 10.0                 # scaling capital adjustment costs
r0                       = 0.04                 # real interest rate
rho                      = 0.035                # time discount rate
sigma                    = 0.9                  # elasticity of inter-temporal substitution
sigL                     = 0.3                  # labor supply elasticity

# note: ages are off-set by 1 year, e.g. age group 1 contains 0-year olds
fag                      = 14                   # first economically active age-group (age 15)
rag0                     = 61.3                 # retirement age group (retirement age 62.3), non-whole numbers allowed
iag0                     = 51                   # first age group giving inter-vivo transfers

ivpc                     = 0.2                  # intervivo transfer received per capita

# some normalizations
N0                       = 100.0                # population
Y0                       = 100.0                # GDP
L0                       = 30.0                 # total labor supply in efficiency units
w0                       = 2.0                  # wage rate
Cons0                    = 0.55 * Y0              # consumption share (calibrated using taul0)

### DEMOGRAPHY ###
for (i in 1:nag) {
  gamv0[i] = 1 - 0.89^(nag - i + 1)                   # some simple profile
}

# survival of last age group is 0
gamv0[nag] = 0

# compute demography
Nv0[1]  = 1
for (i in 2:nag) {
  Nv0[i] = Nv0[i - 1] * gamv0[i - 1]
}

# rescale population
NB0             = 1 / sum(Nv0) * N0
Nv0             = Nv0 / sum(Nv0) * N0

avage0 = sum(Nv0 * (0:(nag - 1))) / N0
report("REPORT: Average age:", avage0)
lifeexpect0 = lifeexpect(gamv0)
report("REPORT: Life-expectancy at age 0:", lifeexpect0[1])
report("REPORT: Life-expectancy at age 65:", lifeexpect0[66])

### AGE PROFILES ###

# indicator for not-retired
notretv0[1:floor(rag0)]   = 1                       # not retired
notretv0[floor(rag0) + 1]   = rag0 - floor(rag0)        # partly retired

# intervivo-transfers
ivv0[iag0:nag]            = -seq(from = ivpc, to = ivpc * 2, length.out = nag - iag0 + 1) # some increasing profile (from ivpc to 2*ivpc)
ivv0[fag:(iag0 - 1)]        = -sum(ivv0[iag0:nag] * Nv0[iag0:nag]) / sum(Nv0[fag:(iag0 - 1)]) * onescol(iag0 - fag)

iv0                       = sum(ivv0 * Nv0)
if (abs(iv0) > 1e-10) stop("ERROR: UNBALANCED INTERVIVO TRANSFERS!")

thetav0                     = zeroscol(nag)                               # labor productivity parameters
theta_peak                  = floor(rag0) - 10                              # assumption: productivity peaks 10 years before retirement
thetav0[fag:theta_peak]     = seq(from = 0.7, to = 1, length.out = (theta_peak - fag + 1))
thetav0[(theta_peak + 1):nag] = seq(from = 1, to = 0.1, length.out = (nag - theta_peak))

ellv0                       = L0 / sum(Nv0 * thetav0 * notretv0) * onescol(nag)   # labor supply

# partition of population
Nc0   = sum(Nv0[1:(fag - 1)])        # number of children
Nw0     = sum(notretv0 * Nv0) - Nc0                   # number of workers
Nr0     = sum((1 - notretv0) * Nv0)             # number of retirees
report("REPORT: Old-age dependency ratio:", sum(Nv0[66:nag]) / sum(Nv0[16:65]))
report("REPORT: Economic dependency ratio:", (Nc0 + Nr0) / Nw0)
report("CHECK: Newborns - deaths:", sum((1 - gamv0) * Nv0) - NB0)
report("CHECK: Children + workers + retriees - pop.:", Nc0 + Nw0 + Nr0 - N0)

### POLICY PARAMETERS ###
tauWv0                   = 0.15 * onescol(nag)            # wage tax rate worker & retiree
tauF0                    = 0.2                          # payroll tax rate
tauC0                    = 0.2                          # consumption tax rate
tauprof0                 = 0.1                          # profit tax rate
pv0                      = 0.65 * sum(w0 * ellv0 * thetav0 * Nv0) / N0 * onescol(nag)  # old-age pension (65% of average wage earnings)
DG0                      = 0.6 * Y0                       # government debt level (60% of GDP)

# cGv0 is used to balance budget in calibration
cGv0_profile             = 0.2 * onescol(nag)
cGv0_profile[1:25]       = seq(from = 0.4, to = 0.2, length.out = 25)
cGv0_profile[55:nag]     = seq(from = 0.2, to = 1.0, length.out = nag - 55 + 1) # some U-shaped profile

# price of consumption and age specific prices and tax rates (but the same for all age groups)
pc0     = 1 + tauC0
tauCv0  = tauC0 * onescol(nag)
pcv0    = pc0 * onescol(nag)
wv0     = w0 * onescol(nag)
rv0     = r0 * onescol(nag)

LS0     = sum(notretv0 * ellv0 * thetav0 * Nv0) # aggregate labor supply
LD0     = LS0
pK0     = (r0 + delta * (1 - tauprof0)) / (1 - tauprof0) # user-cost of capital
K0      = (Y0 - (1 + tauF0) * w0 * LD0) / pK0
Inv0    = delta * K0
J0      = 1 / 2 * psi * K0 * (Inv0 / K0 - delta)^2
alpha   = K0 * pK0 / (K0 * pK0 + LS0 * ((1 + tauF0) * w0))
qTob0   = (1 - tauprof0) * alpha * Y0 / K0 + tauprof0 * delta + (1 - delta) # = 1+r0 # using J0 = 0
TFP0    = Y0 / ((K0^alpha) * (LS0^(1 - alpha)))
# LD0     = ((1-alpha)*TFP0/((1+tauF0)*w0))^(1/alpha)*K0 # also true
TaxF0   = tauprof0 * (Y0 - (1 + tauF0) * w0 * LD0 - J0 - delta * K0)
Div0    = Y0 - (1 + tauF0) * w0 * LD0 - Inv0 - J0 - TaxF0
V0      = (1 + r0) * Div0 / r0

calibfind <- function(xcalib0) {

  retvar = zeroscol(5)

  # rho     <<- xcalib0[1]
  DF0     <<- xcalib0[1]
  cGscale   = xcalib0[2]
  taul0   <<- xcalib0[3]
  ab0     <<- xcalib0[4]
  lambdain  = xcalib0[5]

  abv0[fag:nag]    <<- ab0 / (N0 - Nc0) * onescol(nag - fag + 1) # children do not receive accidental bequest (workers start out with 0 assets)
  taulv0[fag:nag]  <<- taul0
  cGv0             <<- cGv0_profile + cGscale

  # INCOME
  yv0     <<- notretv0 * (wv0 * (1 - tauWv0) * ellv0 * thetav0) + (1 - notretv0) * (1 - tauWv0) * pv0 - taulv0

  # CONSUMPTION FOR ALL AGE GROUPS

  # Euler equation
  lambdav0[fag] <<- lambdain
  for (a in fag:(nag - 1)) {
    lambdav0[a + 1] <<- lambdav0[a] / ((1 / (1 + rho)) * gamv0[a] * (1 + rv0[a]))
  }
  Consv0[fag:nag] <<- (pcv0[fag:nag] * lambdav0[fag:nag])^(-sigma)

  # assets
  Av0[fag] <<- 0
  for (a in (fag + 1):nag) {
     Av0[a]     <<- (1 + rv0[a - 1]) * (Av0[a - 1] + yv0[a - 1] + ivv0[a - 1] + abv0[a - 1] - pcv0[a - 1] * Consv0[a - 1])
  }
  Savv0   <<- Av0 + yv0 + ivv0 + abv0 - pcv0 * Consv0

  # AGGREGATION
  A0      <<- sum(Av0 * Nv0)                                # total assets
  P0      <<- sum((1 - notretv0) * pv0 * Nv0)                   # expend pensions
  CG0     <<- sum(cGv0 * Nv0)                               # government consumption
  Exp0    <<- CG0 + P0                                      # total primary expenditures
  tauW0   <<- sum(tauWv0 * notretv0 * ellv0 * thetav0 * Nv0) / LS0  # average wage tax rate
  Rev0    <<- TaxF0 + (tauF0 * LD0 + tauW0 * LS0) * w0 + taul0 * (Nw0 + Nr0) + tauC0 * Cons0 + sum((1 - notretv0) * tauWv0 * pv0 * Nv0) # total revenues
  PB0     <<- DG0 * r0 / (1 + r0)                               # primary balance
  TB0     <<- -DF0 * r0 / (1 + r0)                              # trade balance

  # EXCESS DEMANDS
  edy0    <<- Cons0 + CG0 + Inv0 + J0 + TB0 - Y0  # goods market
  edl0    <<- LD0 - LS0                   # labor market
  eda0    <<- DG0 + DF0 + V0 - A0             # assets market
  edg0    <<- Rev0 - Exp0 - PB0             # government budget
  ediv0   <<- -iv0                      # intervivo transfers resource constraint
  edab0   <<- sum((1 - gamv0) * Savv0 * Nv0) - ab0 # accidental bequest resource constraint

  retvar[1]   = edy0
  retvar[2]   = edg0
  retvar[3]   = sum(Consv0 * Nv0) - Cons0
  retvar[4]   = edab0
  retvar[5]   = Savv0[nag]

  return(retvar)

}

# MATCH CALIBRATION TARGETS
# xcalib0 = c(0.01, 0.3719, 0.40, 13, 1) # starting guesses for multiroot()
xcalib0 = c(1.0, 0.3719, 0.40, 13, 1) # starting guesses for multiroot()

xout = multiroot(calibfind, xcalib0, rtol = 1e-8)
if (abs(xout$estim.precis) > 1e-6) stop("NEWTON METHOD DID NOT CONVERGE!\n")

### CALIBRATION OF LABOR SUPPLY MARGINS ###

# set parl0 in order to reproduce ell0, FOC ell0
parlv0  = (wv0 * (1 - tauWv0) * thetav0 / pcv0) * (ellv0^(-1 / sigL)) * (Consv0^(-1 / sigma)); parlv0[1:(fag - 1)] = 0
# set parl1 in order to normalize disutility of labor to 0
parlv1  = (sigL / (1 + sigL)) * parlv0 * (ellv0^((1 + sigL) / sigL))
dis_totv0 = (sigL / (1 + sigL)) * parlv0 * (ellv0^((1 + sigL) / sigL)) - parlv1

report("REPORT: Asset-to-output ratio:", A0 / Y0)

checkA0         = sum(Av0 * Nv0) - A0
checkAv0        = Av0[nag] + yv0[nag] + ivv0[nag] + abv0[nag] - pc0 * Consv0[nag] # end of period assets of last age group are zero
checkN0         = sum(Nv0) - N0

chkcalib = c(edy0, edl0, edg0, ediv0, eda0, edab0, checkA0, checkAv0, checkN0)

report("CHECK: Calibration:", sum(chkcalib))

# fill time-dependent variables with calibration values
varsfill    = c("Cons", "DG", "DF", "Inv", "LD", "LS", "K", "N", "NB", "PB", "TFP", "TB", "ab", "pc", "pK", "qTob", "r", "rag", "tauC", "tauF", "tauW", "taul", "tauprof", "w")
fillvars(agedep = F, varsfill)

varsfill_a  = c("A", "Cons", "N", "Sav", "ab", "cG", "ell", "gam", "iv", "lambda", "notret", "p", "tauC", "tauW", "taul", "theta", "r", "w")
fillvars(agedep = T, varsfill_a)

## some optional plots of the calibration

plot(1:nag - 1, Av0, type = "l", xlab = "age", ylab = "assets")

plot(1:nag - 1, Consv0, type = "l", xlab = "age", ylim = c(0, 1), ylab = "")
lines(1:nag - 1, notretv0 * ellv0 * thetav0 * wv0 * (1 - tauWv0), col = "blue")
lines(1:nag - 1, (1 - notretv0) * pv0 * (1 - tauWv0), col = "red")
lines(1:nag - 1, cGv0, col = "green")
legend("topright", legend = c("consumption", "net labor income", "net pension income", "public consumption"), col = c("black", "blue", "red", "green"), lty = 1)
