# ============================================================
# Piecewise Exponential Approximation
# True hazard: h(t) = t^2
# ============================================================

# ------------------------------------------------------------
# 1. Generate right-censored survival data
# ------------------------------------------------------------

set.seed(2026)

n <- 400
censor.rate <- 0.3

# True cumulative hazard:
# H(t) = t^3 / 3
#
# Since H(T) ~ Exp(1),
# T = (3E)^(1/3), E ~ Exp(1)

E <- rexp(n)
T <- (3 * E)^(1/3)

# Independent censoring time
C <- rexp(n, rate = censor.rate)

# Observed survival data
U <- pmin(T, C)
delta <- as.integer(T <= C)

dat <- data.frame(
  U = U,
  delta = delta
)

head(dat)

# Proportion censored
mean(delta == 0)


# ------------------------------------------------------------
# 2. Function for fitting a piecewise exponential model
# ------------------------------------------------------------

fit.pwe <- function(U, delta, cuts, grid) {

  # Partition:
  # [0, cut1), [cut1, cut2), ..., [last cut, Inf)

  v <- c(0, cuts, Inf)

  J <- length(v) - 1

  r <- numeric(J)
  W <- numeric(J)

  for (j in 1:J) {

    # Number of observed events in interval j
    r[j] <- sum(
      delta == 1 &
      U >= v[j] &
      U < v[j + 1]
    )

    # Total observed time at risk in interval j
    W[j] <- sum(
      pmax(
        pmin(U, v[j + 1]) - v[j],
        0
      )
    )
  }

  # Piecewise hazard estimates
  lambda.hat <- ifelse(
    W > 0,
    r / W,
    NA_real_
  )

  # Determine the interval corresponding to each grid point
  index <- findInterval(
    grid,
    c(0, cuts)
  )

  # Estimated hazard function over grid
  h.hat <- lambda.hat[index]

  list(
    r = r,
    W = W,
    lambda.hat = lambda.hat,
    h.hat = h.hat
  )
}


# ------------------------------------------------------------
# 3. Define grid and user-specified cut points
# ------------------------------------------------------------

grid <- seq(
  0,
  2,
  length.out = 200
)

# Fewer intervals:
# 4 intervals = 4 hazard parameters
cuts.few <- c(
  0.75,
  1.25,
  1.75
)


# Fewer intervals:
# 4 intervals = 4 hazard parameters
cuts.few <- c(
 0.25,
0.5, 
0.75,
1,
  1.25,
1.5,
  1.75,
2
)



# More intervals:
# 7 intervals = 7 hazard parameters
cuts.many <- c(
  0.50,
  0.75,
  1.00,
  1.25,
  1.50,
  1.75
)


# ------------------------------------------------------------
# 4. Fit both piecewise exponential models
# ------------------------------------------------------------

fit.few <- fit.pwe(
  U = U,
  delta = delta,
  cuts = cuts.few,
  grid = grid
)

fit.many <- fit.pwe(
  U = U,
  delta = delta,
  cuts = cuts.many,
  grid = grid
)


# ------------------------------------------------------------
# 5. Examine estimated parameters
# ------------------------------------------------------------

# Number of events in each interval
fit.few$r
fit.many$r

# Total time at risk in each interval
fit.few$W
fit.many$W

# Estimated piecewise hazards
fit.few$lambda.hat
fit.many$lambda.hat

# Approximate variance of log(lambda_hat_j)
1 / fit.few$r
1 / fit.many$r


# ------------------------------------------------------------
# 6. Plot:
#    true hazard versus piecewise exponential estimates
# ------------------------------------------------------------

ylim <- range(
  c(
    grid^2,
    fit.few$h.hat,
    fit.many$h.hat
  ),
  na.rm = TRUE
)

plot(
  grid,
  grid^2,
  type = "l",
  lwd = 2,
  ylim = ylim,
  xlab = "t",
  ylab = "Hazard"
)

lines(
  grid,
  fit.few$h.hat,
  type = "s",
  lty = 2,
  lwd = 2
)

lines(
  grid,
  fit.many$h.hat,
  type = "s",
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True hazard",
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(1, 2, 3),
  lwd = 2
)


# ============================================================
# 7. Repeated simulation
# ============================================================

set.seed(2026)

B <- 1000

est.few <- matrix(
  NA_real_,
  nrow = length(grid),
  ncol = B
)

est.many <- matrix(
  NA_real_,
  nrow = length(grid),
  ncol = B
)

for (b in 1:B) {

  # Generate survival times from h(t) = t^2
  E.b <- rexp(n)

  T.b <- (3 * E.b)^(1/3)

  # Generate independent censoring times
  C.b <- rexp(
    n,
    rate = censor.rate
  )

  # Observed data
  U.b <- pmin(
    T.b,
    C.b
  )

  delta.b <- as.integer(
    T.b <= C.b
  )

  # Fit model with fewer intervals
  fit.few.b <- fit.pwe(
    U = U.b,
    delta = delta.b,
    cuts = cuts.few,
    grid = grid
  )

  # Fit model with more intervals
  fit.many.b <- fit.pwe(
    U = U.b,
    delta = delta.b,
    cuts = cuts.many,
    grid = grid
  )

  est.few[, b] <- fit.few.b$h.hat

  est.many[, b] <- fit.many.b$h.hat
}


# ------------------------------------------------------------
# 8. Monte Carlo mean estimated hazard
# ------------------------------------------------------------

mean.few <- rowMeans(
  est.few,
  na.rm = TRUE
)

mean.many <- rowMeans(
  est.many,
  na.rm = TRUE
)

plot(
  grid,
  grid^2,
  type = "l",
  lwd = 2,
  xlab = "t",
  ylab = "Hazard"
)

lines(
  grid,
  mean.few,
  type = "s",
  lty = 2,
  lwd = 2
)

lines(
  grid,
  mean.many,
  type = "s",
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "True hazard",
    "Mean: fewer intervals",
    "Mean: more intervals"
  ),
  lty = c(1, 2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 9. Empirical standard deviation of estimated hazard
# ------------------------------------------------------------

sd.few <- apply(
  est.few,
  1,
  sd,
  na.rm = TRUE
)

sd.many <- apply(
  est.many,
  1,
  sd,
  na.rm = TRUE
)

plot(
  grid,
  sd.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(sd.few, sd.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical SD of estimated hazard"
)

lines(
  grid,
  sd.many,
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 10. Optional: empirical bias
# ------------------------------------------------------------

true.hazard <- grid^2

bias.few <- mean.few - true.hazard

bias.many <- mean.many - true.hazard

plot(
  grid,
  bias.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(bias.few, bias.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical Bias"
)

lines(
  grid,
  bias.many,
  lty = 3,
  lwd = 2
)

abline(
  h = 0,
  lty = 1
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)


# ------------------------------------------------------------
# 11. Optional: empirical MSE
# ------------------------------------------------------------

mse.few <- rowMeans(
  (est.few - true.hazard)^2,
  na.rm = TRUE
)

mse.many <- rowMeans(
  (est.many - true.hazard)^2,
  na.rm = TRUE
)

plot(
  grid,
  mse.few,
  type = "l",
  lty = 2,
  lwd = 2,
  ylim = range(
    c(mse.few, mse.many),
    na.rm = TRUE
  ),
  xlab = "t",
  ylab = "Empirical MSE"
)

lines(
  grid,
  mse.many,
  lty = 3,
  lwd = 2
)

legend(
  "topleft",
  legend = c(
    "Fewer intervals",
    "More intervals"
  ),
  lty = c(2, 3),
  lwd = 2
)