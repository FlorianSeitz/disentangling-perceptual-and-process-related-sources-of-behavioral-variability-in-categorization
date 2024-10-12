# ==========================================================================
# Simulation Effect of Noise on Categorization
# ==========================================================================
if (!require("pacman")) install.packages("pacman")
pacman::p_load(data.table, purrr, Rsolnp, doRNG, MASS, stringr, patchwork, ggplot2, ggtext, truncnorm)
set.seed(123)

# ==========================================================================
structure <- "rb" # ii or rb
nf <- 2 # number of features
s <- 1/10 # modification factor

low <- list(ca = c(92.5, 127.5), # centroid of category A (medium-low)
            cb = c(127.5, 92.5)) # centroid of category B (medium-low)

medium <- list(ca = c(101.5, 118.5), # centroid of category A (medium)
               cb = c(118.5, 101.5)) # centroid of category B (medium)

high <- list(ca = c(106.1, 113.9), # centroid of category A (medium-high)
             cb = c(113.9, 106.1)) # centroid of category B (medium-high)

low <- lapply(low, "*", s); medium <- lapply(medium, "*", s); high <- lapply(high, "*", s)
means <- mget(c("low", "medium", "high"))

# ca_mean <- c(92.5, 127.5) # centroid of category A (medium-low)
# cb_mean <- c(127.5, 92.5) # centroid of category B (medium-low)
# 
# ca_mean <- c(101.5, 118.5) # centroid of category A (medium)
# cb_mean <- c(118.5, 101.5) # centroid of category B (medium)
# 
# ca_mean <- c(106.1, 113.9) # centroid of category A (medium-high)
# cb_mean <- c(113.9, 106.1) # centroid of category B (medium-high)
# 
# ca_mean <- ca_mean * s; cb_mean <- cb_mean * s

var <- 162.5; covar <- 112 # variances and covariances of each category
var <- var * s^2; covar <- covar * s^2
c_cov <- matrix(c(var, covar, covar, var), nrow = nf)

n <- 100 # number of exemplars per category
nreps <- 20 # number of repetitions per transfer stimulus

# Parameters (range): attention w1 (0-1), sensitivity l (0-2), response scaling tau (0-5)
max_l <- 2; max_tau <- 2
n_w1 <- sqrt(1/12) # attention noise (1 SD within range 0-1)
n_l <- sqrt(max_l^2/12) # sensitivity noise (1 SD within range 0-2)
n_tau <- sqrt(max_tau^2/12) # response-scaling noise (1 SD within range 0-5)
n_x <- sqrt(c(.51, 2.75)) # sqrt(var) # perceptual intake noise (1 category SD)
n_s <- sqrt(1/12)

# n_w1 <- n_w1/2; n_l <- n_l/2; n_tau <- n_tau/2; n_x <- n_x/2 # makes noise small
# ==========================================================================

# ==========================================================================
# Specifies goodness-of-fit function
# ==========================================================================
fits <- function(pars, train, trans, value = c("pred", "resp", "ll", "ll_cont"), resps = NULL, noise = "none") {
  r <- 2
  if (value == "pred") {
    trans <- trans[rep(1:nrow(trans), each = nreps)]
  }
  preds <- apply(trans, 1, function(i) {
    w1 <- pars["w1"]
    l <- pars["l"]
    tau <- pars["tau"]
    if(noise == "x") i <- i + runif(n = 2, min = -n_x, max = n_x)
    if(noise == "w") w1 <- min(c(max(c(w1 + runif(1, -n_w1, n_w1), 0)), 1))
    if(noise == "l") l <- min(c(max(c(l + runif(1, -n_l, n_l), 0)), max_l))
    if(noise == "tau") tau <- min(c(max(c(tau + runif(1, -n_tau, n_tau), 0)), max_tau))
    sims <- train[, .(s = exp(-l * (w1 * abs(f1 - i[["f1"]])^r + (1-w1) * abs(f2 - i[["f2"]])^r)^(1/r))), by = c]
    if(noise == "s") sims[, s := pmin(pmax(s + runif(2*n, -n_s, n_s), 0), 1)]
    sims[, .(s = sum(s)), by = c][, 1 / (1 + exp(tau * (s[c == 0] - s[c == 1])))]
  })
  
  if(value == "pred") return(
    rtruncnorm(n = length(preds), a = 0, b = 1, mean = preds, sd = pars["sigma"])
    # return(preds)
  ) 
  if(value == "resp") return(rbinom(length(preds), 1, preds))
  if(value == "ll") return(-sum(dbinom(x = resps, prob = rep(preds, each = nreps), size = 1, log = TRUE)))
  if(value == "ll_cont") return(
    -sum(log(dtruncnorm(x = resps, a = 0, b = 1, mean = rep(preds, each = nreps), sd = pars["sigma"])))
  )
}

# ==========================================================================
# Makes stimuli
# ==========================================================================
samples <- function(ca_mean, cb_mean, structure = "ii", n = 100) {
  # Makes learning stimuli
  ca <- cbind(mvrnorm(n = n, mu = ca_mean, Sigma = c_cov), 1)
  cb <- cbind(mvrnorm(n = n, mu = cb_mean, Sigma = c_cov), 0)
  train <- as.data.table(rbind(ca, cb))
  colnames(train) <- c(paste0("f", 1:nf), "c")
  
  # Makes transfer grid (+/- 4 SDs from overall centroid in 2 SD steps)
  grid_mean <- colMeans(rbind(ca_mean, cb_mean))
  range <- lapply(grid_mean, function(i) {i + c(-1, 1) * 4 * sqrt(var)})
  values <- sapply(1:length(range), function(i) {seq(min(range[[i]]), max(range[[i]]), by = 2 * sqrt(var))})
  trans <- as.data.table(expand.grid(f1 = values[, 1], f2 = values[, 2]))
  
  # rotating the learning stimuli 45 degrees counterclock wise
  if (structure == "rb") {
    theta <- 45 * pi/180
    t_matrix <- matrix(c(cos(theta), -sin(theta), sin(theta), cos(theta)), nrow = 2)
    train[, c("f1", "f2") := as.data.table(
      as.matrix(train[, .(f1, f2)] - grid_mean) %*% t_matrix + grid_mean)][]
  }
  
  return(list(train, trans))
}


# ==========================================================================
# Runs simulation
# ==========================================================================
noises <- "x" # noises <- c("none", "x", "w", "l", "s")
n_iter <- 500

sim_grid <- expand.grid(id = 1:n_iter, overlap = c("low", "medium", "high"))

start_time <- proc.time()
res <- as.data.table(do.call("rbind", lapply(1:1500, function(j) {
  print(paste0("-------------- Iter: ", j, "/", nrow(sim_grid), " --------------"))
  pars <- c(w1 = runif(1), l = runif(1, max = max_l), 
            tau = runif(1, max = max_tau),
            sigma = runif(1, max = .5))
  print(pars)
  curr_means <- means[[sim_grid[j, "overlap"]]]
  stimuli <- samples(ca_mean = curr_means[["ca"]], cb_mean = curr_means[["cb"]], structure = structure)
  # stimuli[[1]] <- stimuli[[1]][, .(f1 = mean(f1), f2 = mean(f2)), by = c] # makes a prototype model
  do.call("rbind", lapply(noises, function(noise) {
    preds <- fits(pars = pars, 
                  train = stimuli[[1]], 
                  trans = stimuli[[2]], 
                  value = "pred", 
                  noise = noise)
    # resps <- rbinom(length(preds), 1, preds)
    # m <- solnp(pars = c(w1 = .5, l = 1, tau = 1), fun = fits, LB = c(0, .001, .001), UB = c(1, max_l, max_tau),
    #            train = stimuli[[1]], trans = stimuli[[2]], value = "ll", resps = resps, noise = "none")
    m <- solnp(pars = c(w1 = .5, l = 1,
                        tau = 1,
                        sigma = .25), fun = fits,
               LB = c(0, .001, .001, .001), UB = c(1, max_l, max_tau, .50),
               train = stimuli[[1]], trans = stimuli[[2]], value = "ll_cont", resps = preds, noise = "none")
    data.table(task = structure,
               overlap = sim_grid[j, "overlap"],
               noise = noise,
               iter = j,
               true.w1 = pars[["w1"]],
               true.l = pars[["l"]],
               true.tau = ifelse("tau" %in% names(pars), pars[["tau"]], as.numeric(NA)),
               true.sigma = ifelse("sigma" %in% names(pars), pars[["sigma"]], as.numeric(NA)),
               fit.w1 = m$pars[["w1"]],
               fit.l = m$pars[["l"]],
               fit.tau = ifelse("tau" %in% names(m$pars), m$pars[["tau"]], as.numeric(NA)),
               fit.sigma = ifelse("sigma" %in% names(m$pars), m$pars[["sigma"]], as.numeric(NA)),
               gof = tail(m$values, 1),
               id = rep(1:nrow(stimuli[[2]]), each = nreps),
               # resps = resps,
               true.preds = preds,
               fit.preds = fits(pars = m$pars, train = stimuli[[1]], trans = stimuli[[2]], value = "pred")
    )
  }))
})))
timetaken(start_time)

fwrite(res, "~/Projects/P1/data/final/rb-continuous-x.csv")
