pacman::p_load(TMB, mvabund, MASS, splm, matrixcalc, tidyverse, RTMB, vegan, Matrix, doParallel, gllvm, ggh4x, gravity)
rm(list=ls())

# Assumptions ------------------
# Data generating process: ZAGA with layer-dependent covariates x_ij^(k), layer-independent covariates w_ij and latent variables u_k, with log link on mu_ij^(k) and dispersion parameter sigma. Zero-inflation parameter pi_ij.

"
GGLLVM / GLAMLE (Jiang, La Vecchia, Rastelli):
   Y_ij^(k) | u_k, x_ij^(k), w_ij  ~ ZAGA(mu_ij^(k), sig, pi_ij)
   log mu_ij^(k) = eta1_ij^k = alpha_ij^T (1, u_k)^T + beta_ij^T x_ij^(k)
   log(sig) = eta2; sig
   logit(pi_ij) = eta3_ij
   Gamma (mu,sig) like in Rigby et al.:
     shape = 1/sig^2
     scale = mu*sig^2
   u_k ~ N(0, I_q)
2-time optimization; first pi then marginal gamma without pi
Laplace approximation integrates out u = (u_1,...,u_K) via random = u
Start from true values + random noise
"

# Load all the necessary functions
files <- list.files(path = "R/Simulations/Functions",
                    pattern = "\\.R$",
                    full.names = TRUE)

# Quietly source the files
invisible(lapply(files, function(file) {
  suppressMessages(suppressWarnings(source(file)))
}))

# Set parameters for the simulation (toy example). Parameters used in the paper are p=10, q=1, K=500, L=1, Lprime=1, MC.size=100, family.ind=4 (ZAGA), seed = 2022.

params1 = list(
  p = 3,              # number of nodes
  q = 1,              # number of latent variables
  K = 0.5e2,          # number of layers
  L = 1,              # number of observable covariates (k dependent)
  Lprime = 1,         # number of observable covariates (k INdependent)
  MC.size = 2,        # number of simulations
  family.ind = 4,     # family: 1=Bernoulli, 2=Poisson, 3=ZIP, 4= ZAGA
  seed = 2022         # seed for reproducibility      
)


# Simulation function:

mc_zaga_xij_bij <- function(par, seed) {
  true_par = gen_params(par)
  getAll(par, true_par, warn = T)
  
  # Starting values (true+noise)
  min1 = -0.05
  set.seed(seed)
  noise = runif(length(a_to_start(alpha)), min = min1, max = -min1)
  noise3 <- matrix(runif(length(z.2), min = min1, max = -min1), nrow = K, ncol = q)
  parameters1 <- list(
    lambda = a_to_start(alpha) + noise,
    u = z.2 + noise3)
  
  if (family.ind == 4){
    noise7 = runif(1, min = min1, max = -min1)
    parameters1$logsig = log(sig + noise7)
  }
  
  if (L != 0) {
    noise4 = matrix(runif(m*L, min = min1, max = -min1), nrow = m, ncol = L)
    parameters1$beta <- beta + noise4
  } else {
    parameters1$beta <- matrix(0, nrow = 0, ncol = m)
  }
  
  if (Lprime > 0) {
    noise2 = runif(r, min = min1, max = -min1)
    noise5 = runif(Lprime, min = min1, max = -min1)
    parameters1$intercept = zeta + noise2
    parameters1$gamma = gamma + noise5
  } else {
    noise2 <- runif(length(alpha[1, 2:m]), min = min1, max = -min1)
    parameters1$intercept = alpha[1, 2:m] + noise2
    parameters1$gamma <- matrix(0, nrow = 1, ncol = Lprime)
  }
  
  # Setup for parallel processing
  no_cores <- detectCores(logical = TRUE)
  cl <- makeCluster(min(no_cores, 4))  # Use 4 cores or fewer if not available
  registerDoParallel(cl)
  ptm <- proc.time()
  
  # Print the number of workers to confirm setup
  cat("Number of cores used:", getDoParWorkers(), "\n")
  
  results_df = data.frame(matrix(NA, MC.size, 2))
  colnames(results_df) <- c("value", "convergence")
  
  # Storage lists for results
  
  results_beta <- vector("list", MC.size)
  results_gamma <- vector("list", MC.size)
  results_z <- vector("list", MC.size)
  results_y <- vector("list", MC.size)
  results_alpha <- vector("list", MC.size)
  results_mu <- vector("list", MC.size)
  results_mu.hat <- vector("list", MC.size)
  results_ey.hat <- vector("list", MC.size)
  results_z_true <- vector("list", MC.size)
  results_zeta <- vector("list", MC.size)
  results_sig <- numeric(MC.size)
  results_pi <- vector("list", MC.size)
  results_pi_SE <- vector("list", MC.size)
  results_gradients <- vector("list", MC.size)
  results_mu_ppml <- vector("list", MC.size)
  coef_ppml <- vector("list", MC.size)
  results_ppml_df <- vector("list", MC.size)
  
  # family.ind = 2 # Family index: 1 for Bernoulli, 2 for Poisson, 3 for ZIP, 4 for ZAGA
  
  # Set up the progress bas
  pb = txtProgressBar(min = 0, max = MC.size, style = 3)
  
  # Run MC simulations
  for (sim in 1:MC.size) {
    
    tryCatch({
      # create Y, X and Z
      graph <- graph_gen_zaga_x_ij_mixK(a = alpha, b = beta, g = gamma, n = K, X = X, z = z, Xprime = (Xprime_mL), seed = seed + sim, pi.par = pi.par, sig = sig)
      y.sim <- graph$y
      results_y[[sim]] <- y.sim
      mu.sim <- graph$mu1
      sig.sim <- graph$sig
      pi.sim <- graph$pi.par
      x.sim <- graph$X
      z2.sim <- graph$lvs
      z.sim <- cbind(rep(1,K),z2.sim)
      layers = paste0("k", 1:K)
      nodes = paste0("e", 1:p)
      EY.sim <- sweep(mu.sim, 2, 1 - pi.sim, FUN = "*")
      df.edge = expand.grid(i = nodes, j = nodes) %>% filter(i!=j) %>% arrange((i)) %>% mutate(edges = paste0(i, "-", j))
      edges = df.edge %>% pull(edges)
      rownames(y.sim) = layers
      colnames(y.sim) = edges
      rownames(EY.sim) = layers; colnames(EY.sim) = edges
      W = rep(1,K) %*% t(Xprime_mL)
      df = data.frame(y = c(y.sim), x1 = c(x.sim[,,1]), w1 = exp(c(W)), layer = rep(paste0("k", 1:K), times = m), edge = rep(edges, each = K))
      
      df = df %>% mutate(i = rep(sapply(1:m*K, function(x) str_split(df$edge, "-")[[x]][1]), each = K), j = rep(sapply(1:m*K, function(x) str_split(df$edge, "-")[[x]][2]), each = K), pi = rep(pi.sim, each = K))
      
      stopifnot(all.equal(df$y[(K+1):(2*K)], as.numeric(y.sim[,2])))
      stopifnot(all.equal(df %>% filter(layer == "k3", edge == "e2-e3") %>% pull(y), y.sim["k3", "e2-e3"]))
      
      #### 1. PI #####
      
      pi.hat <-  colMeans(y.sim == 0)     # length m
      names(pi.hat) <- edges
      results_pi[[sim]] <- pi.hat
      results_pi_SE[[sim]] <- sqrt((pi.hat * (1-pi.hat))/K)
      
      #### 2. GAMMA #####
      # Set up the data
      data_input <- list(y = y.sim, q = q)
      
      if (L != 0) {
        data_input$x <- X
      } else {
        data_input$x <- matrix(0, nrow = K, ncol = 0)
      }
      
      if (Lprime != 0) {
        data_input$xprime <- Xprime_mL
        data_input$H <- H
      } else {
        data_input$xprime <- matrix(0, nrow = m, ncol = 0); data_input$H <- matrix(0, nrow = m, ncol = r)
      }
      
      # Optimize objective function
      objr <- RTMB::MakeADFun(func = nll(glamle_nll_zaga_x_ij_mixK_pos, data_input),
                              parameters = parameters1, 
                              random = "u",
                              silent = T 
      )
      
      out <- nlminb(
        start     = objr$par,
        objective = objr$fn,
        gradient  = objr$gr,
        control   = list(
          eval.max = 5000,
          iter.max = 2000,
          rel.tol  = 1e-10,
          x.tol    = 1e-10
        )
      )
      
      results_df[sim,1] = out$objective
      results_df[sim,2] = out$convergence
      
      # Extract the gradients at the optimized parameters
      grad.sim <- objr$gr(out$par)
      
      # Extract parameter names
      param_names <- names(out$par)
      
      # Extract and store the theta.hat
      par.hat = objr$env$parList(out$par)
      par.hat$alpha.hat = .get_alpha_hat(objr,out,H,q,m)
      par.hat$sig = exp(par.hat$logsig)
      
      mu.hat = .compute_mu(par.hat = par.hat, X = X, Xprime_mL = Xprime_mL, clamp = F)
      ey.hat = sweep(mu.hat, 2, 1 - pi.hat, FUN = "*")
      results_ey.hat[[sim]] <- ey.hat
      # Save results for this simulation
      results_beta[[sim]] <- par.hat$beta
      results_gamma[[sim]] <- par.hat$gamma
      results_z[[sim]] <- par.hat$u
      results_z_true[[sim]] <- z.sim
      results_alpha[[sim]] <- par.hat$alpha.hat
      results_zeta[[sim]] <- par.hat$intercept
      results_mu[[sim]] <- mu.sim
      results_mu.hat[[sim]] <- mu.hat
      results_sig[sim] <- par.hat$sig
      results_gradients[[sim]] <- grad.sim
      
      #### 3. PPMLE #####
      fit_ppml_gravity <- ppml(
        dependent_variable = "y",
        distance = "w1",
        additional_regressors = c("x1", "i", "j"),
        #family = quasipoisson(link = "log"),
        data = df
      )
      
      df = df %>% mutate(ey = c(EY.sim), ey.hat = c(ey.hat), pi.hat = rep(c(pi.hat), each = K), true_mu = c(mu.sim), mu.hat = c(mu.hat))
      stopifnot(all.equal(df %>% filter(layer == "k3", edge == "e2-e3") %>% pull(ey), EY.sim["k3", "e2-e3"]))
      stopifnot(all.equal(pi.hat["e2-e3"], (df %>% filter(layer == "k3", edge == "e2-e3") %>% pull(pi.hat))[1]))
      
      # Save results for this simulation
      df = df %>% mutate(mu_ppml = predict(fit_ppml_gravity, type = "response"))
      results_ppml_df[[sim]] <- df
      coef_ppml[[sim]] = coef(fit_ppml_gravity)
      results_mu_ppml[[sim]] = predict(fit_ppml_gravity, type = "response")
      
      
    }, error = function(e) {
      cat("Error in simulation", sim, ": ", e$message, "\n")
    })
    setTxtProgressBar(pb, value = sim)
  }
  
  
  
  # Stop cluster after use
  stopCluster(cl)
  close(pb)
  
  # Run time
  ptm.end = proc.time()
  cat("Simulation completed. Time taken:\n")
  time.taken = ptm.end - ptm
  return(list(results_beta = results_beta,
              results_gamma = results_gamma,
              results_z = results_z,
              results_y = results_y,
              results_zeta = results_zeta,
              results_alpha = results_alpha,
              results_sig = results_sig,
              results_pi = results_pi,
              results_pi_SE = results_pi_SE,
              results_mu = results_mu,
              results_mu.hat = results_mu.hat,
              results_ey.hat = results_ey.hat,
              results_z_true = results_z_true,
              results_gradients = results_gradients,
              results_ppml_df = results_ppml_df,
              results_mu_ppml = results_mu_ppml,
              coef_ppml = coef_ppml,
              results_df = results_df,
              true_alpha = alpha,
              true_zeta = zeta,
              true_beta = beta,
              true_sig = sig,
              true_pi = pi.par,
              true_gamma = gamma,
              true_X = X,
              true_Xprime = Xprime_mL,
              true_z = z,
              time.taken = time.taken,
              params1 = params1,
              param_names = param_names
  ))
  as.list(environment())
}

# Raw results ------------------
# get raw results for the first set of parameters (toy example)
raw_out = mc_zaga_xij_bij(par = params1, seed = 1)

# Compute performance ------------------
getAll(raw_out, params1)
m = p*(p-1)
beta = t(true_beta)
alpha = true_alpha
gamma = true_gamma

edge_names <- c()
for (i in 1:p) {
  for (j in 1:p) {
    if (i != j) {
      edge_names <- c(edge_names, paste0(i, j))
    }
  }
}

names_alpha <- paste0("Alpha ", 0:(q))

rownames(alpha) <- c(names_alpha)
colnames(alpha) <- c(edge_names)

names_beta <- paste0("Beta", 1:L)

rownames(beta) <- c(names_beta)
colnames(beta) <- c(edge_names)

### Initialize the bias and RMSE matrix BETA
bias_beta = rmse_beta = matrix(0, nrow = nrow(beta), ncol = ncol(beta))

# Subtract the true values from each simulated estimate and accumulate the errors

for (l in 1:L){
  for (ij in 1:m){
    true_b = beta[l,ij]
    beta.est = sapply(1:MC.size, function (x) t(results_beta[[x]])[l,ij])
    bias_beta[l,ij] = mean(beta.est) - true_b
    rmse_beta[l,ij] = caret::RMSE(pred = beta.est, obs = true_b)
  }
}

colnames(bias_beta) <- c(edge_names)
rownames(bias_beta) <-c(names_beta)
colnames(rmse_beta) <- c(edge_names)
rownames(rmse_beta) <-c(names_beta)

### Initialize the bias matrix ALPHA
bias_alpha <- matrix(0, nrow = nrow(alpha), ncol = ncol(alpha))

# Subtract the true values from each simulated estimate and accumulate the errors
for (mat in results_alpha) {
  # bias_alpha <- bias_alpha + abs(mat - alpha)                   # abs(bias)
  bias_alpha <- bias_alpha + (mat - alpha)                        # bias
}

bias_alpha <- bias_alpha / MC.size
bias_alpha <- t(bias_alpha)

colnames(bias_alpha) <- c(names_alpha)
rownames(bias_alpha) <- c(edge_names)

# Calculate RMSE beta
rmse_alpha <- matrix(0, nrow = nrow(alpha), ncol = ncol(alpha))

# Subtract the true values from each simulated estimate and accumulate the errors
for (mat in results_alpha) {
  rmse_alpha <- rmse_alpha + ((mat - alpha)^2)
}

rmse_alpha <- sqrt(rmse_alpha / MC.size)
rmse_alpha <- t(rmse_alpha)
colnames(rmse_alpha) <- c(names_alpha)
rownames(rmse_alpha) <- c(edge_names)

# Extract only the "gamma_est" matrices from each simulation
gamma_est_matrices <- results_gamma

# Calculate the mean across all simulations for each element in the matrices
# `Reduce` combines all matrices by adding them element-wise, then we divide by the number of simulations
gamma_est_mean <- Reduce("+", gamma_est_matrices) / MC.size

# Empty matrix to accumulate squared deviations
gamma_est_variance <- matrix(0, nrow = nrow(gamma_est_mean), ncol = ncol(gamma_est_mean))


# Sum the squared deviations from the mean for each matrix
for (mat in gamma_est_matrices) {
  gamma_est_variance <- gamma_est_variance + (mat - gamma_est_mean)^2
}

# Divide by the number of simulations to get the variance
gamma_est_variance <- gamma_est_variance / (MC.size-1)

### Initialize the bias matrix BETA
bias_gamma <- matrix(0, nrow = 1, ncol = ncol(gamma))

# Subtract the true values from each simulated estimate and accumulate the errors
for (mat in gamma_est_matrices) {
  bias_gamma <- bias_gamma + abs(mat - gamma)
}

bias_gamma <- bias_gamma / MC.size


# Calculate RMSE beta
rmse_gamma <- matrix(0, nrow = 1, ncol = ncol(gamma))

# Subtract the true values from each simulated estimate and accumulate the errors
for (mat in gamma_est_matrices) {
  rmse_gamma <- rmse_gamma + ((mat - gamma)^2)
}

rmse_gamma <- sqrt(rmse_gamma / MC.size)

sig = true_sig
rmse_sig = caret::RMSE(pred = results_sig, obs = true_sig)

pi.par = true_pi
rmse_pi = numeric(length(pi.par))

for (ij in 1:m){
  pi.true1 = pi.par[ij]
  pi.est = sapply(1:MC.size, function(x) results_pi[[x]][ij])
  rmse_pi[ij] = caret::RMSE(pred = pi.true1, obs = pi.est)
}

# Sim Plots -----------
## Density ------------------

### alpha, beta ------------------
lam1 = results_mu[[1]]

mycolors = unlist(sapply(1:4, function(i) rep(cols[i], c(m-1,m-2,m,Lprime)[i])))

# 1) Build a single data.frame of errors for all edges & both alphas & betas
alpha_df <- bind_rows(
  # all alpha0 errors
  lapply(2:m, function(i) {  # exclude the fixed values
    data.frame(
      Error     = sapply(results_alpha, `[`, 1, i) - alpha[1, i],
      Edge      = edge_names[i],
      Parameter = "alpha0"
    )
  }),
  # all alpha1 errors
  lapply(3:m, function(i) {  # exclude the fixed values
    data.frame(
      Error     = sapply(results_alpha, `[`, 2, i) - alpha[2, i],
      Edge      = edge_names[i],
      Parameter = "alpha1"
    )
  }),
  # all beta errors
  lapply(seq_len(m), function(i) {
    data.frame(
      Error     = sapply(results_beta, `[`, i, 1) - true_beta[i,],
      Edge      = edge_names[i],
      Parameter = "beta"
    )
  }) 
)

alpha_df <- alpha_df %>%
  mutate(
    Parameter = factor(Parameter,
                       levels = c("alpha0","alpha1","beta"),
                       labels = c("alpha[0]", "alpha[(2)]","beta"))
  )

p1 <- ggplot(alpha_df, aes(x = Error,
                           colour = Edge,          # edge-specific lines
                           fill   = Parameter)) +  # Parameter-specific fill
  geom_density(alpha    = 0.4,          # more transparent for stacking
               adjust   = 2,            # smoothing
               position = "identity") + # draw them on top of each other
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Parameter,
             nrow          = 1,
             scales        = "free_x",
             strip.position = "top",
             labeller      = label_parsed) +
  scale_fill_manual(name = NULL,
                    values = c(
                      "alpha[0]"   = cols[1],
                      "alpha[(2)]" = cols[2],
                      "beta"       = cols[3]
                    )) +
  scale_colour_grey(start = 0.2, end = 0.6, guide = "none") +
  labs(
    x     = expression(hat(theta) - theta),
    y     = "Density"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.text       = element_text(size = 19, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    text = element_text(size = 20),
    legend.position  = "none",
    axis.text.x = element_text(
      angle = 45, vjust = 1, hjust = 1
    ),
    axis.ticks.length = unit(6, "pt"),
    plot.margin = margin(8, 12, 18, 12)  # extra bottom margin
  ) +
  scale_x_continuous(limits = c(-0.4, 0.4))

print(p1)


### gamma, sigma, pi ------------------

mycolors = unlist(sapply(1:4, function(i) rep(cols[i], c(m-1,m-2,m,Lprime)[i])))

# 1) Build a single data.frame of errors for all edges & both alphas & betas
alpha_df <- bind_rows(
  # all alpha0 errors
  lapply(1:m, function(i) {    # exclude the fixed values
    data.frame(
      Error     = sapply(results_pi, `[`, i) - true_pi[i],
      Edge      = edge_names[i],
      Parameter = "pi"
    )
  }),
  # all gamma errors
  lapply(Lprime, function(i) {  # exclude the fixed values
    data.frame(
      Error     = sapply(results_gamma, `[`, 1, i) - gamma[1, i],
      Edge      = edge_names[i],
      Parameter = "gamma"
    )
  }),
  # all sigma errors
  data.frame(
    Error     = results_sig - true_sig,
    Edge      = edge_names[i],
    Parameter = "sig"
  )
)

alpha_df <- alpha_df %>%
  mutate(
    Parameter = factor(Parameter,
                       levels = c("pi","gamma","sig"),
                       labels = c("pi", "gamma","sigma"))
  )

p1 <- ggplot(alpha_df, aes(x = Error,
                           colour = Edge,          # edge-specific lines
                           fill   = Parameter)) +  # Parameter-specific fill
  geom_density(alpha    = 0.4,          # more transparent for stacking
               adjust   = 2,            # smoothing
               position = "identity") + # draw them on top of each other
  geom_vline(xintercept = 0, linetype = "dashed") +
  facet_wrap(~ Parameter,
             nrow          = 1,
             scales        = "free",
             strip.position = "top",
             labeller      = label_parsed) +
  scale_fill_manual(name = NULL,
                    values = c(
                      "pi"   = cols[1],
                      "gamma" = cols[2],
                      "sigma"       = cols[3]
                    )) +
  scale_colour_grey(start = 0.2, end = 0.6, guide = "none") +
  labs(
    #title = "Densities of (Estimate − True) by Edge and Parameter",
    x     = expression(hat(theta) - theta),
    y     = "Density"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text       = element_text(size = 19, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    text = element_text(size = 20),
    legend.position  = "none",
    axis.text.x = element_text(
      angle = 45, vjust = 1, hjust = 1
    ),
    axis.ticks.length = unit(6, "pt"),
    plot.margin = margin(8, 12, 18, 12)  # extra bottom margin
  ) +
  # each parameters gets their own x limits, otherwise plots are unreadable
  facetted_pos_scales(
    x = list(
      Parameter == "pi"    ~ scale_x_continuous(limits = c(-0.15, 0.15)),
      Parameter == "gamma" ~ scale_x_continuous(limits = c(-0.04, 0.04)),
      Parameter == "sigma" ~ scale_x_continuous(limits = c(-0.01, 0.01))
    )) 

print(p1)

# PPMLE -----------
df_plot = do.call(rbind.data.frame,raw_out$results_ppml_df)

df_diag <- df_plot %>%
  mutate(
    dist_bin = ntile(x1, 20),
    best_bin = ntile(w1, 20),
    is_pos = y > 0
  )

df_cmp = df_diag %>%
  mutate(
    obs_zero = as.integer(y == 0),
    pzero_glamle = pi.hat,
    pzero_ppml_working = exp(-pmin(mu_ppml, 700))
  )

## Zero part (Dimension 1) -----------
layer_chosen = "k1"

zero_cal = bind_rows(
  df_cmp %>%
    transmute(model = "ZAGA-GLAMLE", p_zero = pzero_glamle, obs_zero = pi, edge = edge, layer = layer),
  df_cmp %>%
    transmute(model = "Poisson-PPMLE", p_zero = pzero_ppml_working, obs_zero = pi, edge = edge, layer = layer)
) %>% 
  filter(layer == layer_chosen) %>%
  group_by(model,edge) %>% 
  summarise(
    pred_zero = mean(p_zero, na.rm = TRUE),
    obs_zero = mean(obs_zero, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

cols <- c(
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#009E73", # green
  "#CC79A7", # reddish purple
  "#F0E442", # yellow
  "#56B4E9", # sky blue
  "#E69F00", # orange
  "#000000"  # black
)

p10 = ggplot(zero_cal, aes(x = pred_zero, y = obs_zero, shape = model, color = model)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.7, linetype = "dashed", color = cols[2]) +
  geom_point() +
  geom_line() +
  coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
  scale_color_manual(values = cols[c(3,1)]) + 
  theme_minimal(base_size = 25) +
  theme(legend.position = "top") +
  labs(
    x = "Mean estimated zero probability",
    y = bquote("True"~pi[ij]),
    shape = "",
    color = ""
  )

print(p10)


## Pos part (Dimension 2) -----------
df_cmp = df_cmp %>%
  mutate(
    mu_ppml_pos_working = mu_ppml / (1 - exp(-mu_ppml))
  )

pos_cal = bind_rows(
  df_cmp %>%
    filter(y>0, !is.na(mu.hat), layer == layer_chosen) %>%
    transmute(model = "ZAGA positive part", pred = mu.hat, obs = true_mu, edge = edge, layer = layer),
  df_cmp %>%
    filter(y > 0, layer == layer_chosen) %>%
    transmute(model = "PPML working positive mean", pred = mu_ppml_pos_working, obs = true_mu, edge = edge, layer = layer)
) %>%
  group_by(model, edge, layer) %>%
  summarise(
    pred_mean = mean(pred, na.rm = TRUE),
    obs_mean = mean(obs, na.rm = TRUE),
    n = n(),
    .groups = "drop"
  )

p11 = ggplot(pos_cal, aes(x = pred_mean, y = obs_mean, shape = model, color = model)) +
  geom_abline(intercept = 0, slope = 1, linewidth = 0.7, linetype = "dashed", color = cols[2]) +
  geom_point() +
  geom_line() +
  scale_color_manual(values = cols[c(3,1)]) + 
  theme_minimal(base_size = 25) +
  theme(legend.position = "none") +
  labs(
    x = "Mean estimated positive y",
    y = bquote("True"~mu[ij])
  )

p11

## MAE (Dimension 3) -----------

# compute RMSE and MAE layer-wise
df_rmse_by_layer = df_plot %>% 
  group_by(layer) %>% 
  summarise(
    rmse_zaga = Metrics::rmse(actual = ey, predicted = ey.hat),
    rmse_ppml = Metrics::rmse(actual = ey, predicted = mu_ppml),
    mae_ppml = Metrics::mae(actual = ey, predicted = mu_ppml),
    mae_zaga = Metrics::mae(actual = ey, predicted = ey.hat))

loss_long = df_rmse_by_layer %>%
  select(-c(rmse_zaga,rmse_ppml)) %>%
  pivot_longer(
    cols = -c(layer),
    names_to = c("metric", "model"),
    names_pattern = "(.*)_(zaga|ppml)",
    values_to = "loss"
  ) %>%
  mutate(
    model = recode(model, zaga = "ZAGA-GLAMLE", ppml = "PPMLE"),
    metric = recode(
      metric,
      mae = "Mean absolute error"
    )
  )

p_mae = ggplot(loss_long, aes(x = model, y = loss, color = model)) +
  geom_boxplot(outlier.size = 1) + # outlier.shape = NA to remove outliers
  coord_cartesian(ylim = quantile(loss_long$loss, c(0.01, 0.99), na.rm = TRUE)) +
  scale_color_manual(values = cols[c(3,1)]) +
  theme_minimal(base_size = 25) +
  theme(legend.position = "none") +
  labs(
    x = NULL,
    y = "MAE"
  )

p_mae
