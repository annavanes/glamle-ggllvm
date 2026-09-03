pacman::p_load(TMB, mvabund, MASS, splm, matrixcalc, tidyverse, RTMB, vegan, Matrix, doParallel, gllvm, grDevices, circlize, colorspace, caret)
rm(list=ls())

#  ASSUMPTIONS ------------------
# Data generating process: Poisson with covariates x_ij^(k) and latent variables u_k, with log link. No zero-inflation.
"
GGLLVM / GLAMLE (Jiang, La Vecchia, Rastelli):
   Y_ij^(k) | u_k, x_ij^(k)  ~ Poisson(mu_ij^(k))
   log mu_ij^(k) = alpha_ij^T (1, u_k)^T + beta^T x_ij^(k)
   u_k ~ N(0, I_q)
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

# Set parameters for the simulation (toy example). Parameters used in the paper are p=10, q=1, K=500, L=1, Lprime=0, MC.size=100, family.ind=2 (Poisson), c_shift = log(5), seed = 2001.

params1 = list(
  p = 3,              # number of nodes
  q = 1,              # number of latent variables
  K = 0.5e2,          # number of layers
  L = 1,              # number of observable covariates (k dependent)
  Lprime = 0,         # number of observable covariates (k INdependent)
  MC.size = 2,        # number of simulations
  family.ind = 2,     # family: 1=Bernoulli, 2=Poisson, 3=ZIP, 4= ZAGA
  c_shift = log(5),   # shift for the intercepts to get a higher Poisson mean 
  seed = 2001         # seed for reproducibility      
)

# Simulation function:
# input: par (list of parameters), seed (for reproducibility)
mc_poi_xij_b <- function(par, seed) {
  getAll(par, warn = T)
  m = p*(p-1)         # number of edges for directed graph
  
  # Starting values (true+noise)
  starts = list()
  
  # Draw alphas from Unif(0,1)
  set.seed(seed)
  alpha.temp <- matrix(runif(((q+1)*m),0,1), q+1, m ) # factor loadings
  alpha.temp[1,] <- alpha.temp[1,] + c_shift # shift intercepts to a higher Poisson mean (not all zeros)
  
  # Identification constraints on alpha_(2): lower-triangular + ones on the diagonal
  alpha.temp[lower.tri(alpha.temp)] <- 0
  # Only rows 2:(q+1) are meant to be constrained; leaving row 1 free except for the *extra* choice alpha[1,1]=1 used below, constrained for sign.
  diag(alpha.temp) = 1
  
  alpha <- alpha_rot <- alpha.temp
  
  set.seed(seed)
  
  # generate latent variable from multivariate normal distribution
  z.2 <- MASS::mvrnorm(K, rep(0, q), diag(q))
  z <- cbind(1, z.2)
  
  set.seed(seed)
  # Generate m beta coefficients for the model: beta is m x L matrix
  if (L != 0) {
    beta <- matrix(runif((L)*1,-1,1),1,L)
    X <- array(runif(K * m * L, -1, 1), dim = c(K, m, L))
    
  } else {
    beta <- NULL
    X <- Xc <- NULL
  }
  
  # True parameter vector theta (for checks)
  theta0 <- c(alpha_rot[1, 2:m], a_to_start(alpha_rot), z.2, beta)
  
  # Starting values
  set.seed(1)
  min1 <- -0.05
  noise  <- runif(length(a_to_start(alpha_rot)), min = min1, max = -min1)
  noise2 <- runif(length(alpha_rot[1, 2:m]), min = min1, max = -min1)
  noise3 <- matrix(runif(length(z.2), min = min1, max = -min1), nrow = K, ncol = q)
  noise4 <- if (!is.null(beta)) matrix(runif(L, min = min1, max = -min1),1,L) else NULL
  
  parameters <- list(
    lambda = a_to_start(alpha_rot) + noise,
    intercept = alpha_rot[1, 2:m] + noise2,
    u = z.2 + noise3
  )
  
  if (L != 0) parameters$beta <- beta + noise4
  
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
  results_z <- vector("list", MC.size)
  results_alpha <- vector("list", MC.size)
  results_y <- vector("list", MC.size)
  results_mu <- vector("list", MC.size)
  results_mu.hat <- vector("list", MC.size)
  results_z_true <- vector("list", MC.size)
  results_gradients <- vector("list", MC.size)
  
  # Set up the progress bas
  pb = txtProgressBar(min = 0, max = MC.size, style = 3)
  
  # Run MC simulations
  for (sim in 1:MC.size) {
    
    tryCatch({
      # create Y, X and Z
      graph <- graph_gen_po_x_ij(a = alpha_rot, b = beta, n = K, X = X, z = z, seed =  seed + sim)
      
      y.sim <- graph$y
      results_y[[sim]] <- y.sim
      Lam.sim <- graph$Lam
      x.sim <- graph$X
      z2.sim <- graph$lvs
      z.sim <- cbind(rep(1,K),z2.sim)
      
      # Data for the objective
      # IMPORTANT: if L==0, pass x=NULL (not a 2D matrix) so dim(x)[3] is well-defined logic-wise.
      
      data_input <- list(y = y.sim, q = q)
      if (L != 0) data_input$x <- X
      
      # Optimize objective function
      
      objr <- RTMB::MakeADFun(
        func = nll(glamle_nll_po_x_ij, data_input),
        parameters = parameters,
        random = "u",
        silent = TRUE #, set to FALSE to see trace
      )
      
      # Optimize
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
      
      # Extract the gradient at the optimized parameters
      grad.sim <- objr$gr(out$par)
      
      # Extract parameter names
      param_names <- names(out$par)
      
      # Extract and store the theta.hat
      par.hat = objr$env$parList(out$par)
      par.hat$alpha.hat <- .rebuild_alpha(q = q, m = m,
                                          intercept = par.hat$intercept,
                                          lambda = par.hat$lambda)
      
      par.hat$lam.hat = .compute_mu(par.hat = par.hat, X = X, Xprime_mL = NULL, clamp = F)
      
      # Save results for this simulation
      results_beta[[sim]] <- par.hat$beta
      results_z[[sim]] <- cbind(1, par.hat$u)
      results_z_true[[sim]] <- z.sim
      results_alpha[[sim]] <- par.hat$alpha.hat
      results_mu[[sim]] <- Lam.sim
      results_mu.hat[[sim]] <- par.hat$lam.hat
      results_gradients[[sim]] <- grad.sim
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
              results_z = results_z,
              results_alpha = results_alpha,
              results_y = results_y,
              results_mu = results_mu,
              results_mu.hat = results_mu.hat,
              results_z_true = results_z_true,
              results_gradients = results_gradients,
              results_df = results_df,
              true_alpha = alpha,
              true_beta = beta,
              true_X = X,
              true_z = z,
              time.taken = time.taken,
              params1 = params1,
              param_names = param_names
  ))
  as.list(environment())
}

#  Raw results ------------------
# get raw results for the first set of parameters (toy example)
raw_out = mc_poi_xij_b(par = params1, seed = 1)

#  Compute performance ------------------
getAll(raw_out, params1)
beta = true_beta
names_beta <- paste0("Beta", 1:L)
rownames(beta) <- c(names_beta)

### Initialize the bias and RMSE matrix BETA
bias_beta = rmse_beta = matrix(0, nrow = nrow(beta), ncol = ncol(beta))

# Subtract the true values from each simulated estimate and accumulate the errors

for (l in 1:L){
  true_b = beta[l, ]
  beta.est = sapply(1:MC.size, function (x) t(results_beta[[x]])[l, ])
  bias_beta[l, ] = mean(beta.est) - true_b
  rmse_beta[l, ] = caret::RMSE(pred = beta.est, obs = true_b)
}

rownames(bias_beta) <-c(names_beta)
rownames(rmse_beta) <-c(names_beta)

#  Plots ------------------
##  Beta - density ------------------
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

# Violin plot of beta's
M <- do.call(rbind, results_beta)

beta_df = data.frame(Est = M[,1] - beta[1,1], Covariate = factor(paste("beta",1)))

rmse_df <- data.frame(
  Covariate = factor(paste("beta",1),
                     levels = levels(beta_df$Covariate)),
  rmse      = as.numeric(rmse_beta)
)

p1 <- ggplot(beta_df, aes(x = Est, fill = Covariate)) +
  geom_density(alpha    = 0.4,          # more transparent for stacking
               adjust   = 2,            # existing smoothing
               position = "identity") + # draw them on top of each other
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_x_continuous(limits = c(-0.009, 0.009)) +
  scale_fill_manual(values = c("beta 1" = cols[1])) +
  scale_color_grey(start = 0.2, end = 0.6, guide = "none") +
  labs(
    x     = expression(hat(beta) - beta),
    y     = "Density"
  ) +
  theme_minimal(base_size = 20) +
  theme(
    strip.text       = element_text(size = 14, face = "bold"),
    panel.grid.major = element_line(color = "grey90"),
    legend.position  = "none"
  )

print(p1)

##  Circle plot ------------------
ks <- 1:2       # take the first two layers as an example
nodes <- 1:p    # p is the number of countries (nodes) in the true dataset
y <- results_y[[1]] # choose simulated graph
  
par(
  mfrow = c(1, 2),
  mar = c(0.3, 0.3, 0.3, 0.3),   # inner margins (b l t r)
  oma = c(0, 0, 0, 0),
  xpd = NA,
  cex = 3
)
  
for (k in ks) {
  
  ## build the adjacency matrix for layer k
  M <- direct.matrix(y[k, , drop = FALSE], p)
  rownames(M) <- colnames(M) <- nodes
  
  circos.clear()              # reset everything
  
  chordDiagram(
    M,
    directional = 1,
    transparency = 0.5,
    annotationTrack = c("grid", "name"),   
    direction.type = c("arrows"),     # show arrow heads
    link.arr.type  = "big.arrow",     # arrow style: "triangle", "big.arrow", "arrow"
    link.arr.length = 0.04,           # arrow head length 
    link.arr.width  = 0.02
  )
}
