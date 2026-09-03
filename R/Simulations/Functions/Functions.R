"
For all functions:

nv = p:       number of nodes
m:            number of edges in a directed graph, m = nv(nv-1)
q:            number of latent variables
a=alpha:      factor loadings
z=u:          latent variables
X:            layer-dependent covariates
L:            number of layer-dependent covariates
Xprime_ml(W): layer-INdependent covariates
Lprime:       number of layer-INdependent covariates
eta:          linear predictor, see Models in paper
mu:           mean, see Models in paper
"

# Poisson -------------------
# GLAMLE Poisson functions

# Generate graph
graph_gen_po_x_ij <- function(a, b, n, X = NULL, z = NULL, seed = 1) {
  # a: (q+1) x m            factor loadings / latent weights (alpha)
  # b: m x L or 1 x L       coefficients for observed dyadic covariates X
  # X: n x m x L            observed dyadic covariates (layer-specific)
  # z (u): n x (q+1)            [1, u_k^T] rows; if NULL it is generated with u_k ~ N(0, I_q)
  # n:                      number of layers (K)
  
  set.seed(seed)
  m <- ncol(a)
  q <- nrow(a) - 1L
  L <- if (is.null(X)) 0L else dim(X)[3]
  
  if (is.null(z)) {
    z.2 <- MASS::mvrnorm(n, rep(0, q), diag(q))
    z   <- cbind(1, z.2)
  } else {
    stopifnot(ncol(z) == q + 1L)
    z.2 <- z[, 2:ncol(z), drop = FALSE]
  }
  
  # Linear predictor
  eta <- z %*% a  # n x m
  if (L > 0L) {
    stopifnot(is.array(X), length(dim(X)) == 3L, dim(X)[1] == n, dim(X)[2] == m)
    if (all(is.matrix(b), nrow(b) == m, ncol(b) == L)){
      for (l in seq_len(L)) {
        eta <- eta + X[, , l] * matrix(b[, l], nrow = n, ncol = m, byrow = TRUE)
      }
    } else if (all(is.matrix(b), nrow(b) == 1, ncol(b) == L)){
      for (l in seq_len(L)) {
        eta <- eta + X[, , l] * b[, l]
      }
    }
    
  }
  
  Lam <- exp(eta) # Poisson means
  y   <- matrix(rpois(n * m, lambda = as.vector(Lam)), nrow = n, ncol = m)
  
  list(y = y, Lam = Lam, lvs = z.2, X = X)
}

## Po - bxij -------------------
# Core nll
glamle_nll_po_x_ij <- function(parms, data_input) {
  # parms are provided by RTMB/MakeADFun
  # data_input is a list containing y and (optionally) x
  getAll(data_input, parms, warn = FALSE)
  
  K <- nrow(y)
  m <- ncol(y)
  
  # Robust handling of x (may be missing or a degenerate object when L=0)
  has_x <- exists("x") && !is.null(x) && is.array(x) && length(dim(x)) == 3L && dim(x)[3] > 0L
  L <- if (has_x) dim(x)[3] else 0L
  
  # Rebuild alpha
  alpha <- .rebuild_alpha(q = q, m = m, intercept = intercept, lambda = lambda)
  
  # Linear predictor
  #u_mat <- matrix(u, nrow = K, ncol = q)  # works whether u is vector or matrix
  z1  <- cbind(1, u)            # K x (q+1)
  eta <- z1 %*% alpha           # K x m
  
  if (L > 0L) {
    stopifnot(dim(x)[1] == K, dim(x)[2] == m)
    # stopifnot(is.matrix(beta), nrow(beta) == m, ncol(beta) == L)
    
    if (all(is.matrix(beta), nrow(beta) == m, ncol(beta) == L)) {
      # Efficiently add sum_l x[,,l] * beta[,l]
      for (l in seq_len(L)) {
        eta <- eta + x[, , l] * matrix(beta[, l], nrow = K, ncol = m, byrow = TRUE)
      }
    } else if (all(is.matrix(beta), nrow(beta) == 1, ncol(beta) == L)){
      for (l in seq_len(L)) {
        eta <- eta + x[, , l] * beta[, l]
      }
    }
    
  }
  
  # Poisson log-likelihood (explicit form avoids calling dpois inside AD)
  mu <- exp(eta)
  ll_y <- sum(y * eta - mu - lgamma(y + 1))
  
  # *** CRITICAL: latent Gaussian prior must be included ***
  # Assumption A2 in the paper: u_k ~ N(0, I_q).
  ll_u <- sum(dnorm(u, mean = 0, sd = 1, log = TRUE))
  
  nll <- -(ll_y + ll_u)
  nll
}

# ZAGA -------------------

# generate true parameters for ZAGA
gen_params <- function(par, bij = T, xij = T){
  getAll(par, warn = T)
  out = list()
  
  m = p*(p-1)   # number of edges for direct graph
  
  # Draw alphas from Unif(0,1)
  set.seed(seed)
  alpha.temp <- matrix(runif(((q+1)*m),0,1), q+1, m ) # factor loadings
  
  # Apply constraints on alpha_(2)
  alpha.temp[lower.tri(alpha.temp)] <- 0
  diag(alpha.temp) = 1
  
  set.seed(seed)
  # Generate m beta coefficients for the model: beta is m x L matrix
  # Keep X the same; generate all x_ij,(k)
  if (L != 0){
    if (bij){
      beta <- matrix(runif((L)*m,-0.5,0.5),m,L)
    } else {
      beta <- matrix(runif((L)*1,-0.5,0.5),1,L)
    }
    
    # Generate covariates if L != 0
    if (xij){
      X <- array(runif(K*m*L, min = -1, max = 1), dim = c(K, m, L))
    } else {
      # add if needed to generate X^k 
      print("Add generation of x^k!")
    }
    
  } else {
    beta <- NULL; X <- NULL
  }
  
  # x': generate all x_ij,(l)
  set.seed(seed)
  # Generate L' gamma coefficients for the model: gamma is 1 x L' matrix
  if (Lprime != 0){
    gamma <- matrix(runif((Lprime)*1,-1,1),1,Lprime)
    
    # Generate covariates if L' != 0, constant for all K
    X_prime <- replicate(Lprime, {
      x <- runif(m, min = -1, max = 1)    # an m-vector of Unif(-1,1)
      outer(rep(1, K), x)                 # K x m matrix with columns = x
    }, simplify = "array")
  } else {
    gamma <- NULL; X_prime <- NULL
  }
  
  # X' is m x L' matrix
  Xprime_mL = as.matrix(X_prime[1,,])
  
  # Constraints: w'a0 = 0
  # Build H spanning Null(X'^T)
  qrX  <- qr(Xprime_mL)                        # m x L
  Q    <- qr.Q(qrX, complete = TRUE)           # m x m
  rank <- qrX$rank
  H    <- Q[, (rank+1):ncol(Q), drop = FALSE]  # m x (m-rank), Null(X'^T)
  
  # reconstruct alpha0
  set.seed(seed)
  r = m - rank
  zeta = rnorm(n = r)
  
  alpha0 = H %*% zeta
  alpha = rbind(c(alpha0), alpha.temp[-1,])
  
  set.seed(seed)
  
  # generate latent variable from multivariate normal distribution
  z.2 <- mvrnorm(K,rep(0,q),diag(q))
  z <- cbind(rep(1,K),z.2)
  
  # sigma parameter (fixed; change here however needed)
  sig <- 0.3
  # generate pi from truncated normal; change to however needed
  pi.par = truncnorm::rtruncnorm(n = m, a = 0.01, b = 0.8, mean = 0.5, sd = 0.12) %>% pmin(pmax(., 0.01), 0.6)
  
  # vector of true params theta0
  theta0 = c(zeta, a_to_start(alpha), z.2, beta, gamma, sig, pi.par)
  
  out$m = m
  out$beta = beta
  out$gamma = gamma
  out$X_prime = X_prime
  out$Xprime_mL = Xprime_mL
  out$X = X
  out$alpha = alpha
  out$z = z
  out$z.2 = z.2
  out$theta0 = theta0
  out$zeta = zeta
  out$H = H
  out$r = r
  out$m = m
  out$pi.par = pi.par
  out$sig = sig
  
  return(out)
}

# Draw samples from the zero-inflated gamma (ZAGA) distribution with the
# Rigby/Stasinopoulos parameterization used in GAMLSS (mu, sigma).
# - pi    : zero-inflation probability (0 <= pi <= 1)
# - mu    : mean of the gamma component (>0)
# - sigma : coefficient of variation (>0); shape = 1 / sigma^2, scale = mu * sigma^2
# Scalars are recycled to length n.

rzig_rigby <- function(n, pi, mu, sigma) {
  pi    <- rep_len(pi, n)
  mu    <- rep_len(mu, n)
  sigma <- rep_len(sigma, n)
  
  if (any(pi < 0 | pi > 1)) stop("pi must be in [0, 1].")
  if (any(mu <= 0)) stop("mu must be > 0.")
  if (any(sigma <= 0)) stop("sigma must be > 0.")
  
  shape <- 1 / (sigma^2)
  scale <- mu * (sigma^2)
  
  is_zero <- rbinom(n, size = 1, prob = pi)
  draws   <- rgamma(n, shape = shape, scale = scale)
  draws[is_zero == 1] <- 0
  
  draws
}

graph_gen_zaga_x_ij_mixK <- function(a,b,g,n,X=NULL,z=NULL,Xprime=NULL,seed=1,pi.par,sig, bij = T, xij = T){
  
  ### INPUTS
  # a [(q+1) x m]:          latent variable loadings (alpha)
  # b [L x m]:              coefficients for observed covariates X, layer dependent
  # g [1 x L']:             coefficients for observed covariates X, layer INdependent
  # X [K x m x L]:          observed covariates: setting with x_ij,k
  # X'[m x L]:              observed covariates: setting with w_ij
  # z [K x (q+1)]:          q latent variables and a vector of ones
  # seed:                   set seed for reproducibility
  
  ### OUTPUTS
  # y [K x m]:              simulated Poisson response variables
  # Mat [K x m]:            lambdas for each response variable
  # z.2 [K x q]:            simulated latent variables if none were provided
  # X [K x m x L]:          covariates (same as input)
  # X'[m x L]:              covariates (same as input)
  
  set.seed(seed)
  m = ncol(a)
  
  L = ifelse(is.null(X), 0, dim(X)[3])
  Lprime = ifelse(is.null(Xprime), 0, ncol(Xprime))
  
  # if no latent variables were provided, generate them
  if (is.null(z)) {
    # generate latent variable via multivariate normal distribution
    z.2 <- mvrnorm(n,rep(0,q),diag(q))
    z <- cbind(rep(1,n),z.2)
  } else z.2 = z[,2:ncol(z)]
  
  # pre create matrices
  y <- matrix(data=NA, n,m)
  
  # Reconstruct eta
  # reconstruct the part with covariates dependent on k
  if (L != 0){
    if (xij && bij){
      # eta_ij^k = a_0ij + a_2ij'z^k + b_ij'x_ij^k + g'w_ij
      eta_mat2 <- t(vapply(seq_len(n), function(k) {
        Mk <- matrix(X[k, , , drop = FALSE], nrow = m, ncol = L)  # 1×m×L -> m×L
        rowSums(Mk * b)                                        # length m
      }, numeric(m)))
    } else if (xij && !bij){
      # eta_ij^k = a_0ij + a_2ij'z^k + b'x_ij^k + g'w_ij
      eta_mat2 <- t(vapply(seq_len(n), function(k) {
        Mk <- matrix(X[k, , , drop = FALSE], nrow = m, ncol = L)  # 1×m×L -> m×L
        Mk %*% t(b)                                        # length m
      }, numeric(m)))
    } else {
      # ADD what to do if xbij, xb
      print("Add xbij/xb part!")
    }
    
  } else {
    eta_mat2 = matrix(0, nrow = n, ncol = m)
  }
  
  if (Lprime != 0){
    eta_row3  <- as.numeric(Xprime %*% t(g))               # length m
    eta_mat3  <- matrix(eta_row3, nrow = n, ncol = m, byrow = TRUE)
  } else {
    eta_mat3 = matrix(0, nrow = n, ncol = m)
  }
  
  eta1 = z %*% a + eta_mat2 + eta_mat3
  
  mu1 = exp(eta1)
  
  for (k in 1:n){
    for (ij in 1:m){
      y[k, ij] = rzig_rigby(n = 1, pi = pi.par[ij], mu = mu1[k,ij], sigma = sig)
    }
  }
  
  return(list(y=y, mu1 = mu1, sig = sig, pi.par = pi.par, lvs = z.2, X=X, Xprime = Xprime))
}

# ZAGA nll
glamle_nll_zaga_x_ij_mixK_pos <- function(parms, data_input) {
  # extract all the parameters; data_input is the environment object containing all data
  getAll(data_input, parms, warn=FALSE)
  
  m = ncol(y); K = nrow(y)
  
  # Robust handling of x (may be missing or a degenerate object when L=0)
  has_x <- exists("x") && !is.null(x) && is.array(x) && length(dim(x)) == 3L && dim(x)[3] > 0L
  L <- if (has_x) dim(x)[3] else 0L
  
  # Robust handling of xprime (may be missing or a degenerate object when L'=0)
  has_xp <- exists("xprime") && !is.null(xprime) && is.matrix(xprime) && nrow(xprime) == m
  Lp <- if (has_xp) ncol(xprime) else 0L
  
  # Rebuild alpha
  alpha <- .rebuild_alpha_zeta(q = q, m = m, intercept = intercept, lambda = lambda, H = H)
  
  # Linear predictor
  z1  <- cbind(1, u)            # K x (q+1)
  eta <- z1 %*% alpha           # K x m
  
  if (L>0){
    stopifnot(dim(x)[1] == K, dim(x)[2] == m)
    if (all(dim(beta) == c(1,L))){
      beta = matrix(beta, nrow = m, ncol = L, byrow = T)
    }
    stopifnot(is.matrix(beta), nrow(beta) == m, ncol(beta) == L)
    
    # Efficiently add sum_l x[,,l] * beta[,l]
    for (l in seq_len(L)) {
      eta <- eta + x[, , l] * matrix(beta[, l], nrow = K, ncol = m, byrow = TRUE)
    }
  }
  
  eta_mat3 = matrix(0, nrow = K, ncol = m)
  
  if (Lp>0){
    # reconstruct the gamma *x' matrix
    eta_row3  <- as.numeric(xprime %*% as.numeric(gamma))               # length m
    eta_mat3  <- matrix(eta_row3, nrow = K, ncol = m, byrow = TRUE)
  }
  
  eta1 = eta + eta_mat3
  
  # clamp for numeric stability; comment if no clamping needed
  eta1 = .clamp(eta1)
  
  mu = exp(eta1)
  
  # stabilize logsig before exp (smooth clamp works with AD)
  logsig_s <- 10 * tanh(logsig / 10)
  sig <- exp(logsig_s)
  
  idxp = which(y>0)
  ypos = y[idxp]
  pr = dgamma(ypos, shape = 1/(sig^2),
              scale = mu[idxp]*(sig^2), log = T)
  
  nll <- - sum(sum(pr))
  ## add N(0,1) term to the log likelihood
  nll <- nll - sum(dnorm(u, 0, 1, log = TRUE))
  
  ## Return
  nll
}

# Real data (WTO) -------------------

# ZAGA nll with centering and scaling
glamle_nll_zaga_bx_ij_mixK_pos_cs <- function(parms, data_input) {
  # extract all the parameters; data_input is the environment object containing all data; centered and scaled X; centered and scaled W
  getAll(data_input, parms, warn=FALSE)
  
  m = ncol(y); K = nrow(y)
  
  L <- dim(x)[3]
  Lp <- ncol(ws)
  
  # Compute sigma^2
  sig <- exp(logsig)
  sig2 <- sig * sig
  
  # Rebuild alpha
  alpha <- .rebuild_alpha_zeta2(q = q, m = m, intercept = intercept, lambda = lambda)
  a2 = alpha[-1,,drop=F]
  
  gw <- as.numeric(ws %*% t(gamma))
  
  # Add the z ~ N(0,1)
  nll <- -sum(dnorm(u, 0, 1, log = TRUE))
  
  for (k in seq_len(K)) {
    eta_k <- alpha[1,] + as.numeric(t(a2) %*% u[k, ]) + gw
    
    for (l in seq_len(dim(x)[3])) {
      eta_k <- eta_k + beta[l] * (x[k, , l] + cx[l])
    }
    
    eta_k <- .clamp(eta_k)
    
    idx <- which(y[k, ] > 0)
    mu_k <- exp(eta_k[idx])
    nll <- nll - sum(dgamma(y[k, idx],
                            shape = 1 / sig2,
                            scale = mu_k * sig2,
                            log = TRUE))
    
  }
  
  ## Return
  nll
}

# compute mu with centering and scaling of X and W
.compute_mu_zaga_bx_ij_mixK_pos_cs <- function(parms, data_input) {
  # extract all the parameters; data_input is the environment object containing all data; centered and scaled X; centered and scaled 
  
  getAll(data_input, parms, warn=FALSE)
  
  m = ncol(y); K = nrow(y)
  
  L <- dim(x)[3]
  Lp <- ncol(ws)
  
  # Compute sigma^2
  sig <- exp(logsig)
  sig2 <- sig^2
  
  # Rebuild alpha
  alpha <- .rebuild_alpha_zeta2(q = q, m = m, intercept = intercept, lambda = lambda)
  a2 = alpha[-1,,drop=F]
  
  # Linear predictor
  z1  <- cbind(1, u)            # K x (q+1)
  eta <- z1 %*% alpha           # K x m
  
  if (L>0){
    stopifnot(dim(x)[1] == K, dim(x)[2] == m)
    if (all(dim(beta) == c(1,L))){
      beta = matrix(beta, nrow = m, ncol = L, byrow = T)
    }
    stopifnot(is.matrix(beta), nrow(beta) == m, ncol(beta) == L)
    
    # Efficiently add sum_l x[,,l] * beta[,l]
    for (l in seq_len(L)) {
      eta <- eta + (x[, , l] + cx)* matrix(beta[, l], nrow = K, ncol = m, byrow = TRUE)
    }
  }
  
  eta_mat3 = matrix(0, nrow = K, ncol = m)
  
  if (Lp>0){
    # reconstruct the g.til'w* matrix
    eta_row3  <- as.numeric(ws %*% as.numeric(gamma))               # length m
    eta_mat3  <- matrix(eta_row3, nrow = K, ncol = m, byrow = TRUE)
  }
  
  eta1 = eta + eta_mat3
  
  # clamp for numeric stability; comment if no clamping needed
  eta1 = .clamp(eta1)
  
  mu = exp(eta1)
  ## Return
  return(mu)
}

# Computes Dunn-Smyth residuals (used for starting values)
.ds_resid_zaga = function(y1, mu = NULL, sigma = NULL, nu, seed = 1) {
  set.seed(seed)
  
  K = length(y1)
  u = numeric(K)
  
  idx0 = which(y1 == 0)
  idxp = which(y1 > 0)
  
  # y=0: jump from 0 to nu
  # F(-1)=0; F(0)=nu; u~Unif(0,nu); r=Phi^(-1)(u)
  
  if (length(idx0) > 0) {
    u[idx0] = runif(length(idx0), min = 0, max = .clamp01(nu[idx0], 1e-12))
  }
  
  # y>0: continuous CDF
  # r=Phi^(-1)(F(y;mu,sigma,nu))
  if (length(idxp) > 0) {
    Fy = pZAGA(y1[idxp], mu = mu[idxp], sigma = sigma[idxp], nu = nu[idxp],
               lower.tail = TRUE, log.p = FALSE)
    u[idxp] = .clamp01(Fy, 1e-12)
  }
  
  qnorm(.clamp01(u, 1e-12))
}

# construct a y matrix from the data frame
make_y_matrix <- function(df_final, value_var,
                          year_sel = NULL, agg_sel = NULL,
                          fill_missing = 0,
                          reporter_id_col = "reporter_iso3a",
                          partner_id_col  = "partner_iso3a",
                          drop_self = TRUE) {
  if (missing(value_var) || !is.character(value_var) || length(value_var) != 1L) {
    stop("value_var must be a single character string naming the column to cast (e.g., 'best', 'mfn', 'value').")
  }
  
  df_final <- data.table::copy(df_final)
  data.table::setDT(df_final)
  
  needed <- c("year", "aggregation_level", "product_code",
              "reporter_code", "partner_code",
              reporter_id_col, partner_id_col, value_var)
  missing_cols <- setdiff(needed, names(df_final))
  if (length(missing_cols) > 0L) {
    stop("Missing columns in df_final: ", paste(missing_cols, collapse = ", "))
  }
  
  if (is.null(year_sel)) year_sel <- df_final[, sort(unique(year))][1]
  if (is.null(agg_sel))  agg_sel  <- df_final[, sort(unique(aggregation_level))][1]
  
  value_raw <- df_final[[value_var]]
  value_is_numeric <- is.numeric(value_raw)
  value_is_binary_factor <- is.factor(value_raw) &&
    all(stats::na.omit(as.character(value_raw)) %in% c("0", "1"))
  value_is_binary_logical <- is.logical(value_raw)
  
  if (!(value_is_numeric || value_is_binary_factor || value_is_binary_logical)) {
    stop(sprintf(
      "'%s' must be numeric, logical, or a factor with levels 0/1 to build a matrix.",
      value_var
    ))
  }
  
  # subset to one slice; keep ISO ids for naming
  df_sub <- df_final[
    year == year_sel & aggregation_level == agg_sel,
    .(product_code,
      reporter_id = get(reporter_id_col),
      partner_id  = get(partner_id_col),
      value       = if (value_is_numeric) {
        as.numeric(get(value_var))
      } else if (value_is_binary_logical) {
        as.numeric(get(value_var))
      } else {
        as.numeric(as.character(get(value_var)))
      })
  ]
  
  # optional: drop self pairs at the data level (recommended)
  if (drop_self) df_sub <- df_sub[reporter_id != partner_id]
  
  # sanity: ISO ids should exist
  if (anyNA(df_sub$reporter_id) || anyNA(df_sub$partner_id)) {
    stop("Some ", reporter_id_col, " or ", partner_id_col, " values are NA in the selected slice.")
  }
  
  # node set and fixed dyad ordering
  nodes <- sort(unique(c(df_sub$reporter_id, df_sub$partner_id)))
  
  dyads <- data.table::CJ(reporter_id = nodes, partner_id = nodes, unique = TRUE)
  if (drop_self) dyads <- dyads[reporter_id != partner_id]
  dyad_levels <- dyads[, paste(reporter_id, partner_id, sep = "-")]
  
  prod_levels <- sort(unique(df_sub$product_code))
  
  # build dyad id with fixed ordering
  df_sub[, dyad := paste(partner_id, reporter_id, sep = "-")] # change this to df_sub[, dyad := paste(reporter_id, partner_id, sep = "-")] to reverse order depending on how the data is coded
  df_sub[, product_code := factor(product_code, levels = prod_levels)]
  df_sub[, dyad := factor(dyad, levels = dyad_levels)]
  
  # if duplicates exist per (product_code, dyad), aggregate before casting
  if (value_is_numeric || value_is_binary_factor || value_is_binary_logical){
    df_sub <- df_sub[, .(value = sum(value, na.rm = TRUE)), by = .(product_code, dyad)]
  }
  
  
  # wide reshape: product_code x dyad
  wide <- data.table::dcast(
    df_sub,
    product_code ~ dyad,
    value.var = "value",
    fill = fill_missing,
    drop = FALSE
  )
  
  y <- as.matrix(wide[, -"product_code"])
  rownames(y) <- as.character(wide$product_code)
  y
}

# generate starting values
.generate_start_cs <- function(par, s = 1, family = 2, min1 = -0.05, y = NULL, X1_total = NULL, W_star = NULL, eps = 1e-3, min.pos = 5, bij = T, xij = T, mu_x = NULL, S_x = NULL){
  
  ## ROUND 1 
  
  # compute starting values
  K <- nrow(y); m <- ncol(y)
  has_x <- exists("X1_total") && !is.null(X1_total) && is.matrix(X1_total)
  L <- if (has_x) ncol(X1_total) else 0
  Xt <- if (has_x) aperm(array(t(X1_total), dim = c(L, K, m)), c(2, 3, 1)) else 0L
  has_xprime <- exists("W_star") && !is.null(W_star) && is.matrix(W_star) && nrow(W_star) == m
  Lprime <- if (has_xprime) ncol(W_star) else 0L
  
  ds.residuals = vector("list", m)
  s0 = rep(NA_real_, m)
  beta_mat = matrix(NA, nrow = m, ncol = L)
  sigma_vec_store = numeric(m)
  
  set.seed(s)  # seed once
  parameters1 = list()
  
  if (family == 4){
    ### situation 1: all 0 
    idx0 = which(colSums(y)==0)
    
    if (length(idx0) > 0){
      y1 = y[,idx0[1]]
      # u~Unif(0,nu); r = Phi^(-1)(u)
      # all-zero edge: skip gamlss
      for (ij in idx0){
        ds.residuals[[ij]] = .ds_resid_zaga(y1 = y1, mu = NULL, sigma = NULL, nu = rep(1, K), seed = 1)
        print(paste0("edge ", ij," all 0"))
      }
      # keep s0 as NA so it does not contaminate gamma_hat / a0
    }
    
    ### situation 2: too few positives
    # indices where number of positive entries per edge < min.pos
    idx1 = which((K-colSums(y==0)) < min.pos & (K-colSums(y==0)) > 0)
    
    # indices where 1+ obs are positive
    idx.pos = setdiff(1:m, idx0)
    
    # non-problematic indices 
    idx.good = setdiff(idx.pos,idx1)
    
    if (length(idx.pos) > 0){
      
      for (ij in idx.pos){
        y1 <- y[, ij] #; x1 = X[,ij,]
        idx <- ((ij - 1) * K + 1):(ij * K)
        x1 <- X1_total[idx,,drop=F]
        
        model1 = try(gamlss(y1 ~ x1, family = ZAGA(), method = mixed(), sigma.formula = ~1), silent = TRUE)
        
        if (inherits(model1, "try-error")) {
          # what to do if model fails
        } else {
          # factor analysis
          mu = model1$mu.fv
          sigma_vec = model1$sigma.fv
          sigma = mean(sigma_vec)
          pi_vec = model1$nu.fv
          sigma_vec_store[ij] = sigma
          
          ds.residuals[[ij]] = .ds_resid_zaga(y1, mu = mu, sigma = sigma_vec, nu = pi_vec, seed = NULL)
          
          # extract other coefficients
          cf = coef(model1, what = "mu")
          s0[ij] = cf[1]
          if (L > 0) beta_mat[ij, ] = cf[-1]
        }
        
      }
      
      # back-transform edgewise intercepts
      cx <- if (!is.null(mu_x) && !is.null(S_x)) as.numeric(mu_x / diag(S_x)) else numeric(L)
      
      s01 <- s0
      if (L > 0 && length(cx) > 0) {
        s01 <- s0 - drop(beta_mat %*% cx)
      }
      
      good = which(is.finite(s01))
      bad  = setdiff(seq_len(m), good)
      
      # compute gamma.start using least squares:
      qrX = qr(W_star[good, , drop = FALSE])
      if (qrX$rank < ncol(W_star)) errorCondition("W* is not full rank.")
      
      gamma_hat = qr.coef(qrX, s01[good])
      
      # fill s0 for bad edges using covariate part only (so they do not pollute nullspace component)
      s0_filled = s01
      if (length(bad) > 0) {
        s0_filled[bad] = as.numeric(W_star[bad, , drop = FALSE] %*% gamma_hat)
      }
      
      # compute zeta.start:
      a0  <- s0_filled - as.numeric(W_star %*% gamma_hat)
      z0 <- a0[constr$F1]
      
      # compute alpha(2).start; u.start:
      # factor init: exclude bad edges from factanal, then pad loadings with zeros
      ds.res = do.call(cbind, ds.residuals)
      ds.res = as.matrix(ds.res)
      ds.res[!is.finite(ds.res)] = 0
      a.start = matrix(0, nrow = m, ncol = q)
      
      # try factor analysis on Dunn-Smyth residuals
      fa1 = try(factanal(ds.res, factors = q, scores = "regression"))
      
      if (inherits(fa1, "try-error")){
        ds.good <- do.call(cbind, ds.residuals[good])
        ds.good <- as.matrix(ds.good)
        fa2 = try(factanal(ds.good, factors = q, scores = "regression"))
        if (inherits(fa2, "try-error")){
          print("factanal failed, do PCA")
          # PCA 
          temp.ds <- scale(ds.good, center = TRUE, scale = FALSE)
          sv <- svd(temp.ds)
          
          Uq <- sv$u[, 1:q, drop = FALSE]
          Dq <- diag(sv$d[1:q], nrow = q, ncol = q)
          Vq <- sv$v[, 1:q, drop = FALSE]
          
          # latent-variable starts: K x q
          scores <- sqrt(nrow(temp.ds)) * Uq
          
          # loading starts: m x q
          loadings_good <- Vq %*% Dq / sqrt(nrow(temp.ds))
          
          # a start
          a.start <- matrix(0, nrow = m, ncol = q)
          a.start[good, ] <- loadings_good
          # z(2) start
          z2.start <- clip(scores, -10, 10) # start for z2
        } else {
          a.start <- matrix(0, nrow = m, ncol = q)
          a.start[good, ] <- fa2$loadings
          # z(2) start
          z2.start <- clip(fa2$scores, -10, 10) # start for z2
          print("factanal on ds.good success")
        }
      } else {
        # a start
        a.start = matrix(pmax(pmin(fa1$loadings, 10), -10), nrow = m, ncol = q)
        # z(2) start
        z2.start <- clip(fa1$scores, -10, 10) # start for z2
      }
      
      intercept.start <- clip(t(as.matrix(a0)), -20, 20) # start a0
      a.c1 <- rbind(intercept.start, matrix(a.start,q,m))
      a.c = t(apply_custom_transformation2(t(a.c1))$transformed_A)

    }
    
    # record start for alpha(2) -> disregards intercept, which will be set later
    parameters1$lambda = a_to_start(a.c)
    
    # start for lv's (z/u)
    parameters1$u = z2.start
    # start for log(sigma)
    parameters1$logsig = log(mean(sigma_vec_store[is.finite(sigma_vec_store)]))
  }
  
  if (L > 0) {
    if (any(is.na(beta_mat))){
      
      if (!bij && xij){
        beta_mat = colMeans(beta_mat, na.rm = T)
      } else {
        for (l in seq_len(ncol(beta_mat))) {
          miss <- is.na(beta_mat[, l])
          if (any(miss)) {
            beta_mat[miss, l] <- mean(beta_mat[!miss, l], na.rm = TRUE)
          }
        }
        
      }
    } else if(!bij && xij){
      beta_mat = colMeans(beta_mat, na.rm = T)
    }
    parameters1$beta <- beta_mat
  } else {
    parameters1$beta <- matrix(0, nrow = 0, ncol = m)
  }
  
  if (Lprime >0){
    parameters1$intercept = as.numeric(z0)
    parameters1$gamma = matrix(gamma_hat, nrow = 1)
  } else {
    parameters1$intercept = intercept.start
    parameters1$gamma <- matrix(0, nrow = 1, ncol = Lprime)
  }
  print(paste0("Round 1 done"))

  
  ## ROUND 2 
  # Set up the data
  
  data_input <- list(y = y, q = q, ws = W_star, x = Xt, cx = cx)
  
  objr <- RTMB::MakeADFun(func = nll(glamle_nll_zaga_bx_ij_mixK_pos_cs, data_input),
                          parameters = parameters1, 
                          #random = "u", # does not integrate out lv's for speed - see paper 
                          silent = T,
                          inner.control = list(mgcmax = 1e+200, maxit = 1000)
  )
  
  out <- try(nlminb(
    start     = objr$par,
    objective = objr$fn,
    gradient  = objr$gr,
    control   = list(
      eval.max = 5000,
      iter.max = 2000,
      rel.tol  = 1e-8,
      x.tol    = 1e-10
    )
  ))
  
  if (inherits(out, "try-error")){
    # do nothing; return the parameters1 from ROUND 1
    print(paste0("Round 2 failed"))
  } else {
    # compute start
    par.hat = objr$env$parList(out$par)
    # update starting values from conditional log-likelihood
    parameters1 = par.hat
    print(paste0("Round 2 good"))
  }
  
  
  return(parameters1)
}

# computes scores matrix for the correct var-cov matrix
.compute_scores_matrix_cs_bxij <- function(y_mat, Xt, W_star_pos, cx, par.hat, p){
  # X should be X tilde!
  # Xprime_ml should be W_star ! 
  K = nrow(y_mat); m = ncol(y_mat); q = ncol(par.hat$u)
  L <- dim(Xt)[3]
  
  Lprime <- ncol(W_star_pos)
  
  scores  <- matrix(0, K, p)
  
  for(k in seq_len(K)){
    
    data_input_k <- list(y = y_mat[k,,drop=F], q = q, ws = W_star_pos, x = Xt[k,,,drop=F], cx = cx)
    
    parameters_k = list(
      lambda = par.hat$lambda,
      logsig = par.hat$logsig,
      beta = par.hat$beta,
      intercept = par.hat$intercept,
      gamma = par.hat$gamma,
      u = par.hat$u[k,,drop=F]
    )
    
    obj <- RTMB::MakeADFun(func = nll(glamle_nll_zaga_bx_ij_mixK_pos_cs, data_input_k),
                           parameters = parameters_k,
                           random = "u",
                           silent = T
    )
    
    scores[k, ] <- obj$gr(obj$par) 
  }
  return(scores)
}


# Helpers -------------------

# Helper: rebuild alpha (q+1) x m under the triangular-identification scheme used: for row k=2,...,q+1, alpha[k,1:(k-1)]=0, alpha[k,k]=1.
.rebuild_alpha <- function(q, m, intercept, lambda) {
  alpha <- matrix(0, nrow = q + 1L, ncol = m)
  
  # Intercept row.
  # NOTE: The first entry is fixed by construction, this can be omitted - see paper on more information (only factor loadings need to be constrained under this setup)
  alpha[1, ] <- c(1, intercept)
  
  idx <- 1L
  for (k in 2:(q + 1L)) {
    alpha[k, 1:(k - 1L)] <- 0
    alpha[k, k]          <- 1
    
    n_free <- m - k
    if (n_free > 0L) {
      end_idx <- idx + n_free - 1L
      alpha[k, (k + 1L):m] <- lambda[idx:end_idx]
      idx <- end_idx + 1L
    }
  }
  
  alpha
}

# Helper: rebuild alpha (q+1) x m under the triangular-identification scheme used in the presence of layer-INdependent covariates: a0 = H %*% zeta; for row k=2,...,q+1, alpha[k,1:(k-1)]=0, alpha[k,k]=1.
.rebuild_alpha_zeta <- function(q, m, intercept, lambda, H) {
  alpha <- matrix(0, nrow = q + 1L, ncol = m)
  
  # Intercept row.
  # a0 = H%*% zeta
  alpha[1, ] <- H %*% intercept
  
  idx <- 1L
  for (k in 2:(q + 1L)) {
    alpha[k, 1:(k - 1L)] <- 0
    alpha[k, k]          <- 1
    
    n_free <- m - k
    if (n_free > 0L) {
      end_idx <- idx + n_free - 1L
      alpha[k, (k + 1L):m] <- lambda[idx:end_idx]
      idx <- end_idx + 1L
    }
  }
  
  alpha
}

# Helper: rebuild alpha (q+1) x m under the triangular-identification scheme used in the presence of layer-INdependent covariates: different setup (tailored for wto, more efficient)
.rebuild_alpha_zeta2 <- function(q, m, intercept, lambda) {
  alpha <- RTMB::AD(matrix(0, nrow = q + 1L, ncol = m))
  
  # intercept row
  alpha[1, ] <- constr$G %*% intercept
  
  idx <- 1L
  for (k in 2:(q + 1L)) {
    alpha[k, k] <- 1
    
    n_free <- m - k
    if (n_free > 0L) {
      end_idx <- idx + n_free - 1L
      alpha[k, (k + 1L):m] <- lambda[idx:end_idx]
      idx <- end_idx + 1L
    }
  }
  
  alpha
}

# Extract unconstrained upper-triangular (excluding the diagonal) entries
# from rows 2:(q+1) of the (q+1) x m loading matrix a.
a_to_start <- function(a){
  q <- nrow(a) - 1L
  m <- ncol(a)
  out <- c()
  for(r in seq_len(q)){
    row <- r + 1L
    if(m > row) out <- c(out, a[row, (row+1L):m])
  }
  out
}

# from estimated parameters or starting values, create alpha matrix by adding the constraints to create a matrix
start_to_a <- function(a.start, intercept.start, q, m){
  # repopulate the alpha matrix
  newlam = matrix(0,q+1,m)
  
  # add the intercept
  intrcpt = c(1,intercept.start)
  newlam[1,] = intrcpt
  
  idx = 1
  # add the lambda
  for (k in 2:(q+1)){
    newlam[k, ] = c(rep(0, (k-1)),1,a.start[idx:(m-k)]) 
    idx = idx + (m-k+1)
  }
  return(newlam)
}

# extract estimated alpha.hat when layer-INdependent covariates are present
.get_alpha_hat <- function(obj=NULL, out, H, q, m) {
  if (is.null(obj)){
    a.hat <- start_to_a(a.start = out$lambda, intercept.start = out$intercept[1:(m-q)], q, m)
    a.hat[1,] = H %*% out$intercept
  } else {
    pl <- obj$env$parList(out$par)
    a.hat <- start_to_a(a.start = pl$lambda, intercept.start = out$par[1:(m-q)], q, m)
    # compute the intercept (alpha0)
    a.hat[1,] = H %*% pl$intercept
  }
  
  return(a.hat)
}

# Compute mu from the estimated parameters
.compute_mu <- function(par.hat, X, Xprime_mL, clamp = TRUE){
  "
  Inputs:
  par.hat:          Vector of estimated parameters
  X:                Layer-dependent covariates
  Xprime_mL (W):    Layer-INdependent covariates
  clamp:            For numerical stability; don't use if not needed
  "

  m = ncol(par.hat$alpha.hat); K = nrow(par.hat$u); q = ncol(par.hat$u)
  has_x <- exists("X") && !is.null(X) && is.array(X) && length(dim(X)) == 3L && dim(X)[3] > 0L
  L <- if (has_x) dim(X)[3] else 0L
  
  has_xprime <- exists("Xprime_mL") && !is.null(Xprime_mL) && is.matrix(Xprime_mL) && nrow(Xprime_mL) == m && ("gamma" %in% names(par.hat))
  Lprime <- if (has_xprime) ncol(Xprime_mL) else 0L
  
  Mat <- matrix(data = NA, nrow = K, ncol = m)
  
  # Reconstruct eta (linear predictor)
  if (L != 0){
    if (all(dim(par.hat$beta) == c(m,L))){
      eta_mat2 <- t(vapply(seq_len(K), function(k) {
        Mk <- matrix(X[k, , , drop = FALSE], nrow = m, ncol = L)  # 1×m×L -> m×L
        rowSums(Mk * par.hat$beta)                                        # length m
      }, numeric(m)))
    } else if (all(dim(par.hat$beta) == c(1,L))){
      eta_mat2 <- t(vapply(seq_len(K), function(k) {
        Mk <- matrix(X[k, , , drop = FALSE], nrow = m, ncol = L)  # 1×m×L -> m×L
        (Mk %*% t(par.hat$beta)) # mxL Lx1 -> length m
      }, numeric(m)))
    }
    
  } else {
    eta_mat2 = matrix(0, nrow = K, ncol = m)
  }
  
  if (Lprime != 0){
    eta_row3  <- as.numeric(Xprime_mL %*% t(par.hat$gamma))               # length m
    eta_mat3  <- matrix(eta_row3, nrow = K, ncol = m, byrow = TRUE)
  } else {
    eta_mat3 = matrix(0, nrow = K, ncol = m)
  }
  
  eta = cbind(1,par.hat$u) %*% par.hat$alpha.hat + eta_mat2 + eta_mat3
  
  if (clamp == TRUE) eta = .clamp(eta)
  return(exp(eta))
}

# helper for estimation
nll <- function(f, d) function(p) f(p,d)

# numeric
.clamp <- function(E, c_eta = 30) c_eta * tanh(E / c_eta)

# friendly colors
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

#  matrix transform
direct.matrix <- function(v,p){
  M <- matrix(0,p,p)
  for (i in 1:p) {
    
    M[i,-i] <- v[((i-1)*(p-1)+1):(i*(p-1))]
    
    
  }
  if (!is.null(colnames(v))){
    cols = colnames(v)
  }
  return(M)
}

# a more efficient way to impose constraint w'a0=0
make_a0_constraint <- function(W, tol = 1e-10) {
  W <- as.matrix(W)
  m <- nrow(W)
  p <- ncol(W)
  
  qr_col <- qr(W, tol = tol)
  r <- qr_col$rank
  
  if (r < p) {
    keep_cols <- sort(qr_col$pivot[seq_len(r)])
    W_red <- W[, keep_cols, drop = FALSE]
  } else {
    keep_cols <- seq_len(p)
    W_red <- W
  }
  
  r <- ncol(W_red)
  
  qr_row <- qr(t(W_red), tol = tol)
  C <- sort(qr_row$pivot[seq_len(r)])
  F1 <- setdiff(seq_len(m), C)
  
  WC <- W_red[C,  , drop = FALSE]
  WF <- W_red[F1, , drop = FALSE]
  
  if (qr(WC, tol = tol)$rank < r) {
    stop("Could not find a full-rank row block WC.")
  }
  
  # Full intercept is a0 = G %*% a0_free
  G <- matrix(0, nrow = m, ncol = length(F1))
  G[F1, ] <- diag(length(F1))
  G[C,  ] <- -solve(t(WC), t(WF))
  
  list(
    C = C,
    F1 = F1,
    keep_cols = keep_cols,
    W_red = W_red,
    G = G,
    rebuild_a0 = function(a0_free) drop(G %*% a0_free)
  )
}

# rotate the alpha to get the true alpha with the constraints
apply_custom_transformation2 <- function(A) {
  m <- nrow(A)
  q <- ncol(A)
  if(q != 2) stop("This implementation is for q = 2 only.")
  
  # Step 1: Compute t1 as the minimum-norm solution for t1 * A = (1, 0)
  AtA <- t(A) %*% A       # 2 x 2 matrix
  # Solve for t1: t1^T = A * (AtA)^{-1} * c, with c = (1,0)^T.
  t1_vec <- A %*% solve(AtA, c(1, 0))  
  # t1_vec <- A %*% ginv(AtA, c(1, 0))  
  t1 <- as.vector(t(t1_vec))  # row vector
  
  # Step 2: Compute t2 so that t2 * A gives (*, 1)
  a2 <- A[,2]
  norm_a2_sq <- sum(a2^2)
  t2 <- as.vector(a2 / norm_a2_sq)
  
  # Step 3: Construct T: set the first row to t1, second row to t2,
  # and complete with the standard basis for rows 3 through m.
  Tm <- diag(m)
  Tm[1, ] <- t1
  Tm[2, ] <- t2
  
  # (Optional: Check invertibility.)
  if (abs(det(Tm)) < 1e-8) {
    warning("T is nearly singular. Consider a different completion.")
  }
  
  transformed_A <- Tm %*% A
  
  # make sure it's 0 and not -1e-80 or something
  transformed_A[upper.tri(transformed_A)] <- round(transformed_A[upper.tri(transformed_A)])
  
  return(list(Tm = Tm, transformed_A = transformed_A))
}

# clip starting values of alpha/beta for numeric stability
clip <- function(x, lower = -10, upper = 10) {
  pmin(pmax(x, lower), upper)
}


.clamp01 = function(p, eps = 1e-10) pmin(pmax(p, eps), 1 - eps)