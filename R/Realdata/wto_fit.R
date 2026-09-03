# fit ZAGA GLAMLE on wto data

pacman::p_load(TMB, mvabund, MASS, splm, matrixcalc, tidyverse, RTMB, vegan, Matrix, doParallel, readxl, WDI, purrr, igraph, circlize, countrycode, data.table, scales, fixest, cplm, formatR, flextable, readr, gamlss, scales, gravity, mgcv, ggrepel)

rm(list=ls())

# Assumptions ------------------
"
Fit the model with ZAGA(mu_ij^k, sig, pi_ij)
eta_ij^k = a0_ij + a2_ij'z2^k + b'x_ij^k + g'w_ij
Only estimate the parameters where the colSums(y) >0 and there are at least 30 positive observations; do not estimate the rest!
"

# Data =========================================================

load("data/wto_all2022_agg3.RData")
load("data/countries_to_keep2_eda.RData")

files <- list.files(path = "R/Simulations/Functions",
                    pattern = "\\.R$",
                    full.names = TRUE)

# Quietly source the files
invisible(lapply(files, function(file) {
  suppressMessages(suppressWarnings(source(file)))
}))

# Look into it
colnames(wto)
head(wto)
dim(wto) # [1] 1170432      30

# prepare for analysis
setDT(wto)

# convert binary into factors
bin_cols <- names(wto)[vapply(wto, function(x) {
  is.numeric(x) &&
    all(na.omit(x) %in% c(0, 1)) &&
    length(unique(na.omit(x))) <= 2
}, logical(1))]

wto[, (bin_cols) := lapply(.SD, factor, levels = c(0, 1)),
    .SDcols = bin_cols]

# toy example here with only 5 countries; uncomment to see the one used in the paper 
# wto_subset = wto %>% filter(reporter_iso3a %in% countries_to_keep2, partner_iso3a %in% countries_to_keep2)

wto_subset = wto %>% filter(reporter_iso3a %in% c("EEC", "CHN", "USA","MEX", "NOR"), partner_iso3a %in% c("EEC", "CHN", "USA", "MEX", "NOR"))

nv = length(unique(wto_subset$reporter_iso3a)) # number of nodes (countries)

data=wto_subset

dim(data)
rm(wto,X,Xprime_mL,y) # remove temporary variables
gc()

# layer-dependent covariates (X)
formula = value ~ best
family = "ZAGA"

# layer-INdependent covariates (W)  
w = ~ dist + comlang_ethno + colony + contig
q=1

"
setup parameters:
se: whether to compute SE
sandwich.cov: whether to compute sandwich covariance matrix as in paper
center_x: whether to center layer-dependent covariates
scale_x: whether to scale layer-dependent covariates
center_xprime: whether to center layer-INdependent covariates
scale_xprime: whether to scale layer-INdependent covariates
min1:
eps:
min.pos: minimum number of positive observations for an edge to be included in the estimation
s: 
bij: should the coefficients for layer-dependent covariates be estimated for each edge (bij = T) or be shared across edges (bij = F)
xij: are layer-dependent covariates specific to each layer (xij = T) or shared across layers (xij = F)
"

se = F; sandwich.cov = T; center_x = T; scale_x = T; center_xprime = T; scale_xprime = T; min1 = -0.05; eps = 1e-3; min.pos = 30; s = 1; bij = F; xij = T

#Data Preprocessing=========================================================
# convert to data.table
df = as.data.table(data)
rm(data)
gc()

# extract the name of y
y_var = all.vars(formula)[1]
# construct K*m y matrix
y_mat = make_y_matrix(df, value_var = y_var, fill_missing = NA)
# extract number of layers (K) and edges (m)
K = nrow(y_mat); m = ncol(y_mat)
layers = rownames(y_mat) # goods traded
edges = colnames(y_mat) # directed country pairs: country i -> country j
idx_0 = which(colSums(y_mat) == 0)
idx_few_pos = which(colSums(apply(y_mat,2, function(x) x >0)) > 0 & colSums(apply(y_mat,2, function(x) x >0)) <= min.pos)
idx_pos = setdiff(1:m, union(idx_0, idx_few_pos))
m_pos = length(idx_pos)
edges_pos = edges[idx_pos]

# extract covariates
covariates = all.vars(formula)[-1]
# extract number of layer-dependent covariates (L; X)
L = length(covariates)

# Extract factors and continuous covariates
is_fac <- vapply(df, is.factor, logical(1))

factor_cov <- names(which(is_fac))      # factor covariates
cont_cov <- names(which(!is_fac))       # continuous covariates

L_f <- which(covariates %in% factor_cov)
L_c <- which(covariates %in% cont_cov)
numeric_cov <- rep(F, length(covariates))
numeric_cov[L_c] <- T

if (L > 0){
  # extract layer dependent covariates X: K*m*L
  X <- array(NA, dim = c(K, m, L), dimnames = list(
    rownames(y_mat),      # length K
    colnames(y_mat),      # length m
    covariates  # or "value" if L == 1
  ))
  
  ## if X is already the K x m x L array built in the same order:
  if (length(L_f)>0){
    X_f <- X[, , L_f, drop = FALSE]
  }
  
  if (length(L_c)>0){
    X_c <- X[, , L_c, drop = FALSE]
  }
  
  Xc <- array(NA, dim = c(K, m, L), dimnames = list(
    rownames(y_mat),      # length K
    colnames(y_mat),      # length m
    covariates            # or "value" if L == 1
  ))
  
  if (bij){
    mu_x = matrix(NA, nrow = L, ncol = m)
    S_x = array(NA, dim = c(m, m, L))
  }
  
  for (j in seq_along(covariates)){
    x1 = covariates[j]
    x = make_y_matrix(df, value_var = x1, fill_missing = NA)
    X[,,j] = x
    
    if (numeric_cov[j]){
      if (bij){
        a = .transform_numeric(x, center = center_x, scale = scale_x)
        Xc[,,j] = a$values
        mu_x[j,] = a$center
        S_x[,,j] = diag(a$scale)
      }
    }
    
  }
  if (bij){
    rownames(mu_x) = covariates; colnames(mu_x) = edges
    rownames(S_x) = colnames(S_x) = edges
  } else {
    X1 = matrix(0, nrow = m*K, ncol = L)
    temp2 = expand.grid(edge = 1:m, layer = 1:K) %>% arrange(edge)
    for (i in 1:nrow(temp2)){
      ij = temp2$edge[i]; k= temp2$layer[i]
      X1[i,] = X[k,ij,]
    }
    X1_cont = X1[,numeric_cov,drop=F]
    X1_fact = X1[,!(numeric_cov), drop=F]
    X1c = scale(X1_cont, center = center_x, scale = scale_x)
    mu_x = attr(X1c, "scaled:center")
    sd_x = attr(X1c, "scaled:scale")
    S_x = diag(x = sd_x, nrow = length(sd_x))
    X1_total = cbind(X1c,X1_fact)
    names(mu_x) = rownames(S_x) = colnames(S_x) = covariates
    
    if (length(idx_0) >0){
      rm(X1_total)
      X1 = matrix(0, nrow = m_pos*K, ncol = L)
      temp2 = expand.grid(edge = idx_pos, layer = 1:K) %>% arrange(edge)
      for (i in 1:nrow(temp2)){
        ij = temp2$edge[i]; k= temp2$layer[i]
        X1[i,] = X[k,ij,]
      }
      X1_cont = X1[,numeric_cov,drop=F]
      X1_fact = X1[,!(numeric_cov), drop=F]
      X1c = scale(X1_cont, center = center_x, scale = scale_x)
      mu_x_pos = attr(X1c, "scaled:center")
      sd_x = attr(X1c, "scaled:scale")
      S_x_pos = diag(x = sd_x, nrow = length(sd_x))
      X1_total_pos = cbind(X1c,X1_fact)
      names(mu_x_pos) = rownames(S_x_pos) = colnames(S_x_pos) = covariates
    }
  }
}

# remove temporary variables
rm(Xc,X1_cont, X1_fact, X1, X_c, temp2, x1, X, x, X1c)
gc()

# extract number of layer-INdependent covariates (L'; W)
# extract covariates
covariates_w = all.vars(w)
Lprime = length(covariates_w)

if (Lprime > 0){
  # X' is m x L' matrix
  Xprime_mL = Wc = matrix(NA, nrow = m, ncol = Lprime)
  
  colnames(Xprime_mL) = colnames(Wc) = covariates_w
  
  num_w = rep(F, Lprime)
  Lp_f <- which(covariates_w %in% factor_cov)
  Lp_c <- which(covariates_w %in% cont_cov)
  num_w <- rep(F, length(covariates_w))
  num_w[Lp_c] <- T
  
  for (i in seq_along(covariates_w)){
    c = covariates_w[i]
    x_temp = make_y_matrix(df, value_var = c, fill_missing = NA)
    Xprime_mL[,i] = x_temp[1,]
  }
  cov_w_to_cs = covariates_w[num_w]
  cov_w_fctr = covariates_w[-num_w]
  if (length(cov_w_to_cs) > 0){
    wc1 = scale(Xprime_mL[,cov_w_to_cs], center = center_xprime, scale = scale_xprime)
    mu_w = attr(wc1, "scaled:center")
    s_w = (attr(wc1, "scaled:scale"))
    S_w = diag(x = s_w, nrow = length(s_w))
    names(mu_w) = rownames(S_w) = colnames(S_w) = cov_w_to_cs
  }
  w_star = Xprime_mL[,cov_w_to_cs] / diag(S_w)

  # build W* = [Wf, Wc Sw^(-1)]
  W_star = cbind(Xprime_mL[,cov_w_fctr], w_star)
  colnames(W_star) = c(cov_w_fctr, cov_w_to_cs)
}

# remove temporary variables
rm(w_star, Xprime_mL, x_temp, wc1)
gc()

#Data Inspection ================================================
# countries and their abbreviations
df |> select(reporter_iso3a, reporter_iso2a, reporter_country) |> distinct()
# layers (products) and product codes
df |> select(product_code, product) |> distinct()

# Fit =========================================================
f = 4 # family: ZAGA

## 1. Pi part ####
# compute pi.hat.mle
pi.hat <-  colMeans(y_mat == 0)     # length m

# compute SE for pi
pi_SE <- sqrt((pi.hat * (1-pi.hat))/K)

coef.pi.hat.glm = matrix(NA, nrow = m, ncol = 2)

df = df %>% mutate(value_fac = (value==0), edge = paste0(partner_iso3a, "-", reporter_iso3a))

fit_pi_glm =glm(
  value_fac ~ scale(best) + scale(dist) +
    comlang_ethno + colony + contig +
    partner_iso3a + reporter_iso3a +
    product_code,
  family = binomial,
  data = df
)

df$pi_glm = fit_pi_glm$fitted.values
pi.hat_matrix_temp = make_y_matrix(df, value_var = "pi_glm", fill_missing = NA)
pi.hat_matrix = t(pi.hat_matrix_temp)

## 2. GAMMA #####

W_star_pos = W_star[idx_pos,,drop=F]
y_mat_pos = y_mat[,idx_pos,drop=F]

# generate constraints on W*'a0=0:
constr = make_a0_constraint(W_star_pos)

### Startval ===================================================

if (!("X1_total_pos" %in% ls())){
  X1_total_pos = X1_total
  mu_x_pos = mu_x
  S_x_pos = S_x
} 

parameters1 <- .generate_start_cs(par = NULL, s = 1, family = 4, y = y_mat_pos, X1_total = X1_total_pos, W_star = W_star_pos, bij = bij, xij = xij, mu_x = mu_x_pos, S_x = S_x_pos)

stopifnot(length(parameters1$lambda) == (m_pos-2))
stopifnot(length(constr$rebuild_a0(parameters1$intercept)) == m_pos)


### RTMB setup =================================================
# reconstruct centered and scale Xtilde array
Xt = aperm(array(t(X1_total_pos), dim = c(L, K, m_pos)), c(2, 3, 1))

stopifnot(dim(Xt) == c(K,m_pos,L))

# remove temporary variables and free up memory
rm(X1_total_pos)
gc()

# compute x.tilde offset
cx = drop(mu_x_pos / diag(S_x_pos))

# set up data for model
data_input <- list(y = y_mat_pos, q = q, ws = W_star_pos, x = Xt, cx = cx)

objr <- RTMB::MakeADFun(func = nll(glamle_nll_zaga_bx_ij_mixK_pos_cs, data_input),
                        parameters = parameters1,
                        random = "u",
                        silent = T
)

### Fit ========================================================

fit1 <- try(nlminb(
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

# refit for better convergence

fit2 <- nlminb(
  start = fit1$par,
  objective = objr$fn,
  gradient = objr$gr,
  control = list(
    eval.max = 5000,
    iter.max = 5000,
    rel.tol = 1e-12,
    x.tol = 1e-10
  )
)

fit3 <- optim(
  par = fit2$par,
  fn  = objr$fn,
  gr  = objr$gr,
  method = "BFGS",
  control = list(
    maxit = 5000,
    reltol = 1e-12
  )
)

### Extract params ===========================================
fit1 = fit3
out= list()
grad.sim <- objr$gr(fit1$par)
stopifnot(length(fit1$par) == length(grad.sim))

# Compute the infinity norm of the gradient for diagnostics
(out$Gr$inf.norm = norm(as.matrix(grad.sim), "i"))
# Extract parameter names
param_names <- names(fit1$par)
out$param_names <- param_names
par_index_map <- objr$env$parList(seq_along(fit1$par))
out$par_index_map <- par_index_map

par.hat = objr$env$parList(fit1$par)
if (L>0){
  if (bij){
    rownames(par.hat$beta) = edges
    colnames(par.hat$beta) = covariates
  } else names(par.hat$beta) = covariates
} 
if (Lprime >0) colnames(par.hat$gamma) = c(cov_w_fctr, cov_w_to_cs)

par.hat$a0 = constr$rebuild_a0(par.hat$intercept)
alpha.hat = rbind(par.hat$a0, c(0,1,par.hat$lambda))
colnames(alpha.hat) = edges_pos
rownames(alpha.hat) = c("intercept", paste0("factorL", 1:q))
par.hat$alpha.hat = alpha.hat
par.hat$sig = exp(par.hat$logsig)
mu.hat = .compute_mu_zaga_bx_ij_mixK_pos_cs(parms = par.hat, data_input = data_input)
rownames(mu.hat) = layers; colnames(mu.hat) = edges_pos

# back-transform the parameters to the original scale
a2 = par.hat$alpha.hat[-1,,drop=F]
par.hat.back_transformed = par.hat
beta.hat = par.hat$beta / diag(S_x)
gamma.hat = par.hat$gamma[num_w] / diag(S_w)
par.hat.back_transformed$beta = beta.hat
par.hat.back_transformed$gamma[which(colnames(par.hat$gamma) %in% covariates_w[num_w])] = par.hat$gamma[which(colnames(par.hat$gamma) %in% covariates_w[num_w])] / diag(S_w)

# export and import
m_fit = ncol(mu.hat)
pi.hat_matrix_subset = pi.hat_matrix[rownames(pi.hat_matrix) %in% edges_pos,] %>% t()

DT <- as.data.table(data.frame(edges=edges_pos))

# split each edge into its two countries
edges_fit <- DT[, tstrsplit(edges, "-", fixed = TRUE)]

colnames(edges_fit) = c("exporter", "importer")
## 1. Compute fitted unconditional means:
##    E(Y_ij^k | fitted model) = (1 - pi_ij) * mu_ij^k
EY.hat <- mu.hat * (1-pi.hat_matrix_subset)

df = df %>% mutate(value_fac = (value==0), edge = paste0(partner_iso3a, "-", reporter_iso3a))

W_dist_all = (W_star[,"dist"])*as.numeric(S_w)
W_dist_pos = (W_star_pos[,"dist"])*as.numeric(S_w)
k_id = match(df$product_code, rownames(mu.hat))
e_id = match(df$edge, colnames(mu.hat))
has_fit = !is.na(k_id) & !is.na(e_id)
df$mu.hat[has_fit] = mu.hat[cbind(k_id[has_fit], e_id[has_fit])]
df$ey.hat[has_fit] = EY.hat[cbind(k_id[has_fit], e_id[has_fit])]

e_id = match(df$edge, names(idx_0))
pi1 = !is.na(e_id)
df$ey.hat[pi1] = 0 

e_id = match(df$edge, names(pi.hat))
has_fit = !is.na(e_id)
df$pi.hat[has_fit] = pi.hat[e_id[has_fit]]

### SE =======================================================

gc()

if (se){
  # takes a while to compute
  hess <- numDeriv::jacobian(objr$gr, fit1$par)
  He.inv = ginv(hess)
  # reconstruct D for var-cov matrix (see Supplementary file)
  ## indices in free parameter vector
  idx_a2    <- which(param_names == "lambda")   
  idx_sigma <- which(param_names == "logsig")
  idx_b     <- which(param_names == "beta")
  idx_zeta  <- which(param_names == "intercept")     
  idx_g1    <- which(param_names == "gamma")
  
  ## only numeric gamma entries are rescaled
  idx_g_num <- which(colnames(par.hat$gamma) %in% covariates_w[num_w])
  idx_g     <- idx_g1[idx_g_num]
  
  ## dimensions of reported blocks
  n_a2    <- length(idx_a2)
  n_sigma <- length(idx_sigma)
  n_b     <- length(idx_b)
  n_a0    <- nrow(constr$G)
  n_g     <- length(idx_g1)   
  
  ## build D directly: rows = reported params, cols = free params
  D <- matrix(0, nrow = n_a2 + n_sigma + n_b + n_a0 + n_g,
              ncol = length(param_names))
  
  names_param_rows <- c(rep("a2", n_a2), "logsig", rep("beta",n_b), rep("a0", n_a0), rep("gamma", n_g))
  r <- 1
  
  ## a2_free unchanged
  if (n_a2 > 0) {
    D[r:(r+n_a2-1), idx_a2] <- diag(n_a2)
    r <- r + n_a2
  }
  
  ## sigma unchanged
  if (n_sigma > 0) {
    D[r:(r+n_sigma-1), idx_sigma] <- diag(n_sigma)
    r <- r + n_sigma
  }
  
  ## beta = S_x^{-1} beta_tilde
  if (n_b > 0) {
    D[r:(r+n_b-1), idx_b] <- diag(1 / diag(S_x), n_b)
    r <- r + n_b
  }
  
  ## a0 = G a0,F
  D[r:(r+n_a0-1), idx_zeta] <- constr$G
  r <- r + n_a0
  
  ## gamma block
  
  if (n_g > 0) {
    D[r:(r+n_g-1), idx_g1] <- diag(n_g)
    D[r - 1 + idx_g_num, idx_g] <- diag(1 / diag(S_w), length(idx_g))
  }
  
  if (sandwich.cov){
    # compute var-cov using sandwich method
    # compute Jacobian
    S = .compute_scores_matrix_cs_bxij(y_mat = y_mat_pos, Xt = Xt, W_star_pos = W_star_pos, cx = cx, par.hat = par.hat, p = length(param_names))
    meat  <- crossprod(S)
    V <- He.inv %*% meat %*% He.inv
    V <- 0.5 * (V + t(V))                           # enforce symmetry
    vcov.theta.sandwich <- D %*% V %*% t(D)
    se.theta.sandwich   <- sqrt(diag(vcov.theta.sandwich))
    names(se.theta.sandwich) = names_param_rows
  } else {
    # compute var-cov using inverse Fisher (approximation!)
    vcov.theta.fisher <- D %*% He.inv %*% t(D)
    se.theta.fisher   <- sqrt(diag(vcov.theta.fisher))
    names(se.theta.fisher) = names_param_rows
  }
}


# PPML Fit ====================================================

df_ppml <- df %>% filter(!(is.na(ey.hat)))

fit_ppml_gravity <- ppml(
  dependent_variable = "value",
  distance = "dist",
  additional_regressors = c("comlang_ethno", "colony", "contig","best", "reporter_iso3a", "partner_iso3a"),
  data = df_ppml
)

# Inspect Fit ===================================================

head(df) |> filter(edge == "EEC-CHN") |> select(product, pi_glm, mu.hat, ey.hat)
