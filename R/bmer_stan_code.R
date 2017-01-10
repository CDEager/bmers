## these are functions for internal use in build_bmer_model

# get raw uncommented code
# stanmod is a character string separated with underscores:
#    family name (gaussian/binomial)
#    estimate scale beta (T/F)
#    beta cauchy (T/F)
#    number of random effect groups with slopes
#    number of random effect groups with intercepts only
raw_stan_code <- function(stanmod){
  nm <- stanmod
  stanmod <- str_split(stanmod,"_")[[1]]
  family <- stanmod[1]
  estimate_scale_beta <- as.logical(stanmod[2])
  bc <- as.logical(stanmod[3])
  Qs <- as.numeric(stanmod[4])
  Q0 <- as.numeric(stanmod[5])
  
  gau <- family == "gaussian"

  code <- paste("// bmers (",bmers_version(),") model ",nm,sep="")

  if(Qs>0){
    code <- c(code,
      "",
      "functions {",
      "  matrix vec_to_mat_by_row(int R, int C, vector v){",
      "    matrix[R,C] m;",
      "    for(r in 1:R) m[r] = v[(C*(r-1)+1):(C*r)]';",
      "    return m;",
      "  }",
      "}")
  }
  
  code <- c(code,
    "",
    "data {",
    "  int<lower=0> N;  // number of observations",
    "  int<lower=0> K;  // number of coefficients",
    "",
    "  int<lower=0> nz;  // num non-zero elements in model matrix",
    "  vector[nz] w;  // non-zero elements in model matrix",
    "  int<lower=0> v[nz];  // column indices for w",
    "  int<lower=0> u[N+1];  // row-start indices for non-zero elements in model matrix",
    "")

  if(gau){
    code <- c(code,"  vector[N] y;  // scaled response","")
  } else {
    code <- c(code,"  int<lower=0,upper=1> y[N];  // binary response as integer 0/1","")
  }

  code <- c(code,
    "  int<lower=0> P;  // number of fixed effects",
    "  int<lower=0> G;  // number of random effect groups",
    "  int<lower=0> cindx[G,2];  // coefficient index for random effects")

  for(r in 1:sum(Qs,Q0)){
    code <- c(code,str_replace_all(c(
      "  int<lower=0> M_#;  // number of group # members",
      "  int<lower=0> Q_#;  // number of group # effects per member"),"#",r))
  }

  code <- c(code,"","  // (hyper) priors")

  if(estimate_scale_beta){
    code <- c(code,"  real<lower=0> sc_beta;  // prior scale on scale for beta prior")
  } else {
    code <- c(code,"  real<lower=0> scale_beta;  // prior scale for betas")
  }

  code <- c(code,
    "  real<lower=0> nu_beta;  // degrees of freedom for beta t-dist prior",
    "  real<lower=0> sc_q0;  // prior scale for random intercept standard deviations")
  if(bc) code <- code[-(length(code)-1)]

  if(Qs>0){
    code <- c(code,
      "  real<lower=0> sc_qs;  // prior scale for random slope standard deviations",
      "  real<lower=0> eta_q;  // shape for LKJ prior on random effects correlations")
  }
  
  if(gau) code <- c(code,"  real<lower=0> sc_res;  // prior scale for standard deviation of the residuals")

  code <- c(code,
    "}",
    "",
    "parameters {",
    "  // all parameters sampled on unit scale or with cholesky factors and reparameterized",
    "",
    "  vector[P] beta_raw;")
  if(bc) code[length(code)] <- "  vector<lower=-pi()/2,upper=pi()/2>[P] beta_raw;"

  if(estimate_scale_beta) code <- c(code,"  real<lower=0> scale_beta_raw;")

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "",
        "  matrix[Q_#,M_#] gamma_#_raw;",
        "  vector<lower=0>[Q_#] sigma_#_raw;",
        "  cholesky_factor_corr[Q_#] omega_#_raw;"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "",
        "  vector[M_#] gamma_#_raw;",
        "  real<lower=0> sigma_#_raw;"),"#",p))
      p <- p + 1
    }
  }

  if(gau) code <- c(code,"","  real<lower=0> sigma_res_raw;")

  code <- c(code,
    "}",
    "",
    "transformed parameters {")

  if(estimate_scale_beta) code <- c(code,"  real<lower=0> scale_beta;  // scale of betas' t-distribution")

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "  vector<lower=0>[Q_#] sigma_#;  // standard deviation in the group # effects"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "  real<lower=0> sigma_#;  // standard deviation in the group # intercepts"),"#",p))
      p <- p + 1
    }
  }

  if(gau) code <- c(code,"  real<lower=0> sigma_res;  // standard deviation of the residuals")

  code <- c(code,
    "",
    "  vector[K] coef;  // all coefficients",
    "  vector[N] y_hat;  // fitted values",
    "")

  if(estimate_scale_beta) code <- c(code,"  scale_beta = sc_beta * scale_beta_raw;")
  if(bc){
    code <- c(code,"  for(p in 1:P) coef[p] = scale_beta * tan(beta_raw[p]);")
  } else {
    code <- c(code,"  coef[1:P] = scale_beta * beta_raw;")
  }

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "",
        "  sigma_#[1] = sc_q0 * sigma_#_raw[1];",
        "  sigma_#[2:Q_#] = sc_qs * sigma_#_raw[2:Q_#];",
        "  coef[cindx[#,1]:cindx[#,2]] = to_vector(rep_matrix(sigma_#,M_#) .* (omega_#_raw * gamma_#_raw));"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "",
        "  sigma_# = sc_q0 * sigma_#_raw;",
        "  coef[cindx[#,1]:cindx[#,2]] = sigma_# * gamma_#_raw;"),"#",p))
      p <- p + 1
    }
  }

  if(gau) code <- c(code,"","  sigma_res = sc_res * sigma_res_raw;")

  code <- c(code,
    "",
    "  y_hat = csr_matrix_times_vector(N,K,w,v,u,coef);",
    "}",
    "",
    "model {",
    "  beta_raw ~ student_t(nu_beta,0,1);")
  if(bc) code[length(code)] <- "  beta_raw ~ uniform(-pi()/2,pi()/2);"

  if(estimate_scale_beta) code <- c(code,"  scale_beta_raw ~ normal(0,1);")

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "",
        "  to_vector(gamma_#_raw) ~ normal(0,1);",
        "  sigma_#_raw ~ normal(0,1);",
        "  omega_#_raw ~ lkj_corr_cholesky(eta_q);"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "",
        "  gamma_#_raw ~ normal(0,1);",
        "  sigma_#_raw ~ normal(0,1);"),"#",p))
      p <- p + 1
    }
  }

  if(gau){
    code <- c(code,
      "",
      "  sigma_res_raw ~ normal(0,1);",
      "  y ~ normal(y_hat,sigma_res);")
  } else {
    code <- c(code,"","  y ~ bernoulli_logit(y_hat);")
  }

  code <- c(code,
    "}",
    "",
    "generated quantities {",
    "  vector[N] log_lik;  // log-likelihod",
    "  vector[P] beta;  // fixed effects")

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "  matrix[M_#,Q_#] gamma_#;  // group # effects",
        "  matrix[Q_#,Q_#] omega_#;  // correlation in the group # effects"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "  matrix[M_#,Q_#] gamma_#;  // group # intercepts"),"#",p))
      p <- p + 1
    }
  }

  if(gau){
    code <- c(code,"","  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | y_hat[n],sigma_res);")
  } else {
    code <- c(code,"","  for(n in 1:N) log_lik[n] = bernoulli_logit_lpmf(y[n] | y_hat[n]);")
  }
  
  code <- c(code,"  beta = coef[1:P];")

  p <- 1
  if(Qs>0){
    for(r in 1:Qs){
      code <- c(code,str_replace_all(c(
        "  gamma_# = vec_to_mat_by_row(M_#,Q_#,coef[cindx[#,1]:cindx[#,2]]);",
        "  omega_# = tcrossprod(omega_#_raw);"),"#",p))
      p <- p + 1
    }
  }
  if(Q0>0){
    for(r in 1:Q0){
      code <- c(code,str_replace_all(c(
        "  gamma_#[1:M_#,1] = coef[cindx[#,1]:cindx[#,2]];"),"#",p))
      p <- p + 1
    }
  }

  code <- c(code,"}","")

  code <- paste(code,collapse="\n")

  return(code)
}


# get commented code; it is a bit redundant to run both commented_stan_code and raw_stan_code but they also
# both take essentially no computational time
commented_stan_code <- function(family, random, P, model_name, control){
  for(r in 1:length(random)){
    random[[r]]$s <- random[[r]]$Q > 1
    random[[r]]$num <- r
    random[[r]]$name <- names(random)[r]
  }

  gau <- family == "gaussian"

  estimate_scale_beta <- (control@estimate_scale_beta=="P" & P>=10)  |  control@estimate_scale_beta=="yes"
  
  bc <- control@nu_beta == 1

  code <- c(paste("// Stan mixed effects regression built with bmers version ",bmers_version()," for model '",model_name,"'",sep=""),
    "",
    "// Priors",
    "  // Fixed Effects",
    "    // beta ~ student_t(nu_beta,0,scale_beta)")
  if(bc) code[length(code)] <- "    // beta ~ cauchy(0,scale_beta)"

  if(estimate_scale_beta) code <- c(code,"    // scale_beta ~ half_normal(0,sc_beta)")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
      "",
      "  // Random Effects for grp (Group #)",
      "    // gamma_# ~ mult_normal(0,Sigma_#)",
      "    // Sigma_# = diag(sigma_#) * omega_# * diag(sigma_#)'",
      "    // sigma_#[1] ~ half_normal(0,sc_q0)",
      "    // sigma_#[2:Q_#] ~ half_normal(0,sc_qs)",
      "    // omega_# ~ lkj_corr(eta_q)"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
      "",
      "  // Random Effects for grp (Group #)",
      "    // gamma_# ~ normal(0,sigma_#)",
      "    // sigma_# ~ half_normal(0,sc_q0)"),"#",r$num),"grp",r$name))
    }
  }

  if(gau){
    code <- c(code,
      "",
      "  // Response Prior (Gaussian)",
      "    // sigma_res ~ half_normal(0,sc_res)",
      "    // y ~ normal(Xb + Zg, sigma_res)")
  } else {
    code <- c(code,
      "",
      "  // Response Prior (Logistic)",
      "    // y ~ bernoulli_logit(Xb + Zg)")
  }

  if(any(unlist(lapply(random,function(x) x$s)))){
    code <- c(code,
      "",
      "functions {",
      "  matrix vec_to_mat_by_row(int R, int C, vector v){",
      "    matrix[R,C] m;",
      "    for(r in 1:R) m[r] = v[(C*(r-1)+1):(C*r)]';",
      "    return m;",
      "  }",
      "}")
  }
  
  code <- c(code,
    "",
    "data {",
    "  int<lower=0> N;  // number of observations",
    "  int<lower=0> K;  // number of coefficients",
    "",
    "  int<lower=0> nz;  // num non-zero elements in model matrix",
    "  vector[nz] w;  // non-zero elements in model matrix",
    "  int<lower=0> v[nz];  // column indices for w",
    "  int<lower=0> u[N+1];  // row-start indices for non-zero elements in model matrix",
    "")

  if(gau){
    code <- c(code,"  vector[N] y;  // scaled response","")
  } else {
    code <- c(code,"  int<lower=0,upper=1> y[N];  // binary response as integer 0/1","")
  }

  code <- c(code,
    "  int<lower=0> P;  // number of fixed effects",
    "  int<lower=0> G;  // number of random effect groups",
    "  int<lower=0> cindx[G,2];  // coefficient index for random effects")

  for(r in random){
    code <- c(code,str_replace_all(str_replace_all(c(
      "  int<lower=0> M_#;  // number of grp members",
      "  int<lower=0> Q_#;  // number of grp effects per member"),"#",r$num),"grp",r$name))
  }

  code <- c(code,"","  // (hyper) priors")

  if(estimate_scale_beta){
    code <- c(code,"  real<lower=0> sc_beta;  // prior scale on scale for beta prior")
  } else {
    code <- c(code,"  real<lower=0> scale_beta;  // prior scale for betas")
  }

  code <- c(code,
    "  real<lower=0> nu_beta;  // degrees of freedom for beta t-dist prior",
    "  real<lower=0> sc_q0;  // prior scale for random intercept standard deviations")
  if(bc) code <- code[-(length(code)-1)]

  if(any(unlist(lapply(random,function(x) x$s)))){
    code <- c(code,
      "  real<lower=0> sc_qs;  // prior scale for random slope standard deviations",
      "  real<lower=0> eta_q;  // shape for LKJ prior on random effects correlations")
  }
  
  if(gau) code <- c(code,"  real<lower=0> sc_res;  // prior scale for standard deviation of the residuals")

  code <- c(code,
    "}",
    "",
    "parameters {",
    "  // all parameters sampled on unit scale or with cholesky factors and reparameterized",
    "",
    "  vector[P] beta_raw;")
  if(bc) code[length(code)] <- "  vector<lower=-pi()/2,upper=pi()/2>[P] beta_raw;"

  if(estimate_scale_beta) code <- c(code,"  real<lower=0> scale_beta_raw;")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  matrix[Q_#,M_#] gamma_#_raw;",
        "  vector<lower=0>[Q_#] sigma_#_raw;",
        "  cholesky_factor_corr[Q_#] omega_#_raw;"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  vector[M_#] gamma_#_raw;",
        "  real<lower=0> sigma_#_raw;"),"#",r$num),"grp",r$name))
    }
  }

  if(gau) code <- c(code,"","  real<lower=0> sigma_res_raw;")

  code <- c(code,
    "}",
    "",
    "transformed parameters {")

  if(estimate_scale_beta) code <- c(code,"  real<lower=0> scale_beta;  // scale of betas' t-distribution")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "  vector<lower=0>[Q_#] sigma_#;  // standard deviation in the grp effects"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "  real<lower=0> sigma_#;  // standard deviation in the grp intercepts"),"#",r$num),"grp",r$name))
    }
  }

  if(gau) code <- c(code,"  real<lower=0> sigma_res;  // standard deviation of the residuals")

  code <- c(code,
    "",
    "  vector[K] coef;  // all coefficients",
    "  vector[N] y_hat;  // fitted values",
    "")

  if(estimate_scale_beta) code <- c(code,"  scale_beta = sc_beta * scale_beta_raw;")
  if(bc){
    code <- c(code,"  for(p in 1:P) coef[p] = scale_beta * tan(beta_raw[p]);")
  } else {
    code <- c(code,"  coef[1:P] = scale_beta * beta_raw;")
  }

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  sigma_#[1] = sc_q0 * sigma_#_raw[1];",
        "  sigma_#[2:Q_#] = sc_qs * sigma_#_raw[2:Q_#];",
        "  coef[cindx[#,1]:cindx[#,2]] = to_vector(rep_matrix(sigma_#,M_#) .* (omega_#_raw * gamma_#_raw));"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  sigma_# = sc_q0 * sigma_#_raw;",
        "  coef[cindx[#,1]:cindx[#,2]] = sigma_# * gamma_#_raw;"),"#",r$num),"grp",r$name))
    }
  }

  if(gau) code <- c(code,"","  sigma_res = sc_res * sigma_res_raw;")

  code <- c(code,
    "",
    "  y_hat = csr_matrix_times_vector(N,K,w,v,u,coef);",
    "}",
    "",
    "model {",
    "  beta_raw ~ student_t(nu_beta,0,1);")
  if(bc) code[length(code)] <- "  beta_raw ~ uniform(-pi()/2,pi()/2);"

  if(estimate_scale_beta) code <- c(code,"  scale_beta_raw ~ normal(0,1);")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  to_vector(gamma_#_raw) ~ normal(0,1);",
        "  sigma_#_raw ~ normal(0,1);",
        "  omega_#_raw ~ lkj_corr_cholesky(eta_q);"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "",
        "  gamma_#_raw ~ normal(0,1);",
        "  sigma_#_raw ~ normal(0,1);"),"#",r$num),"grp",r$name))
    }
  }

  if(gau){
    code <- c(code,
      "",
      "  sigma_res_raw ~ normal(0,1);",
      "  y ~ normal(y_hat,sigma_res);")
  } else {
    code <- c(code,"","  y ~ bernoulli_logit(y_hat);")
  }

  code <- c(code,
    "}",
    "",
    "generated quantities {",
    "  vector[N] log_lik;  // log-likelihod",
    "  vector[P] beta;  // fixed effects")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "  matrix[M_#,Q_#] gamma_#;  // grp effects",
        "  matrix[Q_#,Q_#] omega_#;  // correlation in the grp effects"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "  matrix[M_#,Q_#] gamma_#;  // grp intercepts"),"#",r$num),"grp",r$name))
    }
  }

  if(gau){
    code <- c(code,"","  for(n in 1:N) log_lik[n] = normal_lpdf(y[n] | y_hat[n],sigma_res);")
  } else {
    code <- c(code,"","  for(n in 1:N) log_lik[n] = bernoulli_logit_lpmf(y[n] | y_hat[n]);")
  }
  
  code <- c(code,"  beta = coef[1:P];")

  for(r in random){
    if(r$s){
      code <- c(code,str_replace_all(str_replace_all(c(
        "  gamma_# = vec_to_mat_by_row(M_#,Q_#,coef[cindx[#,1]:cindx[#,2]]);",
        "  omega_# = tcrossprod(omega_#_raw);"),"#",r$num),"grp",r$name))
    } else {
      code <- c(code,str_replace_all(str_replace_all(c(
        "  gamma_#[1:M_#,1] = coef[cindx[#,1]:cindx[#,2]];"),"#",r$num),"grp",r$name))
    }
  }

  code <- c(code,"}","")

  code <- paste(code,collapse="\n")

  return(code)
}
