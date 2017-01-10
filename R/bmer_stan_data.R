## these are functions for internal use in build_bmer_model

# assemble list of data to pass to stan::sampling
# and replacement names for parameters
# checked is returned by check_variables
# control is a bmerControl
bmer_stan_data <- function(checked,control){
  mess <- warn <- character()
  
	modmat <- cbind(checked$x,do.call(cbind,args=checked$z))
  random <- checked$varlog$random

	pars <- stan_data <- list()
	stan_data$N <- nrow(modmat)
  stan_data$K <- ncol(modmat)
	stan_data$P <- ncol(checked$x)

  pars$beta <- list(name = "beta", dims = list(Effect = colnames(checked$x)))
  pq <- stan_data$P
  if( (control@estimate_scale_beta=="P" & stan_data$P >= 10)  |  control@estimate_scale_beta=="yes" ){
    pars$scale_beta <- list(name = "scale_beta", dims = list("scale_beta"))
    pq <- pq + 1
  }

  stan_data$G <- length(checked$z)
  coef_indx <- matrix(c(1,stan_data$P),1,2)
  qslopes <- FALSE
  for(g in 1:length(checked$z)){
    stan_data[[paste("M",g,sep="_")]] <- random[[g]]$nlevels
    stan_data[[paste("Q",g,sep="_")]] <- random[[g]]$Q
    coef_indx <- rbind(coef_indx,c(coef_indx[g,2]+1,coef_indx[g,2]+ncol(checked$z[[g]])))
    pars[[paste("gamma",names(random)[g],sep="_")]] <- list(name = paste("gamma",g,sep="_"),
      dims = list(Member = random[[g]]$levels, Effect = random[[g]]$cols))
    pars[[paste("sigma",names(random)[g],sep="_")]] <- list(name = paste("sigma",g,sep="_"),
      dims = list(Effect = random[[g]]$cols))
    pq <- pq + (1 + random[[g]]$nlevels) * random[[g]]$Q
    if(random[[g]]$Q>1){
      qslopes <- TRUE
      pars[[paste("omega",names(random)[g],sep="_")]] <- list(name = paste("omega",g,sep="_"),
        dims = list(Effect1 = random[[g]]$cols, Effect2 = random[[g]]$cols),
        keep = as.vector(lower.tri(diag(random[[g]]$Q))))
      pq <- pq + choose(random[[g]]$Q,2)
    }
  }
  coef_indx <- coef_indx[-1,]
  if(!is.matrix(coef_indx)) coef_indx <- array(coef_indx,dim=c(1,2))
  stan_data$cindx <- coef_indx

  modmat <- extract_sparse_parts(modmat)
	stan_data$nz <- length(modmat$w)
	stan_data$w <- modmat$w
	stan_data$v <- modmat$v
	stan_data$u <- modmat$u

  stan_data$y <- checked$model_frame[,1]
  
  if( (control@estimate_scale_beta=="P" & stan_data$P >= 10)  |  control@estimate_scale_beta=="yes" ){
    stan_data$sc_beta <- control@sc_beta
  } else {
    stan_data$scale_beta <- control@scale_beta
  }
  if(control@nu_beta!=1) stan_data$nu_beta <- control@nu_beta
  stan_data$sc_q0 <- control@sc_q0
  if(qslopes){
    stan_data$sc_qs <- control@sc_qs
    stan_data$eta_q <- control@eta_q
  }
  if(checked$family=="gaussian"){
    stan_data$sc_res <- control@sc_res
    pars$sigma_residual <- list(name = "sigma_res", dims = list("sigma_residual"))
    pq <- pq + 1
  }

  pars$y_hat <- list(name = "y_hat", dims = list(Obs = paste(1:stan_data$N)))
  pars$log_lik <- list(name = "log_lik", dims = list(Obs = paste(1:stan_data$N)))

  attr(pars,"pq") <- pq

	getRowNames <- function(dims){
		g <- list()
		for(d in length(dims):1){
			if(!is.null(dims[[d]])){
				gd <- dims[[d]]
				if(is.data.frame(gd)){
					for(j in 2:ncol(gd)){
						gd[,j] <- paste(colnames(gd)[j],gd[,j],sep="")
					}
					g[[length(g)+1]] <- apply(gd,1,function(x){paste(x[-1],collapse=",")})
				} else if(is.list(gd)){
					g[[length(g)+1]] <- getRowNames(gd)
				} else {
					g[[length(g)+1]] <- gd
				}
			}
		}
		g <- expand.grid(g)
		if(ncol(g)>1) g <- g[,ncol(g):1]
		return(apply(g,1,function(x){paste(x,collapse=",")}))
	}
	
	for(p in 1:length(pars)){
		rows <- getRowNames(pars[[p]]$dims)
		if(length(rows)==1){
			pars[[p]]$rownames <- rows
		} else {
			pars[[p]]$rownames <- paste(names(pars)[p],"[",rows,"]",sep="")
		}
    if(!("keep" %in% names(pars[[p]]))) pars[[p]]$keep <- rep(TRUE,length(pars[[p]]$rownames))
	}
	
	if(pq >= stan_data$N) warn <- c(warn,"WARNING: more unconstrained parameters than observations")
	
	return(list(stan_data = stan_data, pars = pars, messages = mess, warnings = warn))
}
