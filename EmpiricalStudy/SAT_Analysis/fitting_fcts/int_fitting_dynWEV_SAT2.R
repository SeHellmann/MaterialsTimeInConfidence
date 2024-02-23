fittingdynWEVSAT <- function(df, nConds, nRatings, nSATs, fixed, sym_thetas,
                          optim_method, opts,
                          logging, filename,
                          useparallel, n.cores,precision,
                          used_cats, actual_nRatings){
  ## Be sure that the parallel cluster is stopped if anything happens (error or user interupt)
  on.exit(try(stopCluster(cl), silent = TRUE))
  ## Set restrictions on tau
  mint0 <-  min(df$rt)
  maxtau0 <- min(df$RT2)
  #if (mint0 < 0.2) {mint0 <- 0.2}
  tau0 = c(0.1, 0.3, 0.5)#seq(0.1, 0.8, length.out = 3)
  init_grid <- expand.grid(amin = c(0.5, 1, 2),
                           amax = c(1.5, 2.5, 4, 5),               ### a = distance btw. upper and lower bound \in (0,\infty)]
                           vmin = c(0.01, 0.5, 1),                ### vmin = mean drift rate in first condition \in (0,\infty)]
                           vmax = c(1.4, 2, 3,4),              ### vmax = mean drift rate in last condition \in (\vmin,\infty)]
                           sv = c(0.1, 1, 2),                  ### sv = SD of drift rate (normal distr.) \in (0,\infty)]
                           #z = sum((df$response==1)*df$n)/sum(df$n),### z = mean start point (bias) \in [0,1]
                           sz = c(0.1),                             ### sz = range of possible start points (unif ditr.; in units of a-z) \in [0,1]
                           t0 = c(0.05, 0.2, 0.4),                       ### t0 = proportion of minimal motor time of minimal total response time \in [0,1)
                           #st0 = c(0.1,  0.2),                      ### st0 = range of possible motor times (unif. distr.) \in [0, t0/2]
                           tau0 = tau0,                               ### tau = prop of post-decisional accumulation time (between 0 and 1) of minimal rating response time
                           svis = seq(0.01, 0.5, length.out = 2),   ### svis = variability in visibility accumulation process
                           w = seq(0.2, 0.8, length.out = 3),       ### w = weight bewtween evidence and visibility for confidence judgement
                           sigvis = seq(0.01, 1, length.out = 2),   ### sigivis = between trial variability in drift rate of the visibility process
                           lambda = c(0.3, .8, 1.4))                ### lambda = exponent of accumulation time in the denominator of the confidence variable
  # Remove columns for fixed parameters
  init_grid <- init_grid[setdiff(names(init_grid), names(fixed))]
  init_grid <- unique(init_grid)
  init_grid <- subset(init_grid, (vmax > vmin) & (amax > amin))
  ##### span drifts with quadratic distance
  ###  We assume a different V (mean drift rate) for the different conditions --> nConds parameters
  if (nConds==1) {
    init_grid$v1 <- (init_grid$vmin+init_grid$vmax)/2
  } else {
    for (i in 0:(nConds-1)){
      init_grid[paste("v", i+1, sep="")] <- init_grid$vmin+(i/(nConds-1))*(init_grid$vmax-init_grid$vmin)
    }
  }
  if (nSATs==1) {
    init_grid$a <- (init_grid$amin+init_grid$amax)/2
  } else {
    for (i in 0:(nSATs-1)){
      init_grid[paste("a", i+1, sep="")] <- init_grid$amax - (i/(nSATs-1))*(init_grid$amax-init_grid$amin)
    }
  }
  ## Guess suitable confidence thresholds from theoretical distribution of
  ## the confidence measure and proportion of ratings in the data
  # init_grid <- init_grid[sample(1:nrow(init_grid), 900),]
  Rcpp::sourceCpp("fitting_fcts/RNG_WEV_SAT3.cpp")
  init_thetas <- get_thetas_for_init_grid_dynWEV_simulations(init_grid, df, nRatings, fixed, mint0, maxtau0)


  #### 1.1. For Nelder-Mead transform all parameters to real values ####
  if (sym_thetas) {
    init_grid["theta1"] <- init_thetas[,1]
    if (nRatings > 2) {
      for (i in 2:(nRatings-1)) {
        init_grid[paste("dtheta", i, sep="")] <- init_thetas[,i]
      }
      cols_theta <- c("theta1", paste("dtheta", 2:(nRatings-1), sep=""))
    } else {
      cols_theta <- c("theta1")
    }
  } else {
    init_grid[c("thetaUpper1", "thetaLower1")] <- init_thetas[,1]
    if (nRatings > 2) {
      for (i in 2:(nRatings-1)) {
        init_grid[paste(c("dthetaUpper", "dthetaLower"), i, sep="")] <-  init_thetas[,i]
      }
      cols_theta <- c('thetaLower1', paste("dthetaLower", 2:(nRatings-1), sep=""),
                      'thetaUpper1', paste("dthetaUpper", 2:(nRatings-1), sep=""))
    } else {
      cols_theta <- c("thetaLower1", "thetaUpper1")
    }
  }
  parnames <- c('sz', paste("v", 1:nConds, sep=""),paste("a", 1:nSATs, sep=""),
                 'st0', 'sv', 't0', cols_theta,'tau0', 'w', 'svis', 'sigvis', 'lambda')
  inits <- init_grid[, setdiff(parnames, names(fixed))]
  # remove init_grid
  rm(init_grid)


  ## Intermezzo: Setup cluster for parallelization   ####
  if (useparallel) {
    if (is.null(n.cores)) {
      n.cores <- detectCores()-1
    }
    cl <- makeCluster(type="SOCK", n.cores)
    clusterExport(cl, c("df", "restr_tau", "mint0","maxtau0",
                        "nConds","nRatings", "nSATs", "fixed", "simult_conf", "sym_thetas", "precision"), envir = environment())
    clusterExport(cl, "neglikelihood_2DSDSAT_bounded")
    clusterEvalQ(cl, library(dynConfiR))
    clusterEvalQ(cl, library(tidyverse))
    clusterEvalQ(cl, library(minqa))
    
  }


  ### 2. Search initial grid before optimization  ####
  if (logging==TRUE) {
    logger::log_info(paste(length(inits[,1]), "...parameter sets to check"))
    logger::log_info(paste("data got ", nrow(df), " rows"))
    t00 <- Sys.time()
    logger::log_info("Searching initial values ...")
  }


  if (useparallel) {
    logL <-
      parApply(cl, inits, MARGIN=1,
               function(p) try(neglikelihood_dynWEVSAT_bounded(p, df,  nConds, nRatings, nSATs, fixed, mint0, maxtau0, sym_thetas, precision),
                               silent=TRUE))
    #stopCluster(cl)
  } else {
    logL <-
      apply(inits, MARGIN = 1,
            function(p) try(neglikelihood_dynWEVSAT_bounded(p, df,  nConds, nRatings, nSATs, fixed, mint0, maxtau0,  sym_thetas, precision),
                            silent=TRUE))
  }
  logL <- as.numeric(logL)
  inits <- inits[order(logL),]
  if (logging==TRUE) {
    logger::log_success(paste("Initial grid search took...",as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2))," mins"))
  }
  
  
  
                    # sz,  v1, v2,....,, a1, a2,...,      st0, sv, t0,thetaLower1, dthetaLower2.., thetaUpper1... (or theta1,...),  tau0, w,  svis, sigvis, lambda
  lower_optbound <- c(0,  rep(0, nConds),rep(0.01, nSATs), 0,   0,  0, rep(c(-Inf,  rep(0, nRatings-2)), 2-as.numeric(sym_thetas)),    0, 0,  1e-6,    0,      0)[!(parnames %in% names(fixed))]
  upper_optbound <- c(1, rep(Inf,nConds),rep(Inf, nSATs),Inf, Inf,   1, rep(Inf, (2-as.numeric(sym_thetas))*(nRatings-1)),             1, 1,   Inf,  Inf,     Inf)[!(parnames %in% names(fixed))]




  #### 3. Optimization ####
  if (logging==TRUE) {
    logger::log_info("Start fitting ... ")
  }
  if (!useparallel || (opts$nAttempts==1)) {
    noFitYet <- TRUE
    for (i in 1:opts$nAttempts){
      start <- c(t(inits[i,]))
      names(start) <- names(inits)
      for (l in 1:opts$nRestarts){
        start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
        start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
        try(m <- bobyqa(par = start,
                        fn = neglikelihood_dynWEVSAT_bounded,
                        lower = lower_optbound, upper = upper_optbound,
                        data=df,nConds=nConds, nRatings=nRatings,
                        fixed=fixed, mint0=mint0, maxtau0=maxtau0, nSATs=nSATs, 
                        sym_thetas=sym_thetas, precision=precision,
                        control = list(maxfun=opts$maxfun,
                                       rhobeg = min(0.2, 0.2*max(abs(start))),
                                       npt = length(start)+5)))
        ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
        ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
        ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
        ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
        ## rhoend: use default of 1e-6*rhobeg
        if (exists("m") && !inherits(m, "try-error")){
          m$value <- m$fval
        }
        if (logging==TRUE) {
          logger::log_info(paste("Finished attempt No.", i, " restart no. ", l))
        }
        if (!exists("m") || inherits(m, "try-error")){
          if (logging==TRUE) {
            logger::log_error(paste("No fit obtained at attempt No.", i))
            logger::log_error(paste("Used parameter set", paste(start, sep="", collapse=" "), sep=" ", collapse = ""))
          }
          break
        }
        if (exists("m") && is.list(m)){
          if (noFitYet) {
            fit <- m
            noFitYet <- FALSE
            if (logging==TRUE) {
              logger::log_info(paste("First fit obtained at attempt No.", i))
              attempt <- i
              save(logL, inits,  df,fit, attempt,file=filename)
            }
            start <- fit$par
            names(start) <- names(inits)
          } else if (m$value < fit$value) {
            fit <- m
            if (logging==TRUE) {
              logger::log_info(paste("New fit at attempt No.", i, " restart no. ", l))
              attempt <- i
              save(logL, inits,  df,fit, attempt,file=filename)
            }
            start <- fit$par
            names(start) <- names(inits)
          } # end of if better value
        }   # end of if we got a optim-result at all
      }     # end of for restarts
    }       # end of for initial start values
  } else {  # if useparallel
    starts <- inits[(1:opts$nAttempts),]
    temp_parnames <- names(starts)

    
    optim_node <- function(start) { # define optim-routine to run on each node
      noFitYet <- TRUE
      start <- c(t(start))
      for (l in 1:opts$nRestarts){
        start <- start + rnorm(length(start), sd=pmax(0.001, abs(t(t(start))/20)))
        names(start) <- temp_parnames
        start <- pmax(pmin(start, upper_optbound-1e-6), lower_optbound+1e-6)
        m <- try(bobyqa(par = start,
                        fn = neglikelihood_dynWEVSAT_bounded,
                        lower = lower_optbound, upper = upper_optbound,
                        data=df, nConds=nConds, nRatings=nRatings,
                        fixed=fixed, mint0=mint0,  maxtau0=maxtau0, nSATs=nSATs,
                        sym_thetas=sym_thetas, precision=precision,
                        control = list(maxfun=opts$maxfun,
                                       rhobeg = min(0.2, 0.2*max(abs(start))),
                                       npt = length(start)+5)))
        ## rhobeg should be: about 0.1*(greatest expected change in parameters --> <= 1-2 (for a, thetas or v's) )
        ##                   smaller than min(abs(upper-lower)) = min(1, restr_tau)
        ##                   --> so we use min(0.2*restr_tau, 0.2, 0.2*max(abs(par))).
        ##                   Default would be: min(0.95, 0.2*max(abs(par))), respectively 0.2*max(upper_optbound-lower_optbound)
        ## rhoend: use default of 1e-6*rhobeg
        if (exists("m") && !inherits(m, "try-error")){
          m$value <- m$fval
        }
        if (!exists("m") || inherits(m, "try-error")){
          break
        }
        if (exists("m") && is.list(m)){
          if (noFitYet) {
            fit <- m
            noFitYet <- FALSE
            start <- fit$par
            names(start) <- temp_parnames
          } else if (m$value < fit$value) {
            fit <- m
            start <- fit$par
            names(start) <- temp_parnames
          }
        }
      }
      if (exists("fit") && is.list(fit)){
        return(c(fit$value,fit$par))
      } else {
        return(c(m[1] , rep(NA, length(start))))
      } # end of node-function
    }
    clusterExport(cl, c("temp_parnames", "opts", "optim_method","optim_node" ), envir = environment())
    if (optim_method!="Nelder-Mead") {
      clusterExport(cl, c("lower_optbound", "upper_optbound"), envir = environment())
    }
    optim_outs <- parApply(cl, starts,MARGIN=1, optim_node )
    stopCluster(cl)
    optim_outs <- t(optim_outs)
    best_res <- optim_outs[order(optim_outs[,1]),][1,]
    fit <- list(par = best_res[-1], value=best_res[1])
  }   # end of if-else useparallel

  #### 4. Wrap up results ####
  res <-  data.frame(matrix(nrow=1, ncol=0))
  if(exists("fit") && is.list(fit)){
    k <- length(fit$par)
    N <- sum(df$n)
    p <- fit$par
    p <- c(t(p))
    names(p) <- names(inits)
    if (nRatings>2) {
      if (sym_thetas) {
        p[paste("theta", 2:(nRatings-1), sep="")] <- c(t(p['theta1'])) + c(t(cumsum(p[grep(names(p), pattern = "dtheta", value=TRUE)])))
      } else {
        p[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(p['thetaUpper1'])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaUpper", value=TRUE)])))
        p[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(p['thetaLower1'])) + cumsum(c(t(p[grep(names(p), pattern = "dthetaLower", value=TRUE)])))
      }
      p <- p[ -grep(names(p), pattern="dtheta")]
    }
    if (length(fixed)>=1) p <- c(p, unlist(fixed))
    if (!("t0" %in% names(fixed))) p['t0'] <- p['t0']*mint0
    if (!("tau0" %in% names(fixed))) p['tau0'] <- p['tau0']*maxtau0
    # if (!("sz" %in% names(fixed))) p['sz'] <- (min(p['z'], (1-p['z']))*2)*p["sz"] ##

    res <-   data.frame(matrix(nrow=1, ncol=length(p)))
    res[1,] <- p
    names(res) <- names(p)      # a,  z, sz,v1, v2,....,,   st0, sv, t0, thetaLower1,dthetaLower2-4,   thetaUpper1,dthetaUpper2-4,    tau,       w, svis, sigvis

    if (!is.null(used_cats)) {
      # If some rating categories are not used, we fit less thresholds numerically and fill up the
      # rest by the obvious best-fitting thresholds (e.g. +/- Inf for the lowest/highest...)
      res <- fill_thresholds(res, used_cats, actual_nRatings, -1e+24)
      nRatings <- actual_nRatings
      k <- ncol(res)
    }
    if (sym_thetas) {
      parnames <- c(paste("v", 1:nConds, sep=""), paste("a", 1:nSATs, sep=""), 'sv', 'z', 'sz', 't0','st0', paste("theta", 1:(nRatings-1), sep=""), 'tau0', 'w', 'svis','sigvis', 'lambda')
    } else {
      parnames <- c(paste("v", 1:nConds, sep=""), paste("a", 1:nSATs, sep=""), 'sv', 'z', 'sz', 't0','st0', paste("thetaLower", 1:(nRatings-1), sep=""), paste("thetaUpper", 1:(nRatings-1), sep=""), 'tau0', 'w', 'svis','sigvis', 'lambda')
    }
    res <- res[, parnames]

    res$fixed <- paste(c("sym_thetas", names(fixed)), c(sym_thetas,fixed), sep="=", collapse = ", ")
    res$negLogLik <- fit$value
    res$N <- N
    res$k <- k
    res$BIC <-  2 * fit$value + k * log(N)
    res$AICc <- 2 * fit$value + k * 2 + 2*k*(k-1)/(N-k-1)
    res$AIC <- 2 * fit$value + k * 2
    if (logging==TRUE) {
      logger::log_success("Done fitting and autosaved results")
      save(logL, df, res, file=filename)
    }
  }
  return(res)
}


neglikelihood_dynWEVSAT_bounded <-   function(p, data,
                                           nConds, nRatings, nSATs, fixed, mint0, maxtau0, sym_thetas=FALSE, precision=1e-5)
{
  # get parameter vector back from real transformations
  paramDf <-   data.frame(matrix(nrow=1, ncol=length(p)))
  paramDf[1,] <- p
  names(paramDf) <- names(p)
  if (nRatings > 2) {
    if (sym_thetas) {
      paramDf[paste("theta", 2:(nRatings-1), sep="")] <- c(t(paramDf['theta1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dtheta", value=TRUE)])))
    } else {
      paramDf[paste("thetaUpper", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaUpper1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaUpper", value=TRUE)])))
      paramDf[paste("thetaLower", 2:(nRatings-1), sep="")] <- c(t(paramDf['thetaLower1'])) + cumsum(c(t(paramDf[grep(names(paramDf), pattern = "dthetaLower", value=TRUE)])))
    }
    paramDf <- paramDf[ -grep(names(paramDf), pattern="dtheta")]
  }
  if (length(fixed)>=1) {
    paramDf <- cbind(paramDf, as.data.frame(fixed))
  }
  if (!("t0" %in% names(fixed))) paramDf['t0'] <- paramDf['t0']*mint0
  if (!("tau0" %in% names(fixed))) paramDf['tau0'] <- paramDf['tau0']*maxtau0
  #if (!("sz" %in% names(fixed))) paramDf['sz']  <-  (min(p['z'], (1-p['z']))*2)*p["sz"] ##
  
  negloglik <- -LogLikWEVSAT(data, paramDf, "dynaViTE", precision, stop_on_error=FALSE)
  return(negloglik)
}
