get_thetas_for_init_grid_dynWEV_simulations <- function(init_grid, df, nRatings, fixed, mint0, maxtau0) {

  nConds <- length(unique(df$condition))
  nSATs <- length(unique(df$SAT))
  conf_probs <- cumsum(table(df$rating))
  conf_probs <- conf_probs[1:(nRatings-1)]/conf_probs[nRatings]
  
  if (!("t0" %in% names(fixed))) init_grid['t0'] <- init_grid['t0']*mint0
  
  
  get_start_thetas <- function(paramRow) {
    #paramRow['tau'] <- mean(df$RT2)- maxtau0*paramRow['tau0']
    paramRow <- c(paramRow, unlist(fixed, use.names = TRUE))
    
    V <- c(t(paramRow[paste("v", 1:nConds, sep="")]))
    A <- c(t(paramRow[paste("a", 1:nSATs, sep="")]))
    gc()
    res <- expand.grid(condition=1:nConds, SAT=unique(df$SAT)) %>%
      mutate(a=as.numeric(A[SAT]), v=as.numeric(V[condition]))
    paramRow <- as.data.frame(as.list(paramRow))
    conf <- res %>% group_by(SAT, condition) %>%
      reframe(with(paramRow,as.data.frame(RNG_WEV_SAT(1000, .data$a[1], .data$v[1], 
                                        sample(subset(df, SAT==cur_group()$SAT[1] & condition == cur_group()$condition[1])$RT2,
                                               size=1000, replace=TRUE), 
                                        params=c(t0+st0/2, sz, sv,#as.numeric(SV[df[i,]$condition])/s,
                                                 st0, tau0, lambda,
                                                 w, sigvis, svis),
                                        delta = 0.01, maxT =20, stop_on_error=TRUE)))) %>%
      filter(V3 !=0) %>%
      select(V4) 
    thetas <- quantile(conf$V4, probs=conf_probs, names = FALSE)
    c(thetas[1],diff(thetas))
  }
  
  init_thetas <- apply(init_grid, FUN=get_start_thetas, MARGIN=1) # , simplify = TRUE
  init_thetas <- t(init_thetas)
  init_thetas[,2:(nRatings-1)] <- pmax(init_thetas[,2:(nRatings-1)], 0.01)
  init_thetas

}

