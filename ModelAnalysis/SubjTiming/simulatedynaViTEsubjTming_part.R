simulatedynaViTE_subjtimig <- function (paramDf, n=1e+4, simult_conf = TRUE,
                         stimulus = c(-1,1), delta=0.01, maxrt=15,
                         process_results = FALSE)
{

  lambda <- paramDf$lambda

  ## recover parameters from paramDf
  a <- paramDf$a
  z <- paramDf$z
  sz <- paramDf$sz
  t0 <- paramDf$t0  # recalc_t0 (see e.g. ddynaViTE)
  st0 <- paramDf$st0
  tau = paramDf$tau
  Timerate <- paramDf$Timerate
  Timediff <- paramDf$Timediff
  d <- 0

  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1
  if (vary_sv){
    SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  } else {
    SV <- rep(paramDf$sv, nConds)
  }
  vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  if (vary_s){
    S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  } else {
    if ("s" %in% names(paramDf)) {
      S <- rep(paramDf$s, nConds)
    } else {
      S <- rep(1, nConds)
    }
  }
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }


  w = paramDf$w
  svis = paramDf$svis
  sigvis = paramDf$sigvis
  if ("muvis" %in% names(paramDf)) {
    muvis <- rep(paramDf$muvis, nConds)
  } else {
    muvis <- abs(V)
  }

  df <- expand.grid(condition = 1:nConds, stimulus=stimulus)
  ## Produce process outcomes and compute confidence measure
  simus <- data.frame()
  for ( i in 1:nrow(df)) {
    s = as.numeric(S[df[i,]$condition])
    temp <- as.data.frame(RNG_WEV_SubjTime(n=n,
                                params=c(a/s,
                                         as.numeric(V[df[i,]$condition]*df[i,]$stimulus)/s,
                                         t0+st0/2, d, sz, as.numeric(SV[df[i,]$condition])/s,
                                         st0, z, tau,
                                         lambda,
                                         w, # if model=2DSD, this is set to 1
                                         # if model=2DSD, these are set arbitrary
                                         as.numeric(muvis[df[i,]$condition])/s,
                                         sigvis/s, svis/s, Timerate, Timediff),
                                delta = delta, maxT =maxrt, stop_on_error=TRUE))
    names(temp) <- c("rt", "response", "conf", "dec", "vis", "mu")
    temp$conf <- temp$conf * s
    temp$dec  <- temp$dec * s
    temp$vis  <- temp$vis * s

    simus <- rbind(simus,
                   cbind(condition=df[i, "condition"], stimulus=df[i, "stimulus"], temp))
  }
  if (!process_results) simus <- select(simus, -c("dec", "vis", "mu"))

  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
  } else {
    thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
  }

  levels_lower <- cumsum(as.numeric(table(thetas_lower)))
  levels_lower <- levels_lower[-length(levels_lower)]
  levels_upper <- cumsum(as.numeric(table(thetas_upper)))
  levels_upper <- levels_upper[-length(levels_upper)]
  thetas_lower <- unique(thetas_lower)
  thetas_upper <- unique(thetas_upper)

  simus$rating <- 1
  simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                                 breaks=thetas_upper, labels = levels_upper)))
  simus$rating[simus$response==-1] <- as.numeric(as.character(cut(simus$conf[simus$response==-1],
                                                                  breaks=thetas_lower, labels = levels_lower)))

  simus$correct <- as.numeric(simus$response == simus$stimulus)
  if (simult_conf) {
    simus$rt <- simus$rt + tau
  }

  return(simus)
}







simulatedynaViTE_SAT_subjtimig <- function (paramDf, n=1e+4,  df, model = "dynWEV",
                             delta=0.01, maxrt=15, seed=NULL)
{
  lambda <- paramDf$lambda
  ## get parameters from paramDf
  #z <- paramDf$z
  sz <- paramDf$sz
  t0 <- paramDf$t0  # recalc_t0 (see e.g. dWEV)
  st0 <- paramDf$st0
  tau0 = paramDf$tau0
  d <- 0
  Timerate <- paramDf$Timerate
  Timediff <- paramDf$Timediff
  
  nConds <- length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
  if (nConds > 0 ) {
    V <- c(t(paramDf[,paste("v",1:(nConds), sep = "")]))
  } else {
    V <- paramDf$v
    nConds <- 1
  }
  nSATConds <- length(grep(pattern = "^a[0-9]", names(paramDf), value = T))
  if (nSATConds > 0 ) {
    A <- c(t(paramDf[,paste("a",1:(nSATConds), sep = "")]))
  } else {
    A <- paramDf$a
    nSATConds <- 1
  }
  sv <- paramDf$sv
  symmetric_confidence_thresholds <- length(grep(pattern = "thetaUpper", names(paramDf), value = T))<1
  if (symmetric_confidence_thresholds) {
    nRatings <- length(grep(pattern = "^theta[0-9]", names(paramDf)))+1
  } else {
    nRatings <- length(grep(pattern = "^thetaUpper[0-9]", names(paramDf)))+1
  }
  
  if (grepl("2DSD",paramDf$model)) {
    w = 1
    svis = 1
    sigvis = 0
    muvis = rep(0, nConds)
  } else {
    w = paramDf$w
    svis = paramDf$svis
    sigvis = paramDf$sigvis
    if ("muvis" %in% names(paramDf)) {
      muvis <- rep(paramDf$muvis, nConds)
    } else {
      muvis <- abs(V)
    }
  }
  
  #RT2s_sim <- sample(RT2s, size=n, replace=TRUE)
  gc()
  res <- expand.grid(condition=1:nConds, SAT=1:nSATConds) %>%
    mutate(a=A[SAT], v=V[condition],
           SAT=ifelse(SAT==1, "Accuracy", "Speed"))
  
  # a= res[1, ]$a
  # v = res[1, ]$v
  # RT2s <- sample(subset(df, SATnum==1 & condition == 1)$RT2,
  #                size=n, replace=TRUE)
  # RNG_WEV_SAT(n, a, v, RT2s, params=c(t0+st0/2, sz, sv,#as.numeric(SV[df[i,]$condition])/s,
  #                                     st0, tau0, lambda,
  #                                     w, # if model=2DSD, this is set to 1
  #                                     # if model=2DSD, these are set arbitrary
  #                                     # muvis: as.numeric(muvis[df[i,]$condition]),
  #                                     sigvis, svis),
  #             delta = 0.01, maxT =10, stop_on_error=TRUE)
  params <- c(t0+st0/2, sz, sv,#as.numeric(SV[df[i,]$condition])/s,
              st0, tau0, lambda,
              w, # if model=2DSD, this is set to 1
              # if model=2DSD, these are set arbitrary
              # muvis: as.numeric(muvis[df[i,]$condition]),
              sigvis, svis,
              Timerate, Timediff)
  #print(params)
  simus <- res %>% group_by(SAT, condition) %>%
    reframe(as.data.frame(RNG_WEV_SAT_SubjTime(n, .data$a[1], .data$v[1], 
                                      sample(subset(df, SAT==cur_group()$SAT[1] & condition == cur_group()$condition[1])$RT2,
                                             size=n, replace=TRUE), 
                                      params=params,
                                      delta = delta, maxT =maxrt, stop_on_error=TRUE)))
  #simus <- as.data.frame(simus)
  names(simus) <- c("SAT", "condition", "RT2", "rt", "response", "conf")
  if (symmetric_confidence_thresholds) {
    thetas_upper <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("theta",1:(nRatings-1), sep = "")]), Inf)
  } else {
    thetas_upper <- c(-Inf, t(paramDf[,paste("thetaUpper",1:(nRatings-1), sep = "")]), Inf)
    thetas_lower <- c(-Inf, t(paramDf[,paste("thetaLower",1:(nRatings-1), sep="")]), Inf)
  }
  
  levels_lower <- cumsum(as.numeric(table(thetas_lower)))
  levels_lower <- levels_lower[-length(levels_lower)]
  levels_upper <- cumsum(as.numeric(table(thetas_upper)))
  levels_upper <- levels_upper[-length(levels_upper)]
  thetas_lower <- unique(thetas_lower)
  thetas_upper <- unique(thetas_upper)
  
  simus$rating <- 1
  simus$rating[simus$response==1] <- as.numeric(as.character(cut(simus$conf[simus$response==1],
                                                                 breaks=thetas_upper, labels = levels_upper)))
  simus$rating[simus$response==-1] <- as.numeric(as.character(cut(simus$conf[simus$response==-1],
                                                                  breaks=thetas_lower, labels = levels_lower)))
  
  simus$correct <- as.numeric(simus$response == 1)
  return(simus)
}


  
    
