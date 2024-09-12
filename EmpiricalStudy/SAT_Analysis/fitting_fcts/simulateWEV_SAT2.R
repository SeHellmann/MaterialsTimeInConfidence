
## When given vectorised parameters, n is the number of replicates for each parameter set
#' @rdname simulateWEV
#' @export
simulateWEV_SAT <- function (paramDf, n=1e+4,  df, model = "dynWEV",
                         delta=0.01, maxrt=15, seed=NULL)
{
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (!("lambda" %in% names(paramDf))) {
    if (model %in% c("dynaViTE", "2DSDT")) warning("No lambda specified in paramDf. lambda=0 used")
    lambda <- 0
  } else {
    lambda <- paramDf$lambda
  }
  # if (model %in% c("WEVmu", "dynWEV")) model <- "dynaViTE"
  # if (grepl("2DSD", model)) model <- "2DSD"
  # if ("model" %in% names(paramDf)) {
  #   model <- paramDf$model
  #   paramDf$model <- NULL
  # 
  # }
  # if (!(model %in% c("dynaViTE", "2DSD"))) stop("Only models dynaViTE, dynWEV (alias: WEVmu), and 2DSD/2DSDT are allowed!")

  ## recover parameters from paramDf
  #z <- paramDf$z
  sz <- paramDf$sz
  t0 <- paramDf$t0  # recalc_t0 (see e.g. dWEV)
  st0 <- paramDf$st0
  tau0 = paramDf$tau0
  # if ("d" %in% names(paramDf)) {
  #   d <- paramDf$d
  # } else {
  #   d <- 0
  # }

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
  # vary_sv <-   length(grep(pattern = "^sv[0-9]", names(paramDf), value = T))>1
  # if (vary_sv){
  #   SV <- c(t((paramDf[,paste("sv",1:(nConds), sep = "")])))
  # } else {
  #   SV <- rep(paramDf$sv, nConds)
  # }
  # vary_s <-   length(grep(pattern = "^s[0-9]", names(paramDf), value = T))>1
  # if (vary_s){
  #   S <- c(t((paramDf[,paste("s",1:(nConds), sep = "")])))
  # } else {
  #   if ("s" %in% names(paramDf)) {
  #     S <- rep(paramDf$s, nConds)
  #   } else {
  #     S <- rep(1, nConds)
  #   }
  # }
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
                                 sigvis, svis)
  #print(params)
  
  simus <- res %>% group_by(SAT, condition) %>%
    reframe(as.data.frame(RNG_WEV_SAT(n, .data$a[1], .data$v[1], 
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
