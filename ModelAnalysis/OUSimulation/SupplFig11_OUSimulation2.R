## Simulation of mean confidence in an Ornstein-Uhlenbeck based 
# post-decisional accumulation model
# This simulation should demonstrate that the OU-model with time-dependency 
# may produce a double-increase pattern

# This code produces Supplementary Figures 11.

# Preamble and imports    
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())
Rcpp::sourceCpp("SimulateOUConf.cpp")
library(tidyverse)
library(ggpubr)
library(ggh4x)
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
dir.create("figures",showWarnings = FALSE)

#  Parameters as reported by Yu et al. (2015) (https://doi.org/10.1037/xge0000062) in Study2     
#             a,   z, sz,  v,  sv,      k,   tau, t0,  st0    s,  lambda)
params <- c(0.9, 0.5, 0.0, 2,   0, 0.473, 1.0, 0.2, 0.05, 0.1,      0)
targetacc <- c(.578, .662, .791, .942, .99)

if (!file.exists("saved_simulationsOUlambda.RDATA")) {
  # ## First, we find suitable boundary separation and drift rate values
  # ## to produce a reasonable range of accuracy values similar to the ones reported in the study.
  # SSE_accs <- function(X) {
  #   params[1] <- X[1]
  #   V <- X[2:6]
  # 
  #   n <- 3000
  #   simus <- data.frame()
  #   for (i in 1:length(V)) {
  #     params[4] <- V[i]
  #     tempsimus <- cbind(SimulateOUConfidence(n, params),
  #                        condition=i, v=V[i])
  #     simus <- rbind(simus, tempsimus)
  #   }
  #   names(simus) <- c("rt", "response", "conf", "condition", "v")
  #   propnonfinished <- mean(simus$response==0)
  #   if (propnonfinished>.3) return(propnonfinished*20)
  #   simus <- subset(simus, response !=0)
  # 
  #   aggsimus <- simus %>% group_by(condition) %>%
  #     summarise(acc = mean(response==1))
  #   return(sum((aggsimus$acc-targetacc)^2))#+propnonfinished*5)
  # }
  # set.seed(24200801)
  # V  <- c(0.03, 0.06, 0.1, 0.2, 0.3)
  # opt_a <-  optim(c(0.1, V), SSE_accs,
  #                 lower = c(1e-12, rep(0.0001, 5)),
  #                 upper=c(Inf, rep(Inf, 5)), method = "L-BFGS-B")#  control = list(factr=1e12))
  # opt_a
  # a <- opt_a$par[1]
  # params[1] <- a
  # drift_levels <- round(opt_a$par[2:6], 2)
  # dput(drift_levels)
  # a
  
  ## Results:
  params[1] <- 0.13
  drift_levels <- c(0.03, 0.05, 0.1, 0.2, 0.32)
  
  set.seed(2201)
  n = 3*10^6
  params[5] <- 0
  simus <- data.frame()
  Tau <- c(0.1, 0.4, 0.7)
  Lambda <-  c(0.5, 2, 4)
  for (k in 1:length(Tau)) {
    params[7] <- Tau[k]
    for (j in 1:length(Lambda)) {
      params[11] <- Lambda[j]
      for (i in 1:length(V)) {
        params[4] <- V[i]
        tempsimus <- SimulateMeanOUConfidence(n, params)
        tempsimus <- cbind(tempsimus, correct=c(1,0),
                           condition=i, v=V[i], tau = Tau[k], lambda=Lambda[j])
        simus <- rbind(simus, tempsimus)
      }
    }
  }
  names(simus) <- c("response", "meanConf", "count", "correct",  "condition", "v", "tau", "lambda")
  save(simus,drift_levels, params, file="saved_simulationsOUlambda.RDATA")
} else {
  load("saved_simulationsOUlambda.RDATA")
}
meanConf_tau <- simus
meanConf_tau <- meanConf_tau %>%
  mutate(condition = drift_levels[condition],
         lambda=paste("lambda==", format(lambda,nsmall=1), sep=""),
         tau = paste("tau==", format(tau,nsmall=1), sep=""))

nrow(simus) / 5 / 3 / 3
#simus %>% group_by(condition) %>% summarise(mRT = mean(rt), sdRT=sd(rt))

p_MRating_tau <- ggplot(meanConf_tau,
                        aes(x=as.factor(condition), y=meanConf, group = as.factor(correct), 
                            color=as.factor(correct))) +
  geom_line(data=meanConf_tau, linewidth=1.5)+
  geom_point(fill="white", size=1.8)+
  facet_grid2(rows=vars(lambda), cols=vars(tau), 
              scales = "free_y", independent = "y",
              labeller = label_parsed)+ #dir="v"
  scale_x_discrete(name=expression(Mean~drift~rate~nu), 
                     labels = drift_levels)+
  scale_y_continuous(name="Mean confidence (arbitrary units)",
                     labels= function(x) format(x, nsmall=1))+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "",
                     labels=c("Correct", "Wrong")) +
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position="bottom",
        legend.direction = "horizontal", 
        legend.title = element_blank(),
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.margin = margin(0, 0, 0, 0, "cm"),
        legend.box.margin = margin(-0.3, 0, 0, 0, "cm"),
        legend.key.width=unit(1.5,"line"),
        panel.spacing.y=unit(0, "lines"))
p_MRating_tau
# # # Only for manuscript generation
# ggsave("../../../Supplement/figures/simul_2SOU_taulambda.eps",
#        width = 17.62, height=10, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/simul_2SOU_taulambda.tiff", plot=p_MRating_tau,
       width = 17.62, height=15, units="cm",dpi=300)

simus %>% group_by(tau, lambda) %>%
  mutate(Nobs=sum(count)) %>%
  group_by(condition, tau, lambda) %>%
  reframe(acc = count[1]/sum(count)) %>% 
  group_by(condition) %>%
  reframe(MAcc=mean(acc), SDAcc=sd(acc))
