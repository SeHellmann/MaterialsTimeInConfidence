## Visualization of the mean posterior probability of being correct 
# given decision evidence, decision time, and visibility state
# Computation of the KL distance of visibility btw correct and incorrect decisions

# This code produces Figure 15 and
# Supplementary Figures 13 and 14 in the paper. 

# Preamble and imports                                     ----
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

library(dynConfiR)
library(tidyverse)
library(ggplot2)
library(ggh4x)
library(FNN)

odds_fct <- function(RT, X, V, Ds, a, z, sv, tau, svis, sigvis, ...) {
  SigV <- svis^2+(RT+tau)*sigvis^2
  Sc <- sv^2*(RT+tau)+1
  X <- -X
  TSSV <- (RT+tau)*(Sc+SigV)
  exp1 <- exp(-2*V*X/TSSV)
  v1<- rep(0, length(RT))
  v2 <- rep(0, length(RT))
  for (i in 1:length(Ds)) {
    v1 <- v1 + exp(-(TSSV*Ds[i]-V*Sc + X*SigV)^2/(2*Sc*SigV*TSSV))
    v2 <- v2 + exp(-(TSSV*Ds[i]-V*Sc - X*SigV)^2/(2*Sc*SigV*TSSV))
  }
  exp1*v1/v2
}

## Include font style and define colors for plots
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
custom_theme <- function() {
  font = "Times"
  theme_bw() %+replace%
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),    
          #legend.box.margin = margin(0,0,0,0,"cm"),
          text = element_text(size=9, family="Times"),
          axis.text = element_text(size=9, family="Times", color="black"),
          axis.title.y = element_text( angle=90, margin=margin(0,0.2, 0, 0 , unit="cm")),
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank(),
          strip.text = element_text(size=9, margin=margin(0.2, 0.2, 0.2, 0.2, "cm")),
          legend.margin = margin(0, 0, 0, 0, "cm"),#legend.direction = "horizontal",
          legend.key = element_blank(),
          legend.key.width=unit(2,"line"),
          panel.spacing=unit(0, "lines"))
}
dir.create("figures", showWarnings = FALSE)

##############################################################################
###### Figures 15: Simulated mean confidence and                    ##########
###### Suppl. Fig. 14: KL distance for visibility                   ##########
######       for different Vrange and sigvis(&svis)                 ##########
##############################################################################
paramDf <- as.data.frame(list(a=2, z=0.5, sz=0, sv=0,
                              tau=1.5, svis=1, sigvis=0.2,
                              theta1=0.5, lambda=0, t0=0, st0=0, w=0.5))
VRange <- matrix(data = c(0.65,0.65,
                          0.5, 0.8,
                          0.4, 0.9,
                          0.3, 1.1, 
                          0.1, 1.4, 
                            0,   2), ncol = 2, byrow = TRUE)
Vis_Var <- c(0.1, 0.4, 0.8)
vsteps <- 3
if (!file.exists("info_KLdistance_meanConf_simus.RData")) {
  set.seed(2201)
  res <- data.frame()
  meanconf <- data.frame()
  for (l in 1:length(Vis_Var)) {
    for (i in 1:nrow(VRange)) {
      vrange <- VRange[i,]
      vis_var <- Vis_Var[l]
      Ds <- seq(vrange[1],vrange[2], length.out=vsteps)
      paramDf <- as.data.frame(list(a=1.5, z=0.5, sz=0.7, sv=1.5,
                                    tau=0.3, svis=vis_var, sigvis=vis_var,
                                    theta1=0.5, lambda=0, t0=0, st0=0, w=0.5))
      paramDf[,paste("v", 1:vsteps, sep="")] <- Ds
      data <- simulateWEV(paramDf, 1e+6/(2*vsteps), model = "dynWEV", simult_conf = TRUE, process_results = TRUE)
      # ## To compute conditional KL distance for
      # ## visibility state cond. on decision state:
      # # Cut decision state in bins to integrate over
      # data$cut_dec <- cut(data$dec, breaks = 80) 
      # freq_dec <- table(data$cut_dec)
      # temp <- freq_dec[freq_dec > 300]
      # # Mark bins with too few observations
      # drop <- NULL
      # for (i in 1:length(temp)) {
      #   tempdata <- subset(data,cut_dec==names(temp)[i])
      #   if (min(nrow(subset(tempdata, correct==1)),nrow(subset(tempdata, correct==0)))< 20) {
      #     drop <- c(drop, i)
      #   }
      # }
      # if (!is.null(drop)) temp <- temp[-drop]
      # res_KL <- 0
      # for (i in 1:length(temp)) {
      #   tempdata <- subset(data,cut_dec==names(temp)[i])
      #   KLdist <- FNN::KL.dist(tempdata[tempdata$correct==1, "vis"], tempdata[tempdata$correct==0, "vis"])[10]
      #   res_KL <- res_KL + KLdist*temp[i]/sum(temp)
      # }
      KL_vis <- FNN::KL.divergence(data[data$correct==1, "vis"], data[data$correct==0, "vis"], k=15)[15]
      data$odds <- odds_fct(data$rt, data$dec, data$vis, Ds = Ds,
                            paramDf$a, paramDf$z, paramDf$sv, paramDf$tau, paramDf$svis,
                            paramDf$sigvis)
      data$prob <- data$odds/(1+data$odds)
      data$d <- Ds[data$condition]
      temp_meanconf <- data %>% group_by(d, correct) %>%
        summarise(Mconf = mean(prob)) %>%
        mutate(Vrange=diff(vrange), vis_var=vis_var)
      meanconf <- rbind(meanconf, temp_meanconf)
      res <- rbind(res, c(diff(vrange), vis_var, KL_vis, mean(data$correct))) #, KL_dec, res_KL
    }
  }
  names(res) <- c("Vrange", "Visvar", "KL_vis", "Accuracy") #, "KL_dec", "KL_vis_cond_dec"
  save(res, meanconf, file="info_KLdistance_meanConf_simus.RData")
} else {
  load("info_KLdistance_meanConf_simus.RData")
}
p_MeanOptConf <- ggplot(subset(meanconf,Vrange>0.4), aes(x=d, y=Mconf, color=as.factor(correct)))+
  geom_line(linewidth=1.2)+
  xlab("Discriminability d")+ylab("Mean confidence")+
  scale_color_manual(name="",values=two_colors_correct,
                     breaks=c("1","0"), labels=c("Correct", "Incorrect"))+
  facet_nested(rows=c(vars("Range of discriminability"),vars(Vrange)),
               cols=c(vars("Variability in visibility accumulation"), vars(vis_var)),
               strip=strip_nested(by_layer_y = TRUE,
                            text_y = elem_list_text(angle=c(-90,0))))+
  custom_theme()+theme(legend.position = "bottom")
p_MeanOptConf
# # Only for manuscript generation
ggsave("figures/MeanOptimalConf.eps",
       width = 17.62, height=15, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/MeanOptimalConf.tiff", plot = p_MeanOptConf, 
       width = 17.62, height=15, units="cm",dpi=600)

p_KLvisibility <- ggplot(res, aes(x=Vrange, y=KL_vis, 
                                  linetype=as.factor(Visvar),
                                  color=as.factor(Visvar)))+
  geom_line(linewidth=1.2)+
  ylab("KL divergence")+xlab("Range of discriminability values")+
  scale_linetype_discrete(name="Variation in\nvisibility\naccumulation")+
  scale_color_viridis_d(name="Variation in\nvisibility\naccumulation")+
  custom_theme()+theme(legend.position = "right",
                       legend.direction = "vertical")
p_KLvisibility
# ## Only for manuscript generation
# ggsave("C:/Users/go73jec/Documents/Forschung/Manuskripte/TimeInDynWEV/Supplement/figures/KLvisibility.eps",
#        width = 12, height=5, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/KLvisibility.tiff", plot=p_KLvisibility,
       width = 12, height=5, units="cm",dpi=600)




##########################################################################
###### SupplFig 13: Mean Opt Conf Simulation tau and sv   ################
##########################################################################
Tau <- c(0.1, 0.5, 1)
mu_Var <- c(0.4, 0.8, 1.2)
i <- 1
if (!file.exists("info_KLdistance_meanConf_tausv.RData")) {
  set.seed(2201)
  meanconf <- data.frame()
  for (l in 1:length(mu_Var)) {
    for (i in 1:length(Tau)) {
      vrange <- c(0.2, 1.7)
      vsteps <- 4
      vis_var <- 0.1
      tau <- Tau[i]
      mu_var <- mu_Var[l]
      Ds <- seq(vrange[1],vrange[2], length.out=vsteps)
      paramDf <- as.data.frame(list(a=1.5, z=0.5, sz=0, sv=mu_var,
                                    tau=tau, svis=vis_var, sigvis=vis_var,
                                    theta1=0.5, lambda=0, t0=0, st0=0, w=0.5))
      paramDf[,paste("v", 1:vsteps, sep="")] <- Ds
      data <- simulateWEV(paramDf, 1e+6/(2*vsteps), model = "dynWEV", simult_conf = TRUE, process_results = TRUE)
      data$odds <- odds_fct(data$rt, data$dec, data$vis, Ds = Ds,
                            paramDf$a, paramDf$z, paramDf$sv, paramDf$tau, paramDf$svis,
                            paramDf$sigvis)
      data$prob <- data$odds/(1+data$odds)
      data$d <- Ds[data$condition]
      temp_meanconf <- data %>% group_by(d, correct) %>%
        summarise(Mconf = mean(prob)) %>%
        mutate(Vrange=diff(vrange), tau=tau, mu_var=mu_var)
      meanconf <- rbind(meanconf, temp_meanconf)
    }
  }
  save(meanconf, file="info_KLdistance_meanConf_tausv.RData")
} else {
  load("info_KLdistance_meanConf_tausv.RData")
}
p_meanconf_tausv <- ggplot(meanconf, aes(x=d, y=Mconf, color=as.factor(correct)))+
  geom_line(linewidth=1.2)+
  xlab("Discriminability d")+ylab("Mean confidence")+
  scale_color_manual(name="",values=two_colors_correct,
                     breaks=c("1","0"), labels=c("Correct", "Incorrect"))+
  facet_nested(rows=c(vars("Between-trial variability in drift rate"),vars(mu_var)),
               cols=c(vars("Postdecisional accumulation time"), vars(tau)))+
  custom_theme()+theme(legend.position = "bottom")
# # Only for manuscript generation
# ggsave("C:/Users/go73jec/Documents/Forschung/Manuskripte/TimeInDynWEV/Supplement/figures/MeanOptimalConftausv.eps",
#       width = 15, height=15, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/MeanOptimalConftausv.tiff", plot=p_meanconf_tausv,
       width = 17.62, height=15, units="cm",dpi=600)
