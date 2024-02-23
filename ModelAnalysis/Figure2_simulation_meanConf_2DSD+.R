### Simulation of the 2DSD+ model 
# This simulation should demonstrate that the 2DSD+ model (2DSD with a 
# time dependent confidence measure) is able to predict a double-increase model
# under certain parameter sets

# The Code produces Figure 2 in the paper
# and Suppl. Figure 3

# Preamble and imports                                     ----
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

library(dynConfiR)
library(tidyverse)
library(ggh4x)
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
dir.create("figures",showWarnings = FALSE)

# Figure 1                                                ----
set.seed(2201)
drift_levels <- seq(0, 3.2, length.out=5)
paramDf_tau <- expand.grid(sv=0.1, 
                           a=1.8, z=0.5, sz=0, t0=0, st0=0, 
                           tau=c(0.1, 0.4, 0.7, 1),
                           theta1=0, theta2=0.5, theta3=1, theta4=1.5,
                           lambda=c(0,0.5, 1, 2), 
                           model="2DSDT")
paramDf_tau[,paste("v", 1:5, sep="")] <- 
  rep(drift_levels, each=nrow(paramDf_tau))

simus_tau <- paramDf_tau %>% group_by(lambda, tau) %>%
  reframe(simulateRTConf(cbind(cur_group(), pick(everything())), 
                         n=10^5), .groups = "drop")

meanConf_tau <- simus_tau %>% group_by(lambda, tau, condition, correct) %>% 
  summarise(meanConf = mean(conf), .groups = "drop")
meanConf_tau <- meanConf_tau %>%
  mutate(condition = drift_levels[condition],
         lambda=paste("lambda==", format(lambda,nsmall=1), sep=""),
         tau = paste("tau==", format(tau,nsmall=1), sep=""))

unique(simus_tau$stimulus)
nrow(simus_tau) / 6 / 2 / 5
simus_tau %>% group_by(condition) %>% summarise(mRT = mean(rt), sdRT=sd(rt))
p_MRating_tau <- ggplot(meanConf_tau,
                        aes(x=condition, y=meanConf, group = as.factor(correct), 
                            color=as.factor(correct))) +
  geom_line(data=meanConf_tau, linewidth=1.5)+
  geom_point(fill="white", size=1.8)+
  facet_grid2(rows=vars(lambda), cols=vars(tau), 
              scales = "free_y", independent = "y",
              labeller = label_parsed)+ #dir="v"
  scale_x_continuous(name=expression(Mean~drift~rate~nu), 
                     breaks=drift_levels, labels = drift_levels)+
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
# # Only for manuscript generation
# ggsave("figures/simul_2DSD_taulambda.eps",
#        width = 17.62/2, height=8, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/simul_2DSD_taulambda.tiff", plot=p_MRating_tau,
       width = 17.62, height=15, units="cm",dpi=600)




# Supplementary Figure 1                                  ----
set.seed(2201)
drift_levels <- seq(0, 3.2, length.out=5)
paramDf_a <- expand.grid(sv=c(0.1, 0.5, 1, 1.5), 
                          a=c(0.5, 1,1.5, 2), z=0.5, sz=0, t0=0, st0=0, 
                         tau=0.1,
                         theta1=0, theta2=0.5, theta3=1, theta4=1.5,
                         lambda=2, 
                         model="2DSD")
paramDf_a[,paste("v", 1:5, sep="")] <- rep(drift_levels, each=nrow(paramDf_a))

simus_a <- paramDf_a %>% group_by(a, sv) %>%
  reframe(simulateRTConf(cbind(cur_group(), pick(everything())), n=10^5))

meanConf_a <- simus_a %>% group_by(a, sv, condition, correct) %>% 
  summarise(meanConf = mean(conf), .groups = "drop")
meanConf_a <- meanConf_a %>% mutate(condition = drift_levels[condition],
                                    sv = paste("s[nu]==", sv, sep=""),
                                    a = paste("a==", a, sep=""))

unique(simus_a$stimulus)
nrow(simus_a) / 6 / 2 / 5
simus_a %>% group_by(condition) %>% summarise(mRT = mean(rt), sdRT=sd(rt))
p_MRating_a <- ggplot(meanConf_a,
                      aes(x=condition, y=meanConf, group = as.factor(correct), 
                          color=as.factor(correct))) +
  geom_line(data=meanConf_a, linewidth=1.5)+
  geom_point(fill="white", size=1.8)+
  facet_grid2(rows=vars(sv), cols=vars(a), scales = "free_y", independent = "y",
              labeller = label_parsed)+ #dir="v"
  scale_x_continuous(name=expression(Mean~drift~rate~nu), breaks=drift_levels, labels = drift_levels)+
  scale_y_continuous(name="Mean confidence (arbitrary units)",
                     labels= function(x) format(x, nsmall=1))+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "",
                     labels=c("Correct", "Wrong")) +
  theme_bw() + expand_limits(y=1.2)+
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
# Only for manuscript generation
ggsave("C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/simul_2DSD_asv.eps",
       width = 17.62, height=15, units="cm",dpi=1200, device = cairo_ps)
ggsave("figures/simul_2DSD_ssv.tiff", plot = p_MRating_a,
       width = 17.62, height=15, units="cm",dpi=600)
