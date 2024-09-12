############################################################---
#####    Additional script for the analysis of dynWEV with ----
#####         scaled parameters                            ----
#__________________________________________________________----

# Sebastian Hellmann, 23.02.2024

## Structure:
# Preamble and imports                                   
# A  Load previous results
# B  Scale dynWEV parameters and compare with dynaViTE parameters
# Not: C  Predict distribution and compare to other dynWEV predicitons

# Preamble and imports                                     ----
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

{
  library(Rmisc)  # for summarySEwithin
  #library(Hmisc)
  #library(ggpubr) # for ggarrange
  library(dynConfiR)  # for fitting and predicting the models
  library(tidyverse)
  #library(gridExtra)
  #library(grid)
  library(BayesFactor) 
  #library(xtable)   # for latex tables
  library(ggh4x)  # for facet_nested
  library(parallel)
  # To use Times as font on Windows (otherwise ignore the ggplot warnings)
  windowsFonts(Times=windowsFont("Times New Roman"))
}


# A  Load previous results
load("../modelfits_threeexps.RData")

experiment_colors = c("#AA4499","#117733", "#DDCC77")
par_labels <- c('nu[1]', 'nu[2]', 'nu[3]', 'nu[4]', 'nu[5]', 's[nu]', 'a', 'z', 's[z]', 
                't[0]', 's[t0]', 'tau', 'w', 'sigma[V]', 's[V]',
                paste0('theta[', -1,'~', 1:4,']'),#'',
                paste0('theta[', 1,'~',  1:4,']'))
par_levels <- c('v1', 'v2', 'v3', 'v4', 'v5', 'sv', 'a', 'z', 'sz',
                't0', 'st0', 'tau', 'w', 'sigvis', 'svis',
                paste0('thetaLower',1:4),
                paste0('thetaUpper',1:4))
parfits_long_dynaViTEdynWEV <- fits %>% filter(model %in% c("dynaViTE", "dynWEV")) %>%
  select(model, participant, experiment, t0:sigvis, v4, v5 ) %>%
  pivot_longer(cols = t0:v5, names_to ="parameter", values_to="fit") %>%
  pivot_wider(names_from = model, values_from=fit) %>%
  mutate(parameter = factor(parameter, levels = par_levels,labels = par_labels))
DDMparfits_long_dynaViTEdynWEV <- parfits_long_dynaViTEdynWEV  %>%
  filter(!grepl("theta|V|w|tau", parameter)) 

est_cors <- subset(DDMparfits_long_dynaViTEdynWEV , abs(dynWEV) < 1e+6 & abs(dynaViTE) < 1e+6) %>%
  group_by(parameter) %>%
  reframe(cor = cor(dynWEV, dynaViTE),
          x_min=max(dynWEV), 
          y_max=min(dynaViTE))

DDM_sscaled_pars <- DDMparfits_long_dynaViTEdynWEV %>% 
  filter(!grepl("z|t", parameter) & !is.na(dynaViTE))
sscaled_pars_slope <-  lm(dynaViTE~dynWEV-1, data=DDM_sscaled_pars)$coefficients
slopes <- data.frame(parameter=unique(DDMparfits_long_dynaViTEdynWEV$parameter), 
                     slope=c(NA, sscaled_pars_slope)[(!grepl("z|t", unique(DDMparfits_long_dynaViTEdynWEV$parameter)))+1])
round(sscaled_pars_slope, 2)




# B  Scale dynWEV parameters and compare with dynaViTE parameters
sscaled_pars_slope
## Following parameters are affected by scaling s 
scaled_dynWEV_parameters <- c('v1', 'v2', 'v3', 'v4', 'v5', 'sv', 'a',
                              'sigvis', 'svis',
                              paste0('thetaLower',1:4),
                              paste0('thetaUpper',1:4))
scaled_dynWEV_fits <- subset(fits, model=="dynWEV") %>% 
  mutate(across(scaled_dynWEV_parameters, ~ .x*sscaled_pars_slope)) %>%
  mutate(model="dynWEV(scaled)")
# temp2 <- subset(fits, model=="dynWEV")
# head(temp1)
# head(temp2)


## Suppl. Figure 10: Cor DDM parameters btw dynaViTE and scaled dynWEV   ----

parfits_long_dynaViTEdynWEV <- rbind(fits, scaled_dynWEV_fits) %>% filter(model %in% c("dynaViTE", "dynWEV(scaled)")) %>%
  select(model, participant, experiment, t0:sigvis, v4, v5 ) %>%
  pivot_longer(cols = t0:v5, names_to ="parameter", values_to="fit") %>%
  pivot_wider(names_from = model, values_from=fit) %>%
  mutate(parameter = factor(parameter, levels = par_levels,labels = par_labels))

DDMparfits_long_dynaViTEdynWEV_scaled <- parfits_long_dynaViTEdynWEV  %>%
  filter(!grepl("theta|V|w|tau", parameter)) 

est_cors <- subset(DDMparfits_long_dynaViTEdynWEV_scaled , abs(`dynWEV(scaled)`) < 1e+6 & abs(dynaViTE) < 1e+6) %>%
  group_by(parameter) %>%
  reframe(cor = cor(`dynWEV(scaled)`, dynaViTE),
          x_min=max(`dynWEV(scaled)`), 
          y_max=min(dynaViTE))

DDM_sscaled_pars <- DDMparfits_long_dynaViTEdynWEV_scaled %>% 
  filter(!grepl("z|t", parameter) & !is.na(dynaViTE))
sscaled_pars_slope_corrected <-  lm(dynaViTE~`dynWEV(scaled)`-1, data=DDM_sscaled_pars)$coefficients
slopes <- data.frame(parameter=unique(DDMparfits_long_dynaViTEdynWEV$parameter), 
                     slope=c(NA, sscaled_pars_slope_corrected)[(!grepl("z|t", unique(DDMparfits_long_dynaViTEdynWEV$parameter)))+1])
round(sscaled_pars_slope, 2)
ggplot(filter(DDMparfits_long_dynaViTEdynWEV_scaled, abs(`dynWEV(scaled)`)<1e+6 & abs(dynaViTE) < 1e+6), 
       aes(x=`dynWEV(scaled)`, y=dynaViTE, color=experiment))+
  geom_smooth(data=filter(DDMparfits_long_dynaViTEdynWEV_scaled, abs(`dynWEV(scaled)`)<1e+6 & abs(dynaViTE) < 1e+6),  
              method="lm", alpha=0.5, linewidth=0.9, 
              aes(x=`dynWEV(scaled)`, y=dynaViTE), formula = 'y~x-1', inherit.aes = FALSE)+
  facet_wrap(.~parameter, scales = "free", labeller = label_parsed, 
             drop = TRUE, ncol=4)+
  geom_abline(col="black", linewidth=1, alpha=0.2)+
  geom_abline(data=slopes, aes(slope=slope, intercept=0),col="red",linewidth=1, alpha=0.2)+
  scale_color_manual(name="", values=experiment_colors)+
  geom_point(alpha=0.5)+ scale_shape_discrete("")+
  geom_text(data = est_cors, #mutate(est_cors, experiment="Experiment 1"),
            aes(x=x_min, y=y_max,
                label=paste0(format(round(cor, 2), nsmall=2))),#"rho == ",
            color="black",
            parse=TRUE, hjust=1, vjust=0)+
  theme_bw()+ylab("Fitted parameter for dynaViTE")+xlab("Scaled parameter fits for dynWEV")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = c(0.78, 0.28), legend.justification = c(0, 1),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9))
ggplot(filter(DDM_sscaled_pars, abs(`dynWEV(scaled)`)<1e+6 & abs(dynaViTE) < 1e+6), 
       aes(x=`dynWEV(scaled)`, y=dynaViTE, color=experiment))+
  geom_smooth(data=filter(DDM_sscaled_pars, abs(`dynWEV(scaled)`)<1e+6 & abs(dynaViTE) < 1e+6),  
              method="lm", alpha=0.5, linewidth=0.9, 
              aes(x=`dynWEV(scaled)`, y=dynaViTE), formula = 'y~x-1', inherit.aes = FALSE)+
  geom_abline(col="red", linewidth=1, alpha=0.2)+
  #geom_abline(data=slopes, aes(slope=slope, intercept=0),col="red",linewidth=1, alpha=0.2)+
  scale_color_manual(name="", values=experiment_colors)+
  geom_point(alpha=0.5)+ scale_shape_discrete("")+
  geom_text(data = filter(est_cors,!grepl("z|t", parameter)), #mutate(est_cors, experiment="Experiment 1"),
            aes(x=x_min, y=y_max,
                label=paste0(format(round(cor, 2), nsmall=2))),#"rho == ",
            color="black",
            parse=TRUE, hjust=1, vjust=0)+
  facet_wrap(.~parameter, scales = "free", labeller = label_parsed, 
             drop = TRUE, ncol=4)+
  theme_bw()+ylab("Fitted parameters for dynaViTE")+xlab("Scaled parameters for dynWEV")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = c(0.73, 0.29), legend.justification = c(0, 1),
        legend.text = element_text(size=8),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9),
        legend.background = element_blank())
dir.create("figures", showWarnings = FALSE)
ggsave("figures/parametersDynWEVdynaViTE_scaled.tiff",
       width = 17.625, height=10, units="cm",dpi=600)
# # Only for manuscript generation
ggsave("../../../Supplement/figures/parametersDynWEVdynaViTE_scaled.eps",
       width = 15, height=8, units="cm",dpi=600, device = cairo_ps)



# # C  Predict distribution and compare to other dynWEV predicitons
# scaled_dynWEV_fits$s <- sscaled_pars_slope
# if (!file.exists("raw_predictions_scaleddynWEV.RData")) {
#   preds_confidence_scaleddynWEV <- predictConfModels(scaled_dynWEV_fits, maxrt = 18, subdivisions = 250, 
#                                         simult_conf = TRUE, parallel=TRUE)
#   preds_RT_scaleddynWEV <- predictRTModels(scaled_dynWEV_fits, maxrt = 13, subdivisions = 200, minrt = 0, 
#                               simult_conf=TRUE, scaled=TRUE, 
#                               DistConf=preds_confidence_scaleddynWEV, parallel=TRUE)
#   preds_confidence_scaleddynWEV$experiment <- if_else(preds_confidence_scaleddynWEV$participant > 300, 
#                                          "Shekhar & Rahnev (2021)\nExperiment 4", 
#                                          if_else(preds_confidence_scaleddynWEV$participant > 200, 
#                                                  "Hellmann et al. (2023)\nExperiment 1", 
#                                                  "Hellmann et al. (2023)\nExperiment 2"))
#   preds_RT_scaleddynWEV$experiment <- if_else(preds_RT_scaleddynWEV$participant > 300, 
#                                  "Shekhar & Rahnev (2021)\nExperiment 4", 
#                                  if_else(preds_RT_scaleddynWEV$participant > 200, 
#                                          "Hellmann et al. (2023)\nExperiment 1", 
#                                          "Hellmann et al. (2023)\nExperiment 2"))
#   save(preds_confidence_scaleddynWEV, preds_RT_scaleddynWEV, file="raw_predictions_scaleddynWEV.RData")
# } else {
#   load("raw_predictions_scaleddynWEV.RData")
# }
# load("../collected_raw_data.RData")
# load("../raw_predictions.RData")
# Ns_part <- Data %>% group_by(participant) %>% 
#   summarise(N=n(), MinRT = min(rt), .groups = "drop")  %>%
#   select(participant, N)
# Preds_RTdens_dynWEV <-  preds_RT %>% 
#   left_join(Ns_part, by="participant") %>%
#   ungroup() %>% filter(model=="dynWEV") %>%
#   select(-participant) %>%
#   group_by(rating, condition, stimulus, response, experiment, correct, rt) %>%
#   summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data), .groups = "drop") %>%
#   mutate(condition = factor(condition+ifelse(experiment=="Hellmann et al. (2023)\nExperiment 2", 5, 0)+
#                               ifelse(experiment=="Shekhar & Rahnev (2021)\nExperiment 4", 10, 0),
#                             levels = 1:13, 
#                             labels = levels(Data$condition)))
# 
# Preds_RTdens_dynWEV_scaled <-  preds_RT_scaleddynWEV %>% 
#   left_join(Ns_part, by="participant") %>%
#   ungroup() %>% 
#   select(-participant) %>%
#   group_by(rating, condition, stimulus, response, experiment, correct, rt) %>%
#   summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data), .groups = "drop") %>%
#   mutate(condition = factor(condition+ifelse(experiment=="Hellmann et al. (2023)\nExperiment 2", 5, 0)+
#                               ifelse(experiment=="Shekhar & Rahnev (2021)\nExperiment 4", 10, 0),
#                             levels = 1:13, 
#                             labels = levels(Data$condition)))
# Preds_RTdens_dynWEV_scaled$class <- "scaled"
# Preds_RTdens_dynWEV <-Preds_RTdens_dynWEV %>%mutate(class="fitted") %>%
#   rbind(Preds_RTdens_dynWEV_scaled)
# 
# ggplot(Preds_RTdens_dynWEV, aes(x=rt, y=dens, shape=class, color=condition))+
#   geom_line()+
#   facet_nested(cols=c(vars(stimulus), vars(response)),rows=vars(rating), 
#                scales="free", independent = "y")
# 
# temp <- Preds_RTdens_dynWEV %>% 
#   select(-densscaled) %>%
#   pivot_wider(names_from = class, values_from = dens)
# all.equal(temp$fitted, temp$scaled)
# 
# 
# 
# Preds_RTdens_dynWEV2 <- Preds_RTdens_dynWEV %>%
#   group_by(class, correct, rating, rt, experiment, condition) %>%
#   summarise(dens=mean(dens))
# ggplot(Preds_RTdens_dynWEV2, aes(x=rt, y=dens, linetype=class, color=condition))+
#   geom_line()+ ylab("Predicted defective response time density")+xlab("Response time")+
#   xlim(c(0, 8))+
#   facet_nested(cols=c(vars(experiment), vars(correct)),rows=vars(rating), 
#                scales="free", independent = "y")
