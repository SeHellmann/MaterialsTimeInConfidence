############################################################-
#####    Main script for the empirical analyses in          -
#####       Time in dynamical confidence models             -
#__________________________________________________________----

# Sebastian Hellmann, 21.9.2023

## Structure:
# Preamble and imports                                   
# A  Model fitting and Data Preparation                  
## 1. Read experimental data                             
## 2. Fit the four models to the three experimental data 
## 3. Compute predicted distributions from model fits    
## 4. Aggregate data and predictions                     
# B  Generate visualizations                             
# C  Quantitative comparison                             
# D  dynaViTE and dynWEV parameter analysis              
# E  Parameter and model recovery                       
# F  Discussion and Supplement                           


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
}



#__________________________________________________________----
# A  Model fitting and Data Preparation                    ----
## 1. Read experimental data                               ----
if (!file.exists("collected_raw_data.RData")) {
  load("data/Data_KonfMaskTime.RData")
  Data_KonfMask <- Data
  load("data/Data_RDK.RData")
  Data_RDK <- Data %>% select(c("participant", "condition", "correct", "rt", "stimulus", "response", "rating")) %>%
    mutate(experiment="Hellmann et al. (2023)\nExperiment 2")
  Data_KonfMask <- Data_KonfMask %>% select(c("participant", "condition", "correct", "rt", "stimulus", "response", "rating")) %>%
    mutate(experiment="Hellmann et al. (2023)\nExperiment 1",
           participant = participant + 200)
  load("data/Shekar2021.RData")
  Data_Shekhar <- Data %>% select(c("participant", "condition", "correct", "rt", "stimulus", "response", "rating")) %>%
    mutate(experiment="Shekhar & Rahnev (2021)\nExperiment 4",
           participant = participant + 300)
  summary(Data_KonfMask)
  summary(Data_RDK)
  summary(Data_Shekhar)
  cond_levels <- c(sort(unique(Data_KonfMask$condition)),
                   sort(unique(Data_RDK$condition)),
                   sort(unique(Data_Shekhar$condition)))
  Data <- rbind(Data_RDK,Data_KonfMask, Data_Shekhar) %>% 
    mutate(condition = factor(condition,
                              levels = cond_levels, 
                              labels = c("8.3", "16.7", "33.3", "66.7", "133.3",
                                         "1.6", "3.2", "6.4", "12.8", "25.6",
                                         "4.5", "6", "8")))
  save(Data, Data_KonfMask, Data_RDK, Data_Shekhar, 
       file="collected_raw_data.RData")
} else {
  load("collected_raw_data.RData")
}
## 2. Fit the four models to the three experimental data   ----
if (!file.exists("modelfits_threeexps.RData")) {
  fits_RDK <- fitRTConfModels(Data_RDK, 
                              models=c("dynaViTE", "2DSDT", "dynWEV", "2DSD"), 
                              nRatings=5, fixed=list(sym_thetas=FALSE), 
                              restr_tau = "simult_conf", logging=TRUE,
                              opts=list(nAttempts=4),
                              parallel = "both", n.cores=c(5,4))
  fits_KonfMask <- fitRTConfModels(Data_KonfMask, 
                                   models=c("dynaViTE", "2DSDT", "dynWEV", "2DSD"), 
                                   nRatings=5, fixed=list(sym_thetas=FALSE), 
                                   restr_tau = "simult_conf", logging=TRUE,
                                   opts=list(nAttempts=4),
                                   parallel = "both", n.cores=c(5,4))
  fits_Shekhar <- fitRTConfModels(Data_Shekhar,
                                  models=c("dynaViTE", "2DSDT", "dynWEV", "2DSD"), 
                                  nRatings=5, fixed=list(sym_thetas=FALSE), 
                                  restr_tau = "simult_conf", logging=TRUE,
                                  opts=list(nAttempts=4),
                                  parallel = "both", n.cores=c(5,4))
  
  fits_Shekhar$experiment = "Shekhar & Rahnev (2021)\nExperiment 4"
  fits_KonfMask$experiment="Hellmann et al. (2023)\nExperiment 1"
  fits_RDK$experiment = "Hellmann et al. (2023)\nExperiment 2"
  fits_Shekhar$v4 <- NA
  fits_Shekhar$v5 <- NA
  
  fits <- rbind(fits_Shekhar, fits_RDK, fits_KonfMask)
  save(fits, file="modelfits_threeexps.RData")
} else {
  load("modelfits_threeexps.RData")
}
rm(list=c("Data_KonfMask", "Data_RDK", "Data_Shekhar"))

## 3. Compute predicted distributions from model fits      ----
if (!file.exists("raw_predictions.RData")) {
  preds_confidence <- predictConfModels(fits, maxrt = 18, subdivisions = 250, 
                                        simult_conf = TRUE, parallel=TRUE)
  preds_RT <- predictRTModels(fits, maxrt = 13, subdivisions = 200, minrt = 0, 
                              simult_conf=TRUE, scaled=TRUE, 
                              DistConf=preds_confidence, parallel=TRUE)
  preds_confidence$experiment <- if_else(preds_confidence$participant > 300, 
                                         "Shekhar & Rahnev (2021)\nExperiment 4", 
                                         if_else(preds_confidence$participant > 200, 
                                                 "Hellmann et al. (2023)\nExperiment 1", 
                                                 "Hellmann et al. (2023)\nExperiment 2"))
  preds_RT$experiment <- if_else(preds_RT$participant > 300, 
                                 "Shekhar & Rahnev (2021)\nExperiment 4", 
                                 if_else(preds_RT$participant > 200, 
                                         "Hellmann et al. (2023)\nExperiment 1", 
                                         "Hellmann et al. (2023)\nExperiment 2"))
  save(preds_confidence, preds_RT, file="raw_predictions.RData")
} else {
  load("raw_predictions.RData")
}
## 4. Aggregate data and predictions                       ----
if (!file.exists("collected_Data_Fits_Predicts.RData")) {
  ### a) Data confidence rating distribution               ----
  Data <- Data %>% group_by(experiment, participant, condition) %>%
    mutate(nrows=n())
  Data_RatingDist_part <- Data %>% 
    group_by(experiment, participant, condition, rating, correct) %>% 
    summarise(p = n()/(mean(nrows)),.groups = "drop") %>%
    full_join(y = merge(distinct(Data[,c("correct","rating", "experiment")]),
                        distinct(Data[,c("condition", "participant", "experiment")]), by="experiment"),
              by = c("experiment", "participant", "condition", "rating", "correct")) %>%
    mutate(p = ifelse(is.na(p), 0, p)) 

  ### Sanity Checks:
  # sum(Data_RatingDist_part$p)
  # table((Data_RatingDist_part %>% group_by(condition, participant) %>% summarise(p = sum(p)))$p)
  
  # For the plots we won't differentiate between 
  # stimulus directions and participants
  Data_RatingDist_corr_cond <- Data_RatingDist_part %>% 
    group_by(experiment, correct, condition, rating) %>% 
    summarise(p = mean(p), .groups = "drop")
  #sum(Data_RatingDist_corr_cond$p)
  
  ### b) Data mean confidence rating                       ----
  Data_MRating_corr_cond <- Data %>% 
    group_by(participant, experiment, condition, correct) %>%
    summarise(MRating = mean(rating)) 
  Data_MRating_corr_cond_SE <- Data_MRating_corr_cond %>%
    group_by(experiment) %>%
    reframe(summarySEwithin(pick(everything()),measurevar = "MRating", 
                            idvar = "participant", withinvars = c("condition", "correct"))) %>%
    select(experiment, condition, correct, sewithin = se) 
  Data_MRating_corr_cond <- Data_MRating_corr_cond %>%
    group_by(experiment, condition, correct) %>%
    summarise(MRating=mean(MRating, na.rm=TRUE), .groups="drop") %>%
    mutate(correct=as.factor(correct)) %>%
    left_join(Data_MRating_corr_cond_SE, by=c("condition", "correct", "experiment"))
  ### c) Data response time quantiles                      ----  
  # Reaction Time Quantiles of the Data grouped by rating and accuracy  
  Data_RTQuants_corr_rating <- Data %>%
    group_by(experiment, rating, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  # Reaction Time Quantiles of the Data grouped by condition and accuracy  
  Data_RTQuants_corr_cond <- Data %>%
    group_by(experiment, condition, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  
  
  ### d) Prediction confidence rating distribution         ----
  Preds_RatingDist_corr_cond <- preds_confidence %>%
    group_by(model, experiment, rating, correct, condition) %>%
    summarise(p = mean(p), .groups="drop") %>%
    mutate(model = factor(model, levels = c("dynaViTE", "2DSDT", "dynWEV", "2DSD"),
                          labels = c("dynaViTE", "2DSDT", "dynWEV", "2DSD")), 
           condition = factor(condition+ifelse(experiment=="Hellmann et al. (2023)\nExperiment 2", 5, 0)+
                                ifelse(experiment=="Shekhar & Rahnev (2021)\nExperiment 4", 10, 0),
                              levels = 1:13, 
                              labels = levels(Data$condition)))
  # # Sanity checks:
  # Preds_RatingDist_corr_cond %>% group_by(model, experiment,condition) %>%
  #    summarise(p=sum(p))
  
  ### e) Prediction mean confidence rating                 ----
  # This is good to visualize a folded-X- and double-increase-pattern 
  Preds_MRating_corr_cond <- Preds_RatingDist_corr_cond %>% 
    group_by(model, experiment, condition, correct) %>%
    summarise(MRating = sum(p*rating)/sum(p), .groups = "drop")
  
  ### f) Prediction aggregated RT densities                ----
  RT_dist  <- preds_RT %>%  group_by(model, experiment, participant, correct, rating, condition, rt) %>%
    summarise(dens = mean(dens), densscaled = mean(densscaled), .groups = "drop")
  #   Combine and Aggregate RT densities 
  Ns_part <- Data %>% group_by(participant) %>% 
    summarise(N=n(), MinRT = min(rt), .groups = "drop")  %>%
    select(participant, N)
  Preds_RTdens_corr_cond_rating <-  RT_dist %>% 
    left_join(Ns_part, by="participant") %>%
    ungroup() %>%
    select(-participant) %>%
    group_by(rating, condition, model, experiment, correct, rt) %>%
    summarise(dens = sum(dens*N)/nrow(Data), densscaled = sum(N*densscaled)/nrow(Data), .groups = "drop") %>%
    mutate(model = factor(model, levels = c("dynaViTE", "2DSDT", "dynWEV", "2DSD"),
                          labels = c("dynaViTE", "2DSDT", "dynWEV", "2DSD")), 
           condition = factor(condition+ifelse(experiment=="Hellmann et al. (2023)\nExperiment 2", 5, 0)+
                                ifelse(experiment=="Shekhar & Rahnev (2021)\nExperiment 4", 10, 0),
                              levels = 1:13, 
                              labels = levels(Data$condition)))
  ### g) Prediction response time quantiles                ----
  Preds_RTQuants_corr_rating <- Preds_RTdens_corr_cond_rating %>% 
    group_by(model, experiment, rt, correct, rating) %>%
    summarise(dens = mean(dens), .groups = "drop") %>%
    group_by(model, experiment, correct, rating) %>% 
    do(PDFtoQuantiles(.)) %>%  
    ungroup() %>%
    left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
                                 model, experiment, rating, correct), 
                        p_correct=mean(p), .groups = "drop"), 
              by = c("model", "experiment", "correct", "rating"))
  # Preds_RTQuants_corr_cond <-   Preds_RTdens_corr_cond_rating %>%
  #   group_by(model, experiment, rt, correct, condition) %>%
  #   summarise(dens = sum(dens), .groups = "keep") %>%
  #   group_by(model, experiment, correct, condition) %>% 
  #   do(PDFtoQuantiles(.)) %>%  
  #   left_join(summarise(group_by(Preds_RatingDist_corr_cond, 
  #                                model, experiment, condition, correct), 
  #                       p_correct=sum(p), .groups="drop"),
  #             by = c("model", "experiment", "correct", "condition"))
  ### h) Remove raw prediction df's from environment       ----
  rm(list=c("preds_RT", "RT_dist", "preds_confidence"))
  ### i) Save aggregated data and predictions              ----
  save(Data,
       Data_RatingDist_part, Data_RatingDist_corr_cond,
       Data_MRating_corr_cond, Data_RTQuants_corr_rating, 
       Data_RTQuants_corr_cond,
       fits, 
       Preds_RatingDist_corr_cond, Preds_MRating_corr_cond ,
       Preds_RTdens_corr_cond_rating,Preds_RTQuants_corr_rating, 
       file="collected_Data_Fits_Predicts.RData")
} else {
  load("collected_Data_Fits_Predicts.RData")
}

#__________________________________________________________----
# B  Generate visualizations                               ----
## Include font style and define colors for plots          ----
#windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")

## Figure 6: Plot of Mean Ratings                          ----
pd <- position_dodge(0.02)
Preds_MRating_corr_cond_plot <- mutate(Preds_MRating_corr_cond,
                                       confidence ="Confidence measure") 
Data_MRating_corr_cond_plot <- merge(Data_MRating_corr_cond, 
                                     expand.grid(model=c("dynaViTE", "dynWEV", "2DSDT", "2DSD"), 
                                                 confidence ="Confidence measure"))


p_MRating <- ggplot(Data_MRating_corr_cond_plot,
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  # geom_ribbon(data=Preds_MRating_corr_cond, 
  #             aes(ymin=MRating-sewithin, ymax=MRating+sewithin, fill=as.factor(correct)), alpha=0.5)+ #
  geom_line(data=Preds_MRating_corr_cond_plot, aes(color=as.factor(correct)), linewidth=0.8)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.25,position =pd, linewidth=0.6) +
  geom_point(position = pd, fill="white", size=1.8)+
  facet_nested(rows=vars(model), cols=vars(experiment), scales = "free_x", drop=TRUE)+ #, dir="v"
  scale_x_discrete(name="Stimulus-onset-asynchrony [ms]                  Motion coherence [%]                                     Contrast [%]            ")+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_fill_manual(values= two_colors_correct, breaks=c(1,0),
                    name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed", labels=c("Correct", "Wrong"))  +
  scale_y_continuous(limits = c(1, 5), breaks=1:5, labels = c("0", "25", "50", "75", "100"), name="Mean confidence")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3), fill="none")+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(angle=90),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom", 
        legend.margin = margin(0,0.3,0,0.3,"cm"),
        legend.direction = "horizontal", #legend.position = "bottom"
        legend.box="horizontal",
        legend.key = element_blank(), legend.spacing = unit(0,"line"),
        legend.key.width=unit(1,"line"),
        panel.spacing=unit(0, "lines"))
p_MRating
dir.create("figures", showWarnings = FALSE)
ggsave("figures/meanRating1.tiff",
       width = 17.62, height=9/0.7, units="cm",dpi=600)
ggsave("figures/meanRating1.eps",
       width = 17.62, height=14, units="cm",dpi=1200, device = cairo_ps)

## Figure 7: RTQuantiles accross correct X rating          ----
Data_RTQuants_corr_rating_plot <- Data_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))%>% 
  filter(p %in% c(0.1, 0.5, 0.9))
Preds_RTQuants_corr_rating_plot <- Preds_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))),
         model = factor(model, levels=c("dynaViTE", "dynWEV", "2DSDT", "2DSD")))%>% 
  filter(p %in% c(0.1, 0.5, 0.9))
ggplot()+
  geom_line(data=mutate(Preds_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                        rating = as.factor(rating)),
            aes(x=rating, y=log(q), group=as.factor(p),color=correct), linewidth=0.7)+
  geom_point(data=mutate(Data_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                         rating = as.factor(rating)),
             aes(x=rating, y=log(q), shape=correct),
             size=1.2, fill="white")+
  scale_color_manual(values= two_colors_correct, breaks=c("Correct", "Wrong"),
                     name = "Predicted",
                     labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c("Correct", "Wrong"),
                     name = "Observed",
                     labels=c("Correct", "Wrong"))  +
  scale_x_discrete(name="Confidence [% scale width]", breaks=1:5,
                   labels=c("0-20","20-40", "40-60","60-80", "80-100"))+
  scale_y_continuous(breaks = log(c(1, 2, 4, 7)),
                     # labels = paste("log(", c(1.5, 2.5, 4, 7), ")", sep=""), 
                     labels = paste(c(1, 2, 4, 7), sep=""), 
                     name="Reaction time quantiles [s] (log scaled)")+
  facet_nested(model ~experiment+correct)+ #,dir="v"
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        legend.text = element_text(size=9),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.text.x = element_text(size=7, family="Times", color="black", angle=90),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9, family = "Times", angle=0),
        strip.text.y.right = element_text(angle = -90),
        legend.box = "horizontal",
        legend.position = "bottom",# legend.position = c(.70, .005), legend.justification = c(0,0),
        legend.direction = "horizontal", legend.spacing.y = unit(0, "lines"),
        legend.margin =margin(0,0,0,1, "cm"), legend.box.spacing = unit(0.2,"lines"),
        #legend.title = element_blank(),
        legend.key = element_blank(),
        legend.key.width=unit(1.5,"line"),
        panel.spacing=unit(0, "lines"))

ggsave("figures/RTQuantsConf1.eps",
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ggsave("figures/RTQuantsConf.tiff",
       width = 15, height=13, units="cm",dpi=600)   # Filling a whole power point slide
# ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
#        width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
# ggsave("../../Draft/figures/results/RTQuantsConf1.eps",
#        width=18.288, height=15.24, dpi=1200, units="cm", device=cairo_ps)


#__________________________________________________________----
# C  Quantitative comparison                               ----
compute_bf <- function(x,y) {
  drop <- is.na(x) | is.na(y)
  x <- x[!drop]
  y <- y[!drop]
  bf <- ttestBF(x,y,  paired = TRUE, rscale=1)
  bfdir <- ttestBF(x,y,  paired = TRUE, rscale=1, nullInterval = c(-Inf, 0))
  res <-  data.frame(bf = as.numeric(extractBF(bfdir)["bf"][2,1]), Lower = NA, Upper = NA)
  if (! is.na(res$bf)) {
    posterior <- quantile(posterior(bf, iterations=100000)[,"delta"], probs = c(.025, .975))
    res$Lower <- posterior[1]
    res$Upper <- posterior[2]
  }
  res
}
## Table 2: dynaViTE vs. dynWEV and 2DSDT vs. 2DSD         ----
BIC_BF_time_dep <-  fits %>%
  mutate(time_dep = ifelse(model %in% c("dynaViTE", "2DSDT"), "general", "restricted"),
         model = ifelse(model %in% c("dynaViTE", "dynWEV"), "dynaViTE", "2DSDT")) %>%
  select(c("participant", "model","time_dep",  "BIC", "experiment")) %>%
  pivot_wider(id_cols=c("participant", "experiment", "model"), names_from = time_dep, values_from = BIC) %>%
  group_by(experiment, model) %>%
  #filter(!(model=="dynWEV" & participant == 207)) %>%
  summarise(MeanBIC = mean(restricted-general),
            N = n(),
            SDBIC = sd(restricted-general),
            SERBIC = sd(restricted-general)/sqrt(n()),
            compute_bf(restricted, general))%>% ungroup() %>%
  mutate(model= factor(model)) %>% 
  arrange(experiment, desc(model))
BIC_BF_time_dep

## Figure 8: Comparison dynaViTE vs. all other models      ----
descr_BIC <- fits[,c("participant", "model", "BIC", "experiment")] %>% 
  group_by(experiment) %>%
  reframe(Rmisc::summarySEwithin(pick(everything()), measurevar="BIC", withinvars = "model", idvar="participant")) %>%
  mutate(model = factor(model, levels=rev(c("dynaViTE", "dynWEV", "2DSDT", "2DSD"))))
plot_BF <- fits[,c("participant", "model", "BIC", "experiment")] %>%
  group_by(model,participant, experiment) %>%  
  filter(model!="dynaViTE") %>% 
  left_join(fits[fits$model=="dynaViTE",c("participant","BIC", "experiment")],
            by=c("participant", "experiment")) %>%
  group_by(model, experiment) %>%
  summarise(compute_bf(BIC.x, BIC.y))%>% ungroup() %>%
  mutate(model= factor(model)) %>%
  mutate(BF_power = floor(log(bf, base=10)),
         BF_coeff = format(round(bf/10^BF_power, digits = 1),nsmall=1),
         BF_label=ifelse(BF_power > 1, paste0("BF==",format(BF_coeff, nsmall=1), "*x*10^~", BF_power),
                         ifelse(BF_power >= 0, paste0("BF==",format(round(bf,digits=2), nsmall=2)),
                                paste0("BF==",format(round(bf, digits=3), nsmall=3)))))%>%
  left_join(summarise(group_by(descr_BIC, experiment), 
                      minBIC = min(BIC),
                      rangediffBIC = diff(range(BIC)),
                      .groups="drop")) %>%
  select(experiment, model, BF_label, rangediffBIC, minBIC) 
my_breaks <- function(x) { 
  if (max(x) > - 3990) {
    c(-4060, -4020, -3980)
  } else if (max(x) > -5600) {
    c(-5800, -5600, -5400)
  } else {
    c(-11360, -11330, -11300)
  }
}

p_BIC <- ggplot(descr_BIC, aes(x=model, y=-BIC, group=experiment))+
  geom_line()+ #stat="identity", 
  geom_errorbar(aes(ymin=-BIC-se, ymax=-BIC+se), width=0.2)+
  facet_nested(.~experiment, scales="free_x", axes="x", independent = "x")+
  #facet_nested(experiment~., scales="free_y", axes="y", independent = "y")+
  geom_text(data=plot_BF, aes(label=BF_label, y=-minBIC), parse=TRUE,
            hjust=0, vjust=0.5, size=12/.pt)+
  geom_point(data=plot_BF, aes(y=-minBIC+0.8*rangediffBIC), alpha=0)+
  scale_y_continuous(breaks = my_breaks)+
  coord_flip()+
  xlab("")+ylab("-BIC")+
  theme_bw() +
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=14, family="Times"),
        axis.text = element_text(size=14, family="Times", color="black"),
        axis.text.y = element_text(hjust = 1), #( angle=90, hjust = 0.5),
        # panel.grid.minor = element_blank(),  # switch off minor gridlines
        # panel.grid.major = element_blank(),
        strip.text = element_text(size=14),
        legend.text = element_text(
          margin = margin(l = 4, r=4, t=2, b=2,unit = "pt")))
ggsave("figures/BIClines_againstdynavite.eps", plot=p_BIC,
       width=17.62, height=7, dpi=1200, units="cm", device=cairo_ps)
ggsave("figures/BIClines_againstdynavite.jpg", plot=p_BIC,
       width = 22, height=7, units="cm",dpi=1200)

#__________________________________________________________----
# D  dynaViTE and dynWEV parameter analysis                ----
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
## Table 3: Mean proportion of tau, t0, and decision time  ----
mean_posttime <- data.frame()
for (i in 1:4 ) {
  curmodel <- unique(fits$model)[i]
  temp <- Data %>% left_join(select(subset(fits, model==curmodel), participant, tau, t0)) %>%
    group_by(experiment) %>%
    reframe(prop_post = round(mean(tau/rt),2),
            prop_t0 = round(mean(t0/rt),2),
            prop_dec = 1-prop_post-prop_t0)  
  mean_posttime <- rbind(mean_posttime, cbind(model=curmodel, temp))
  #if (i==1)  mean_posttime[i,] <- c(curmodel, round(temp$prop_post,2))
}
mean_posttime <- mean_posttime[,c(2,1,3,4,5)] %>%
  mutate(model = factor(model, levels=c("dynaViTE", "dynWEV", "2DSDT", "2DSD")))  %>%
  arrange(experiment, model) %>%
  mutate(experiment=ifelse(model=="dynaViTE", experiment, " ")) %>%
  filter(!grepl("2DSD", model)) %>%
  mutate(experiment = str_replace(experiment, "\n", " "))
write.table(mean_posttime, file="figures/data_accuracies.txt",
            sep = ",", quote = FALSE, row.names = F)

## Figure 9: Trade-off between tau and t0                  ----

diff_dynaViTE_vs_dynWEV_in_tau_and_t0 <- parfits_long_dynaViTEdynWEV %>% 
  filter(parameter %in% c("tau", "t[0]")) %>%
  mutate(diff = dynWEV-dynaViTE) %>% 
  pivot_wider(names_from = "parameter", values_from = "diff", 
              id_cols = c("participant", "experiment"))%>%
  mutate(experiment = str_replace(experiment, "\n", " "))
p_taut0 <- ggplot(diff_dynaViTE_vs_dynWEV_in_tau_and_t0, 
       aes(x=tau, y=`t[0]`, color=experiment))+
  geom_abline(slope= -1, col="black", alpha=0.7)+
  geom_point()+ 
  theme_bw()+
  scale_color_manual(name="",values = experiment_colors )+
  scale_y_continuous(name=parse(text="Delta~t[0]"))+xlab(parse(text="Delta~tau"))+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = "bottom", legend.direction = "vertical",
        legend.margin = margin(-0.5, 0, 0, 0, "cm"),
        #legend.position = c(0.99, 0.99), legend.justification = c(1, 1),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9))
p_taut0
ggsave("figures/differencetaut0.tiff", plot=p_taut0,
       width = 6.49, height=6.5, units="cm",dpi=600)
ggsave("figures/differencetaut0.eps", plot=p_taut0,
       width = 6.49, height=6.5, units="cm",dpi=600, device = cairo_ps)


## Figure 10: Cor DDM parameters btw dynaViTE and dynWEV   ----
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
ggplot(filter(DDMparfits_long_dynaViTEdynWEV, abs(dynWEV)<1e+6 & abs(dynaViTE) < 1e+6), 
       aes(x=dynWEV, y=dynaViTE, color=experiment))+
  geom_smooth(data=filter(DDMparfits_long_dynaViTEdynWEV, abs(dynWEV)<1e+6 & abs(dynaViTE) < 1e+6),  
              method="lm", alpha=0.5, linewidth=0.9, 
              aes(x=dynWEV, y=dynaViTE), formula = 'y~x', inherit.aes = FALSE)+
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
  theme_bw()+ylab("Fitted parameter for dynaViTE")+xlab("Fitted parameter for dynWEV")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = c(0.78, 0.28), legend.justification = c(0, 1),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9))
ggsave("figures/parametersDynWEVdynaViTE.tiff",
       width = 17.625, height=12, units="cm",dpi=600)
ggsave("figures/parametersDynWEVdynaViTE.eps",
       width = 17.625, height=12, units="cm",dpi=600, device = cairo_ps)

#__________________________________________________________----
# E  Parameter and model recovery                          ----
if (!file.exists("results_mimikry_analysis.RData")) {
  #### Simulate and fit                                    ####
  part_Ns <- bind_rows(Data_KonfMask[,c("participant")], Data_RDK[,c("participant")], Data_Shekhar[,c("participant")]) %>% 
    group_by(participant) %>%
    summarise(Ntrials = n())
  generative_models_for_comparison <- c("dynaViTE", "dynWEV")#, "2DSDT", "2DSD")
  
  # Bootstrap sample from fitted parameter sets
  Nsims = 50
  set.seed(50)
  planned_analysis <- expand.grid(n_simu = 1:Nsims,
                                  model=generative_models_for_comparison) %>%
    mutate(participant = sample(fits$participant, replace=TRUE,
                                size=length(generative_models_for_comparison)* Nsims)) %>%
    left_join(part_Ns)
  # Define a function to be parallelized over for simulation and fitting
  fun_sim_fit <- function(simu_row) {
    # Get information for currect simulation
    gen_model <- simu_row$model
    cur_participant <- as.numeric(simu_row$participant)
    cur_trials <- as.numeric(simu_row$Ntrials)
    cur_n_simu <- as.numeric(simu_row$n_simu)
    paramDf <- subset(fits, participant == cur_participant & model==gen_model)    
    paramDf <- paramDf[ , colSums(is.na(paramDf))==0]
    nConds <-  length(grep(pattern = "^v[0-9]", names(paramDf), value = T))
    cur_participant <- as.numeric(
      paste0("1", str_pad(cur_n_simu, 2, side="left", "0"),
             str_pad(cur_participant, 3, side="left", "0")))
    set.seed(cur_participant+cur_trials+100*cur_n_simu+20*str_length(gen_model))
    print(paste("Starting participant:", cur_participant, "model:", gen_model, "simulation:", cur_n_simu))
    
    outnames <- colnames(fits)
    out <- data.frame(matrix(NA, nrow=0, ncol=length(outnames)))
    colnames(out) <- outnames
    
    # Generate articifial data 
    sim_df <- simulateRTConf(paramDf = paramDf, model = gen_model,
                             n = ceiling(cur_trials/(2*nConds)),  # number of trials per cond and stimulus  
                             delta=0.01,        # discretization step for simulation (in sec.) 
                             maxrt=30, simult_conf = TRUE)         # maximum decision time for each trial
    sim_df <- sim_df %>% filter(response !=0) # if maxrt is exceeded; response is set to 0
    
    sim_df$participant <- cur_participant
    save(sim_df, file=paste0("simulated_data_part", cur_participant, "_model", gen_model,"_simulation", cur_n_simu,".RData"))
    if (grepl("dyn", gen_model)) {
      models <- c("dynaViTE", "dynWEV")
    } else {
      models <- c("2DSDT", "2DSD")
    }
    # Fit models and combine output
    res_fit_gen_model <- fitRTConfModels(sim_df, models=models, nRatings=5, 
                                         restr_tau="simult_conf",
                                         logging=TRUE,
                                         opts=list(nAttempts=4),
                                         parallel = FALSE)
    true_negLogLik <- -LogLikWEV(sim_df, paramDf = paramDf, model = gen_model, simult_conf = TRUE)
    
    out[1:2, colnames(res_fit_gen_model)] <- res_fit_gen_model
    out$gen_model <- gen_model 
    out$true_negLogLik <- true_negLogLik
    
    print(paste("Finished participant:", cur_participant, "model:", gen_model, "simulation:", cur_n_simu))
    return(out)
  }
  # Change working directory for autosaving results
  dir.create("identification_analysis")
  setwd("identification_analysis")
  t00 <- Sys.time()
  
  planned_analysis <- split(planned_analysis, seq(nrow(planned_analysis)))
  n.cores <- parallel::detectCores() - 1
  cl <- makeCluster(n.cores, "SOCK", outfile = "")
  registerDoSNOW(cl)
  clusterEvalQ(cl, library(dynConfiR))
  clusterEvalQ(cl, library(tidyverse))
  clusterExport(cl, c("planned_analysis", "fun_sim_fit", "fits"))
  mimikry_collected_results <- parLapplyLB(cl, planned_analysis,fun_sim_fit)
  stopCluster(cl)
  
  mimikry_collected_results <- do.call(rbind, mimikry_collected_results)
  mimikry_collected_results <- 
    mimikry_collected_results %>%
    left_join(distinct(select(fits, c("participant", "experiment"))))
  setwd("..")
  save(file = "results_mimikry_analysis.RData", mimikry_collected_results)
  
  print("Finished model mimikry analysis!")
  print(paste("Mimikry analysis took...",
              as.character(round(as.double(difftime(Sys.time(),t00,units = "mins")), 2)),
              " mins"))
} else {
 #### If recovery analysis was done before, load results   ####
  load("results_mimikry_analysis.RData")
}
## Analyse recovery results                                ----
# Identification accuracy for all simulations 
# (also thos with lambda=0)
mimikry_collected_results %>%
  group_by(simu, gen_model) %>%
  summarise(BIC = model[which(BIC==min(BIC))],
            AIC = model[which(AIC==min(AIC))],
            AICc = model[which(AICc==min(AICc))]) %>%
  group_by(gen_model) %>%
  reframe(accBIC = mean(BIC==gen_model),
          accAIC = mean(AIC==gen_model),
          accAICc = mean(AICc==gen_model))

# Fitted lambda for dynWEV-generated data
mimikry_collected_results %>%
  filter(gen_model=="dynWEV" & model=="dynaViTE") %>%
  summarise(meanlambda=round(mean(lambda),2),
            sdlambda=round(sd(lambda),2)) 

# Get true parameters
true_pars <- mimikry_collected_results %>%
  select(participant, experiment, simu, model=gen_model) %>%
  distinct() %>%
  filter(model=="dynaViTE") %>%
  left_join(fits) 
# Number of correctly or incorrectly identified simulations
# generated by dynaViTE depending on the true generating lambda
mimikry_collected_results %>% filter(gen_model =="dynaViTE") %>%
  group_by(simu,  participant, experiment,N) %>%
  summarise(identified =model[which(BIC==min(BIC))], .groups = "drop") %>%
  left_join(select(true_pars, true_lambda=lambda,simu, participant, experiment)) %>% 
  group_by(true_lambda==0, identified) %>%
  summarise(N = n())

## Figure 12: Proportions of correctly identified models  ----
props_identified <- mimikry_collected_results %>%
  left_join(select(true_pars, true_lambda=lambda,simu, participant, experiment)) %>%
  filter(gen_model=="dynWEV" | true_lambda!=0) %>%
  pivot_longer(cols=c("BIC", "AIC", "AICc"), values_to = "IC", names_to = "criterion") %>%
  group_by(simu, gen_model, criterion) %>%
  summarise(best_model = model[which(IC==min(IC))]) %>%
  group_by(gen_model, criterion, best_model) %>%
  reframe(N=n()) %>% arrange(gen_model, criterion)%>% 
  group_by(gen_model, criterion) %>%
  mutate(prop = N/sum(N)) %>% ungroup() 
props_correct_id <-   
  props_identified %>% mutate(correct = gen_model==best_model) %>%
  group_by(criterion) %>% mutate(Ntot=sum(N)) %>%
  group_by(criterion, correct) %>%
  reframe(prop=sum(N)/Ntot[1])
pie_theme <-    theme_void()+  
  theme(
    #plot.background =  element_rect(fill="white", color = "white"),
    plot.margin = margin(0, 0, 0, 0, "cm"),
    text = element_text(size=9, family="Times"),
    legend.text = element_text(size=9),
    strip.text = element_text(size=9, family = "Times", angle=0),
    strip.text.y.right = element_text(angle = -90),
    legend.direction = "vertical", legend.position = "bottom", 
    legend.box = "horizontal", legend.spacing.y = unit(0.5, "lines"),
    legend.margin =margin(0,0,0,0, "cm"), legend.box.spacing = unit(0.2,"lines"),
    #legend.title = element_blank(),
    legend.key = element_blank(),
    legend.key.width=unit(1.5,"line"),
    panel.spacing=unit(0, "lines"))

strips <- strip_nested(text_y = list(element_blank()),
                       background_y = list(element_blank()), by_layer_y = TRUE)
gen_model_identification <- ggplot(props_identified , aes(x="", y=prop, fill=best_model))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  facet_nested(rows=vars(criterion), 
               cols=c(vars("Generative model"),vars(gen_model)),
               strip = strips)+
  geom_text(data=subset(props_identified , gen_model==best_model),
            aes(y = 0.5*prop, label = format(round(prop, 2), nsmall=2)), 
            color = "black", size=9/.pt) +
  scale_fill_brewer(palette="Set2", name="Model favored by IC")+
  pie_theme
total_identification <- ggplot(props_correct_id , aes(x="", y=prop, fill=as.factor(correct)))+
  geom_bar(stat="identity", width=1, color="white")+
  coord_polar("y", start=0)+
  facet_nested(rows=c(vars(criterion)), 
               cols=c(vars("Total"), vars("")))+
  geom_text(data=subset(props_correct_id, correct),
            aes(y = 0.5*prop, label = format(round(prop, 2), nsmall=2)), 
            color = "black", size=9/.pt) +
  scale_fill_brewer(palette="Set1", name="Identified",
                    labels=c("Incorrect", "Correct"))+
  pie_theme
ggpubr::ggarrange(gen_model_identification, total_identification, 
                  nrow=1, widths = c(.64, .36)) + 
  ggpubr::bgcolor("white")+
  ggpubr::border("white")
ggsave("figures/Identification_Accuracy.tiff",
       width = 8, height=9, units="cm",dpi=600)   # Filling a whole power point slide


## Figure 12: True and recovered dynaViTE parameters       ----
experiment_colors = c("#AA4499","#117733", "#DDCC77")
par_labels <- c('nu[1]', 'nu[2]', 'nu[3]', 'nu[4]', 'nu[5]', 's[nu]', 'a', 'z', 's[z]', 
                't[0]', 's[t0]', 'tau', 'w', 'sigma[V]', 's[V]', 'lambda',
                paste0('theta[', -1,'~', 1:4,']'),#'',
                paste0('theta[', 1,'~',  1:4,']'))
par_levels <- c('v1', 'v2', 'v3', 'v4', 'v5', 'sv', 'a', 'z', 'sz',
                't0', 'st0', 'tau', 'w', 'sigvis', 'svis','lambda',
                paste0('thetaLower',1:4),
                paste0('thetaUpper',1:4))
comp_true_rec_dynaViTE <- mimikry_collected_results %>% 
  filter(gen_model =="dynaViTE" & model=="dynaViTE") %>%
  select(par_levels, participant, simu, experiment) %>%
  mutate(true_rec="recovered") %>%
  rbind(mutate(select(true_pars,simu, participant, experiment,
                      par_levels), true_rec="true"))
comp_true_rec_dynaViTE <- comp_true_rec_dynaViTE %>% 
  pivot_longer(cols = par_levels, names_to ="parameter") %>%
  pivot_wider(names_from = true_rec, values_from=value) %>%
  mutate(parameter = factor(parameter, levels = par_levels,labels = par_labels))

CCC <- function(x, y) { ## Function to compute the 
  # concordance correlation coefficient 
  NAinds <- is.na(x) | is.na(y)
  x <- x[!NAinds]; y <- y[!NAinds]
  k <- length(y)
  sy2 <- var(y) * (k - 1)/k
  sx2 <- var(x) * (k - 1)/k
  est <- 2 * cov(x, y) *(k-1)/k /(sx2 + sy2 + (mean(y) - mean(x))^2)
  return(est)
}
est_cors <- subset(comp_true_rec_dynaViTE, 
                   abs(true) < 1e+6 & abs(recovered) < 1e+6) %>%
  group_by(parameter) %>%
  reframe(est = CCC(true, recovered),
          x_min=max(true), 
          y_max=min(recovered))
comp_true_rec_dynaViTE <- comp_true_rec_dynaViTE %>%
  mutate(experiment = str_replace(experiment, "\n", " "))
ggplot(filter(comp_true_rec_dynaViTE, abs(true)<1e+6 & abs(recovered) < 1e+6), 
       aes(x=true, y=recovered, color=experiment))+
  geom_smooth(data=filter(comp_true_rec_dynaViTE, abs(true)<1e+6 & abs(recovered) < 1e+6),  
              method="lm", alpha=0.5, linewidth=0.9, 
              aes(x=true, y=recovered), formula = 'y~x', inherit.aes = FALSE)+
  facet_wrap(.~parameter, scales = "free", labeller = label_parsed, 
             drop = TRUE, ncol=4)+
  geom_abline(col="black", linewidth=1, alpha=0.2)+
  scale_color_manual(name="", values=experiment_colors)+
  geom_point(alpha=0.5)+ scale_shape_discrete("")+
  geom_text(data = est_cors, #mutate(est_cors, experiment="Experiment 1"),
            aes(x=x_min, y=y_max,
                label=paste0(format(round(est, 2), nsmall=2))),#"rho == ",
            color="black",
            parse=TRUE, hjust=1, vjust=0)+
  theme_bw()+ylab("Recovered parameter")+xlab("True generating parameter")+
  theme(plot.margin = margin(0, 0.5, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.position = "bottom",# legend.justification = c(0, 1),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        strip.text = element_text(size=9))
ggsave("figures/parameterstruerecovered.tiff",
       width = 17.625, height=21, units="cm",dpi=600)
ggsave("figures/parameterstruerecovered.eps",
       width = 17.625, height=21, units="cm",dpi=600, device = cairo_ps)


#__________________________________________________________----
# F  Discussion and Supplement                             ----
## Table 4: Accuracy and w parameter                       ----
fitted_w <- fits %>% filter(model%in% c("dynWEV","dynaViTE")) %>%
  group_by(model,  experiment) %>%
  summarise(Meanw=mean(w), 
            SDw=sd(w)) %>%
  mutate(w = paste(format(round(Meanw, 2), nsmall=2), 
                   "(", format(round(SDw, 2), nsmall=2), ")", sep="")) %>%
  select(-Meanw, -SDw) %>%
  pivot_wider(id_cols=experiment, values_from = w, names_from = model)
Data_Acc <- Data %>% group_by(participant, experiment, condition) %>%
  summarise(Acc = mean(correct), .groups="drop") %>%
  group_by(experiment) %>%
  reframe(summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("condition"))) %>%
  rename(sewithin = se) 
Accuracies <- Data_Acc %>% 
  select(experiment, condition, Acc, sewithin) %>%
  mutate(Acc=paste(format(round(Acc, 2), nsmall=2),"(", format(round(sewithin, 2), nsmall=2), ")", sep="")) %>%
  select(-sewithin) %>%
  mutate(condition=as.character(condition)) %>%
  rbind(data.frame(experiment=rep("Shekhar & Rahnev (2021)\nExperiment 4",2), condition=rep("--", 2), Acc=rep("--", 2))) %>%
  pivot_longer(cols = c("condition", "Acc"), values_to="value") %>%
  mutate(names2=rep(rep(1:5,each=2), 3)) %>%
  pivot_wider(id_cols=c("experiment", "name"), values_from = "value", names_from = names2)
Accuracies <- Accuracies %>% right_join(fitted_w)
write.table(Accuracies, file="figures/data_accuracies.txt",
            sep = ",", quote = FALSE, row.names = F)
## Suppl Fig 6: Plot of Fitted Accuracy                    ----
Data_Acc <- Data %>% group_by(participant, experiment, condition) %>%
  summarise(Acc = mean(correct), .groups="drop") %>%
  group_by(experiment) %>%
  reframe(summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("condition"))) %>%
  rename(sewithin = se) 

Preds_Acc <- Preds_RatingDist_corr_cond %>% 
  group_by(model, experiment, condition) %>%
  summarise(Acc = sum(p*correct)/sum(p), .groups="drop") 
p_Acc <- ggplot(Data_Acc,
                aes(x=condition, y=Acc)) +
  geom_line(data=Preds_Acc, aes(x = condition, linetype="Predicted", group=model), linewidth=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=Acc-sewithin, ymax=Acc+sewithin), colour="black", width=.2, linewidth=0.4) +
  geom_point(fill="white", aes(shape="Observed"))+
  facet_grid(rows=vars(model), 
             cols=vars(experiment), scales = "free_x", drop=TRUE)+ #, dir="v"
  ylab("Mean accuracy")+
  scale_x_discrete(
    name="Stimulus-onset-asynchrony [ms]                  Motion coherence [%]                                     Contrast [%]            ")+
  scale_linetype_manual(name="", values=1) +
  scale_shape_manual(values=c(21),name = "")  +
  scale_y_continuous(breaks=c(0.5, 0.7, 0.9), name="Mean Accuracy")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        axis.text.x = element_text(angle=90,hjust = 0.5, vjust = 0.5),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom",
        # legend.position=c(.80, .005), legend.justification = c(0,0),
        legend.margin = margin(-0.1, 0, 0, 0, "cm"),#legend.direction = "horizontal",
        legend.box.margin = margin(-0.1, 0, 0, 0, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
p_Acc

dir.create("figures", showWarnings = FALSE)
ggsave("C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/modelAccuracy1.eps",
       width = 17.62, height=17, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/modelAccuracy1.tiff",
       width = 17.62, height=9/0.75, units="cm",dpi=600)

## Suppl Table 1: Descriptive parameter fits dynaViTE      ----
meanFits <- fits %>% 
  select(-c("participant", "k", "N", "AICc","BIC", "AIC", "negLogLik","fixed")) %>%
  group_by(experiment, model) %>% 
  summarise(across(.cols = everything(),
                   ~paste(round(mean(.x[abs(.x)<1e+24], na.rm=TRUE),2), 
                          " (", round(sd(.x[abs(.x) <1e+24], na.rm = TRUE), 2), ")", sep=""))) %>%
  select(c("experiment", "model", paste0("v", 1:5), "sv", "a", "z", "sz", paste0("thetaUpper", 1:4), paste0("thetaLower", 1:4),
           "t0", "st0", "tau", "w","sigvis", "svis", "lambda"))
names(meanFits)[3:7] <- paste0("$\\nu_", 1:5, "$")
names(meanFits)[8] <- paste0("$s\\nu$")
names(meanFits)[11] <- paste0("$s_z$")
names(meanFits)[12:15] <- paste0("$\\theta_{1,", 1:4, "}$")
names(meanFits)[16:19] <- paste0("$\\theta_{-1,", 1:4, "}$")
names(meanFits)[20:21] <- c("$t_0$", "$s_{t0}$")
names(meanFits)[22] <- paste0("$\\tau$")
names(meanFits)[24: 25] <- c("$\\sigma_V$", "$s_V$")
names(meanFits)[26] <- paste0("$\\lambda$")
names(meanFits)
meanFits[meanFits=="NA (NaN)"]
meanFits[meanFits=="NA (NA)"]
meanFits[meanFits=="NaN (NaN)"]
meanFits[meanFits=="NaN (NA)"] <- "--"
meanFits[meanFits$model %in% c("dynWEV", "2DSD"), "$\\lambda$"] <- "--"

# # Use this to produce a wide table
# meanFits <- column_to_rownames(meanFits, "model")
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)
# print(xtable(meanFits, align = c("l", "l", rep("c", ncol(meanFits)-1))), type="latex", file="../Supplement/figures/tablepars1.tex",sanitize.text.function=function(x){x}, include.rownames=FALSE)

# Use this to produce a long table with only fits for dynaViTE
meanFits_long <- meanFits %>% 
  # Report only dynaViTE model  
  subset(model=="dynaViTE") %>% 
  select(-model) %>%
  pivot_longer(cols = 2:(ncol(meanFits)-1), names_to = "Parameter") %>%
  pivot_wider(id_cols = "Parameter", names_from = c("experiment"), values_from = value) 
meanFits_long

# Use flextable for word output
if (!"flextable" %in% installed.packages()[,"Package"]) {
  install.packages("flextable")
} 
library(flextable)
Fits_wordtable <- flextable(meanFits_long) %>%
  align(align="left", j=1) %>%
  align(align="center", j=2:4)
save_as_docx("Parameter table"=Fits_wordtable,
             path="../Draft/tablepars.docx")
print(Fits_wordtable, preview="docx")

# Use xtable for latex output
# Use xtable for latex output
if (!"xtable" %in% installed.packages()[,"Package"]) {
  install.packages("xtable")
}
library(xtable)
# Use this to produce a long table with all model fits
meanFits_long <- meanFits %>% 
  # Report only dynaViTE model  
  pivot_longer(cols = 3:(ncol(meanFits)-1), names_to = "Parameter") %>%
  pivot_wider(id_cols = "Parameter", names_from = c("experiment", "model"), values_from = value) 
meanFits_long
Fits_textable <- xtable(meanFits_long, align = c("l", "l", rep("c", ncol(meanFits_long)-1)),
                        caption="\\raggedright Mean and standard deviation of parameter fits")
names(Fits_textable) <- 
  c("Model", str_split_i(names(Fits_textable)[2:ncol(Fits_textable)], "_", 2))
experimentlabels <- str_replace(str_replace(unique(meanFits$experiment),"\n", " "), "&", "\\\\&")
addtorow <- list()
addtorow$pos <- list(-1, -1, c(0:(nrow(Fits_textable)-1)))
addtorow$command <- c('\\toprule  & \\multicolumn{12}{c}{Experiment}\\\\ \\cmidrule{2-13} ',
                      paste0(" & ", 
                             paste0("\\multicolumn{4}{c}{",experimentlabels , "}", collapse = " & "),
                             " \\\\ \\cmidrule{2-5} \\cmidrule{6-9} \\cmidrule{10-13}"),
                      '[1mm]')
print(Fits_textable, type="latex", 
      sanitize.text.function=function(x){x}, 
      include.rownames=FALSE, 
      add.to.row=addtorow, 
      hline.after = c(0, nrow(Fits_textable)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
print(Fits_textable, type="latex", 
      file="C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/tablepars2.tex",
      sanitize.text.function=function(x){x}, 
      include.rownames=FALSE, 
      add.to.row=addtorow, 
      hline.after = c(0, nrow(Fits_textable)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
# ### If we want to report all model, use this code:   
# ### Use this to produce a long table
# meanFits_long <- meanFits %>% pivot_longer(cols = 3:ncol(meanFits), names_to = "Parameter") %>%
#   pivot_wider(id_cols = "Parameter", names_from = c("experiment", "model"), values_from = value) %>%
#   select(c(1, 5, 3, 4, 2, 9, 7, 8, 6, 13, 11, 12, 10))
# meanFits_long

## Suppl Fig 7-9: Fitted response distributions            ----
model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")
for (i in 1:3) {
  cur_experiment <- sort(unique(fits$experiment))[i]
  p_ratingdist_KonfMask <- ggplot(data=filter(Data_RatingDist_corr_cond, experiment==cur_experiment), aes(x=rating, y=p))+
    geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
    geom_point(data=filter(Preds_RatingDist_corr_cond, experiment==cur_experiment),
               aes(shape=model, fill=model), 
               position=position_dodge(1))+
    facet_nested(condition~experiment+correct,
                 labeller = labeller(correct=c("0"="Wrong", "1"="Correct"),
                                     condition=function(x) paste(x, "ms")))+
    scale_x_continuous(name = "Confidence [% scale width]", breaks = 1:5,
                       labels = c("0-20","20-40", "40-60", "60-80", "80-100")) +
    scale_y_continuous(name = "Probability", breaks = c(0, 0.4, 0.8))+
    scale_shape_manual(name = "Model prediction",values=c(21, 22,23,24)) +
    scale_fill_manual(name = "Model prediction",values=model_colors) +
    theme_bw() +
    theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
          text = element_text(size=9, family="Times"),
          axis.text = element_text(size=9, family="Times", color="black"),
          axis.text.x = element_text(angle=90),
          axis.title.y = element_text( angle=90),
          panel.grid.minor = element_blank(),  # switch off minor gridlines
          panel.grid.major = element_blank(),
          legend.key = element_blank(),
          legend.position = "bottom", legend.text = element_text(size=12, family = "Times"),
          strip.text.y = element_text(angle = 0, size=9, family = "Times"),
          panel.spacing=unit(0, "lines"))+
    coord_cartesian(xlim = c(0.5,5.5),ylim=c(0, 0.9),  clip = 'off')
  p_ratingdist_KonfMask+
    theme(plot.margin = margin(0, 0, 0, 0, "cm"),
          legend.box = "vertical", legend.spacing.y = unit(-.2, "cm"))
  ggsave(paste0("figures/RespDist", 
                str_replace_all(cur_experiment, c("\n"="", " "="", "\\("="", "\\)"="", "\\."="")),
                ".tiff"),
         width=17.3, height=18, units="cm",dpi=1200)
  ggsave(paste0("C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/RespDist", 
                str_replace_all(cur_experiment, c("\n"="", " "="", "\\("="", "\\)"="", "\\."="")),
                ".eps"),
         width=17.3, height=18, units="cm",dpi=1200, device = cairo_ps)
  
}

# The end                                                  ####
