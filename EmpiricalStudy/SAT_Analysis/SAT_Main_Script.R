############################################################-
#####             Main script for the                       -
#####       speed-accuracy-tradeoff analyses in             -
#####       Time in dynamical confidence models             -
#__________________________________________________________----

# Sebastian Hellmann, 23.02.2024

## Structure:
# Preamble and imports                                   
# A   Read data and fit models 
## 1. Read experimental data                             
## 2. Fit the four models to the experimental data 
## 3. Simulate predicted distributions from model fits    
## 4. Aggregate data and predictions                     
# B  Generate visualizations                             
# C  Quantitative comparison (--> see main script)
# D  Discussion and Supplement                           


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
# A  Read data and fit models                              ----
## 1. Read and prepare experimental data                   ----
if (!file.exists("Preprocessed_SATdata.RDATA")) {
  SATdata <- read_csv2(file = "LineStudyData.CSV")
  names(SATdata)
  SATdata <- SATdata[,c(1, 3, 6, 7, 8, 10, 11)]
  names(SATdata) <- c("participant", "SAT", "condition", "correct",
                      "rt", "confidence", "RT2")
  
  SATdata <- SATdata %>% 
    mutate(
      rt = rt/1000, RT2 = RT2/1000,
      rating = as.numeric(as.factor(confidence)))
  unique(SATdata$confidence)
  
  SATdata %>%group_by(participant) %>%
    reframe(maxrt = max(rt),
            minrt=min(rt))
  
  cutoffs <- SATdata %>%group_by(participant) %>%
    reframe(maxrt = mean(rt)+4*sd(rt),
            minrt=0.3)#median(rt)-1*sd(rt))
  ggplot(SATdata, aes(x=rt))+
    geom_density(bw=0.02)+
     geom_vline(data=cutoffs, aes(xintercept=maxrt), color="red")+
     geom_vline(data=cutoffs, aes(xintercept=minrt), color="red")+
    facet_wrap(.~participant, scales = "free")
  cutoffs <- SATdata %>%group_by(participant) %>%
    reframe(maxrt = mean(RT2)+4*sd(RT2),
            minrt=0.15)#median(rt)-1*sd(rt))
  ggplot(SATdata, aes(x=RT2, col=as.factor(rating)))+
    geom_density(bw=0.1)+
    geom_vline(data=cutoffs, aes(xintercept=maxrt), color="red")+
    geom_vline(data=cutoffs, aes(xintercept=minrt), color="red")+
    facet_wrap(.~participant, scales = "free")
  temp <- SATdata %>% group_by(participant, rating) %>% 
    reframe(N=n())
  temp <- SATdata %>% group_by(participant, rating, condition) %>% 
    reframe(N=n())
  Nrowges <- group_by(SATdata, participant) %>% reframe(nrow=n())
  SATdata <- SATdata %>%
    group_by(participant) %>%
    filter(rt < mean(rt)+4*sd(rt) & RT2 < mean(RT2)+4*sd(RT2) & rt >.3 & RT2 > .15) %>%
    ungroup()
  Nrowfiltered <- group_by(SATdata, participant) %>% reframe(nrow=n())
  1-Nrowfiltered$nrow/Nrowges$nrow
  SATdata %>%group_by(participant) %>%
    reframe(maxrt = max(rt),
            minrt=min(rt),
            minRT2 = min(RT2), 
            maxRT2 = max(RT2))
  temp <- SATdata %>% group_by(participant, rating, correct) %>% 
    reframe(N=n())
  temp <- SATdata %>% group_by(participant, rating, condition) %>% 
    reframe(N=n())
  
  save(SATdata, file="Preprocessed_SATdata.RDATA")
} else {
  load("Preprocessed_SATdata.RDATA")
}

## 2. Fit the four models to the experimental data   ----
if (!file.exists("SATfits.RDATA")) {
  source("fitting_fcts/fitRTConf_SAT.R")
  source("fitting_fcts/fitRTConfModels_SAT.R")
  fits <- fitRTConfModelsSAT(SATdata, models=c("dynaViTE", "dynWEV", "2DSD", "2DSDT"), 
                             fixed = list(sym_thetas=TRUE, z=.5, st0=0),
                             logging=TRUE, parallel="models")#,opts=list(nAttempts=2, nRestarts=2, maxfun=500))
  save(fits, file="SATfits.RDATA")
} else {
  load("SATfits.RDATA")
}
## 3. Simulate predicted distributions from model fits      ----
SATdata <- SATdata %>% mutate(SATnum = as.numeric(as.factor(SAT)))
if (!file.exists("collected_Data_Fits_Predicts.RData")) {

  if (!file.exists("raw_simulations.RData")) {
    Rcpp::sourceCpp("fitting_fcts/RNG_WEV_SAT3.cpp")
    source("fitting_fcts/simulateWEV_SAT2.R")
    n <- 1e+5
    sim_combine <- function(paramDf) {
      RT2s <- subset(SATdata, participant==paramDf$participant)$RT2
      paramDf$tau0 <- min(paramDf$tau0, min(RT2s-1e-3))
      df <- subset(SATdata, participant==paramDf$participant)
      simu <- simulateWEV_SAT(paramDf, n, df, model=paramDf$model, maxrt = 30)
      return(simu)
      #return(cbind(model=paramDf$model, participant=paramDf$participant, simu))
    }
    set.seed(2345)
    simulations <- fits %>% group_by(model, participant) %>%
      do(sim_combine(.))
    #simulations <- simulations %>% mutate(SAT=ifelse(SAT==1, "Accuracy", "Speed"))
    
    save(simulations, file="raw_simulations.RData")
  } else {
    load("raw_simulations.RData")
  }
## 4. Aggregate data and predictions                       ----
  ### a) Data confidence rating distribution               ----
  Data <- SATdata %>% #filter(participant !=5) %>%
    group_by(SAT, participant, condition) %>%
    mutate(nrows=n(), 
           condition = factor(condition, levels=1:6, 
                              labels=c("32.27", "32.59", "33.23", "33.87", "34.51","35.15")), 
           SAT = as.factor(SAT))
  Data_RatingDist_part <- Data %>%# filter(participant !=5) %>%
    group_by(SAT, participant, condition, rating, correct) %>% 
    summarise(p = n()/(mean(nrows)),.groups = "drop") %>%
    full_join(y = merge(distinct(Data[,c("correct","rating", "SAT")]),
                        distinct(Data[,c("condition", "participant", "SAT")]), by="SAT"),
              by = c("SAT", "participant", "condition", "rating", "correct")) %>%
    mutate(p = ifelse(is.na(p), 0, p)) 

  # ### Sanity Checks:
  # sum(Data_RatingDist_part$p)
  # table((Data_RatingDist_part %>% group_by(SAT,condition, participant) %>% summarise(p = sum(p)))$p)
  # # 7*2*6 = 72
  
  # For the plots we won't differentiate between 
  # stimulus directions and participants
  Data_RatingDist_corr_cond <- Data_RatingDist_part %>% 
    group_by(SAT, correct, condition, rating) %>% 
    summarise(p = mean(p), .groups = "drop")
  #sum(Data_RatingDist_corr_cond$p)
  
  ### b) Data mean confidence rating                       ----
  Data_MRating_corr_cond <- Data %>% #filter(participant !=5) %>%
    group_by(participant, SAT, condition, correct) %>%
    summarise(MRating = mean(rating), .groups = "drop") 
  Data_MRating_corr_cond_SE <- Data_MRating_corr_cond %>%
    #group_by(SAT) %>%
    reframe(summarySEwithin(pick(everything()),measurevar = "MRating", 
                            idvar = "participant", withinvars = c("SAT", "condition", "correct"))) %>%
    select(SAT, condition, correct, sewithin = se) 
  Data_MRating_corr_cond <- Data_MRating_corr_cond %>%
    group_by(SAT, condition, correct) %>%
    summarise(MRating=mean(MRating, na.rm=TRUE), .groups="drop") %>%
    mutate(correct=as.factor(correct)) %>%
    left_join(Data_MRating_corr_cond_SE, by=c("condition", "correct", "SAT"))
  ### c) Data response time quantiles                      ----  
  # Reaction Time Quantiles of the Data grouped by rating and accuracy  
  Data_RTQuants_corr_rating <- Data %>%#filter(participant !=5) %>%
    group_by(SAT, rating, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  # Reaction Time Quantiles of the Data grouped by condition and accuracy  
  Data_RTQuants_corr_cond <- Data %>%#filter(participant !=5) %>%
    group_by(SAT, condition, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  
  
  ### d) Prediction confidence rating distribution         ----
  Preds <- simulations %>%# filter(participant !=5) %>%
    group_by(model, SAT, participant, condition) %>%
    mutate(nrows=n(), 
           condition = factor(condition, levels=1:6, 
                              labels=c("32.27", "32.59", "33.23", "33.87", "34.51","35.15")), 
           SAT = as.factor(SAT),
           model = factor(model, levels = c("dynaViTE", "2DSDT", "dynWEV", "2DSD"),
                                 labels = c("dynaViTE", "2DSD+", "dynWEV", "2DSD"))) %>%
    ungroup()
  rm(simulations)
  Preds_RatingDist_part <- Preds %>% 
    group_by(model, SAT, participant, condition, rating, correct) %>% 
    summarise(p = n()/(mean(nrows)),.groups = "drop") %>%
    full_join(y = merge(distinct(Preds[,c("model", "participant", "correct","rating")]),
                        distinct(Preds[,c("model", "condition", "participant", "SAT")]), by=c("model", "participant")),
              by = c("model", "SAT", "participant", "condition", "rating", "correct")) %>%
    mutate(p = ifelse(is.na(p), 0, p)) 
  
  ### Sanity Checks:
  # subset(Preds_RatingDist_part, is.na(rating))
  # sum(Preds_RatingDist_part$p)
  # table((Preds_RatingDist_part %>% group_by(model, SAT,condition, participant) %>% summarise(p = sum(p)))$p)
  
  # For the plots we won't differentiate between 
  # stimulus directions and participants
  Preds_RatingDist_corr_cond <- Preds_RatingDist_part %>% 
    group_by(model, SAT, correct, condition, rating) %>% 
    summarise(p = mean(p), .groups = "drop")
  #sum(Preds_RatingDist_corr_cond$p)
  
  ### e) Preds mean confidence rating                       ----
  Preds_MRating_corr_cond <- Preds %>% 
    group_by(model, participant, SAT, condition, correct) %>%
    summarise(MRating = mean(rating)) #mean(conf)) 
  Preds_MRating_corr_cond_SE <- Preds_MRating_corr_cond %>%
    group_by(model) %>%
    reframe(summarySEwithin(pick(everything()),measurevar = "MRating", 
                            idvar = "participant", withinvars = c("SAT", "condition", "correct"))) %>%
    select(model, SAT, condition, correct, sewithin = se) 
  Preds_MRating_corr_cond <- Preds_MRating_corr_cond %>%
    group_by(model, SAT, condition, correct) %>%
    summarise(MRating=mean(MRating, na.rm=TRUE), .groups="drop") %>%
    mutate(correct=as.factor(correct)) %>%
    left_join(Preds_MRating_corr_cond_SE, by=c("model", "condition", "correct", "SAT"))
  ### f) Preds response time quantiles                      ----  
  # Reaction Time Quantiles of the Preds grouped by rating and accuracy  
  Preds_RTQuants_corr_rating <- Preds %>%
    group_by(model, SAT, rating, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  # Reaction Time Quantiles of the Preds grouped by condition and accuracy  
  Preds_RTQuants_corr_cond <- Preds %>%
    group_by(model, SAT, condition, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  
  fits <- mutate(fits, model=factor(model, levels=rev(c("dynaViTE", "dynWEV", "2DSDT", "2DSD")),
                                    labels=rev(c("dynaViTE", "dynWEV", "2DSD+", "2DSD"))))
  
  ### g) Remove raw prediction df's from environment       ----
  rm(list=c("Preds"))
  gc()
  ### i) Save aggregated data and predictions              ----
  save(Data,
       Data_RatingDist_part, Data_RatingDist_corr_cond,
       Data_MRating_corr_cond, Data_RTQuants_corr_rating, 
       Data_RTQuants_corr_cond,
       fits, 
       Preds_RatingDist_part, Preds_RatingDist_corr_cond,
       Preds_MRating_corr_cond, Preds_RTQuants_corr_rating, 
       Preds_RTQuants_corr_cond,
       file="collected_Data_Fits_Predicts.RData")
} else {
  load("collected_Data_Fits_Predicts.RData")
  gc()
}

#__________________________________________________________----
# B  Generate visualizations                               ----
## Include font style and define colors for plots          ----
windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")

## Figure 5: Plot of Mean Ratings                          ----
pd <- position_dodge(0.02)
Preds_MRating_corr_cond_plot <- mutate(Preds_MRating_corr_cond,
                                       confidence ="Confidence measure") 
Data_MRating_corr_cond_plot <- merge(Data_MRating_corr_cond, 
                                     expand.grid(model=c("dynaViTE","dynWEV", "2DSD+", "2DSD"), 
                                                 confidence ="Confidence measure"))


p_MRating <- ggplot(Data_MRating_corr_cond_plot,
                    aes(x=condition, y=MRating, group = as.factor(correct), shape=as.factor(correct))) +
  geom_line(data=Preds_MRating_corr_cond_plot, aes(color=as.factor(correct)), linewidth=0.8)+
  geom_errorbar(aes(ymin=MRating-sewithin, ymax=MRating+sewithin), colour="black", width=.25,position =pd, linewidth=0.6) +
  geom_point(position = pd, fill="white", size=1.8)+
  facet_nested(rows=vars(model), cols=vars(SAT), scales = "free_x", drop=TRUE)+ #, dir="v"
  scale_x_discrete(name="Longer line length [mm]")+
  scale_color_manual(values= two_colors_correct, breaks=c(1,0),
                     name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_fill_manual(values= two_colors_correct, breaks=c(1,0),
                    name = "Predicted", labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c(1,0),
                     name = "Observed", labels=c("Correct", "Wrong"))  +
  scale_y_continuous(limits = c(1.5, 6.3), breaks=2:6, 
                     labels = seq(60, 100, by=10), name="Mean confidence")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3), fill="none")+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        #axis.text.x = element_text(angle=90),
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
ggsave("figures/meanRating1_SAT.tiff",
       width = 17.62*0.8, height=9/0.7, units="cm",dpi=600)
ggsave("figures/meanRating1_SAT.eps",
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
         model = factor(model, levels=c("dynaViTE", "dynWEV", "2DSD+", "2DSD")))%>% 
  filter(p %in% c(0.1, 0.5, 0.9)) #%>%  filter(model !="dynWEV")
ggplot()+
  geom_line(data=mutate(Preds_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                        rating = as.factor(rating)),
            aes(x=rating, y=log(q), group=as.factor(p),color=correct), linewidth=0.7)+
  geom_point(data=mutate(Data_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                         rating = as.factor(rating)),
             aes(x=rating, y=log(q), shape=correct),
             size=1.2, fill="white")+
  geom_line(data=mutate(Data_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                         rating = as.factor(rating)),
             aes(x=rating, y=log(q), group=as.factor(p)),
             linewidth=0.5)+
  scale_color_manual(values= two_colors_correct, breaks=c("Correct", "Wrong"),
                     name = "Predicted",
                     labels=c("Correct", "Wrong")) +
  scale_shape_manual(values=c(21,17),breaks=c("Correct", "Wrong"),
                     name = "Observed",
                     labels=c("Correct", "Wrong"))  +
  scale_x_discrete(name="Confidence", breaks=1:6, 
                   labels = seq(50, 100, by=10))+
  scale_y_continuous(breaks = log(c(0.5, 1, 2, 4)),
                     labels = paste(c(0.5, 1, 2, 4), sep=""), 
                     name="Reaction time quantiles [s] (log scaled)")+
  facet_nested(model ~SAT+correct)+ #,dir="v"
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

ggsave("figures/RTQuantsConf1_SAT.eps",
       width=17.62, height=15, dpi=1200, units="cm", device=cairo_ps)
ggsave("figures/RTQuantsConf_SAT.tiff",
       width = 15, height=15, units="cm",dpi=600)   # Filling a whole power point slide
# ggsave("../../Draft/figures/results/RTQuantsConf1_SAT.eps",
#        width=17.62, height=22, dpi=1200, units="cm", device=cairo_ps)
# ggsave("../../Draft/figures/results/RTQuantsConf1_SAT.eps",
#        width=18.288, height=15.24, dpi=1200, units="cm", device=cairo_ps)


#__________________________________________________________----
# C  Quantitative comparison                               ----
# --> See main script
#__________________________________________________________----
# D  Discussion and Supplement                             ----
## Table 6: Accuracy and w parameter                       ----
fitted_w <- fits %>% filter(model%in% c("dynWEV","dynaViTE")) %>%
  group_by(model) %>%
  summarise(Meanw=mean(w), 
            SDw=sd(w)) %>%
  mutate(w = paste(format(round(Meanw, 2), nsmall=2), 
                   "(", format(round(SDw, 2), nsmall=2), ")", sep="")) %>%
  select(-Meanw, -SDw) %>%
  pivot_wider(values_from = w, names_from = model)
Data_Acc <- Data %>% group_by(participant, SAT, condition) %>%
  summarise(Acc = mean(correct), .groups="drop") %>%
  reframe(summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("SAT", "condition"))) %>%
  rename(sewithin = se) 
Accuracies <- Data_Acc %>% 
  select(SAT, condition, Acc, sewithin) %>%
  # mutate(Acc=paste(format(round(Acc, 2), nsmall=2),"(", format(round(sewithin, 2), nsmall=2), ")", sep="")) %>%
  mutate(Acc=format(round(Acc, 2), nsmall=2)) %>%
  select(-sewithin) %>%
  mutate(condition=as.character(condition)) %>%
  pivot_wider(id_cols=c("SAT"), values_from = "Acc", names_from = condition)
Accuracies <- Accuracies %>% cbind(fitted_w)
dir.create("figures", showWarnings = FALSE)
write.table(Accuracies, file="figures/data_accuracies_SAT.txt",
            sep = ",", quote = FALSE, row.names = F)

## Suppl Fig 5: Plot of Fitted Accuracy                    ----
Data_Acc <- Data %>% group_by(participant, SAT, condition) %>%
  summarise(Acc = mean(correct), .groups="drop") %>%
  reframe(summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("SAT", "condition"))) %>%
  rename(sewithin = se) 

Preds_Acc <- Preds_RatingDist_corr_cond %>% 
  group_by(model, SAT, condition) %>%
  summarise(Acc = sum(p*correct)/sum(p), .groups="drop") 
p_Acc <- ggplot(Data_Acc,
                aes(x=condition, y=Acc)) +
  geom_line(data=Preds_Acc, aes(x = condition, linetype="Predicted", group=model), linewidth=1)+
  #geom_line(linetype="dashed", alpha=0.5,position = pd)+
  geom_errorbar(aes(ymin=Acc-sewithin, ymax=Acc+sewithin), colour="black", width=.2, linewidth=0.4) +
  geom_point(fill="white", aes(shape="Observed"))+
  facet_grid(rows=vars(model), 
             cols=vars(SAT), scales = "free_x", drop=TRUE)+ #, dir="v"
  ylab("Mean accuracy")+
  scale_x_discrete(
    name="Longer line length [mm]")+
  scale_linetype_manual(name="", values=1) +
  scale_shape_manual(values=c(21),name = "")  +
  scale_y_continuous(breaks=c(0.6, 0.7, 0.8, 0.9, 1.0), name="Mean Accuracy")+
  theme_bw() +
  guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        #axis.text.x = element_text(angle=90,hjust = 0.5, vjust = 0.5),
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
ggsave("../../../Supplement/figures/modelAccuracy1_SAT.eps",
       width = 17.62, height=17, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/modelAccuracy1_SAT.tiff",
       width = 17.62, height=9/0.75, units="cm",dpi=600)

## Suppl Table 2: Descriptive parameter fits dynaViTE      ----
meanFits <- fits %>% 
  select(-c("participant", "k", "N", "AICc","BIC", "AIC", "negLogLik","fixed")) %>%
  group_by(model) %>% 
  summarise(across(.cols = everything(),
                   ~paste(round(mean(.x[abs(.x)<1e+24], na.rm=TRUE),2), 
                          " (", round(sd(.x[abs(.x) <1e+24], na.rm = TRUE), 2), ")", sep=""))) %>%
  select(c("model", paste0("v", 1:6), "sv", "a1","a2", "sz", paste0("theta", 1:5),
           "t0", "tau0", "w","sigvis", "svis", "lambda"))
names(meanFits)[2:7] <- paste0("$\\nu_", 1:6, "$")
names(meanFits)[8] <- paste0("$s\\nu$")
names(meanFits)[c(9, 10)] <- paste0("$a_", 1:2, "$")
names(meanFits)[11] <- paste0("$s_z$")
names(meanFits)[12:16] <- paste0("$\\theta_{", 1:5, "}$")
names(meanFits)[17] <- c("$t_0$")
names(meanFits)[18] <- paste0("$\\tau_0$")
names(meanFits)[20:21] <- c("$\\sigma_V$", "$s_V$")
names(meanFits)[22] <- paste0("$\\lambda$")
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

# # Use this to produce a long table with only fits for dynaViTE
# meanFits_long <- meanFits %>% 
#   # Report only dynaViTE model  
#   subset(model=="dynaViTE") %>% 
#   select(-model) %>%
#   pivot_longer(cols = 1:(ncol(meanFits)-1), names_to = "Parameter") 
# meanFits_long
# 
# # Use flextable for word output
# if (!"flextable" %in% installed.packages()[,"Package"]) {
#   install.packages("flextable")
# } 
# library(flextable)
# Fits_wordtable <- flextable(meanFits_long) %>%
#   align(align="left", j=1) %>%
#   align(align="center", j=2:4)
# save_as_docx("Parameter table"=Fits_wordtable,
#              path="../Draft/tablepars.docx")
# print(Fits_wordtable, preview="docx")


# Use xtable for latex output
if (!"xtable" %in% installed.packages()[,"Package"]) {
  install.packages("xtable")
}
library(xtable)
# Use this to produce a long table with all model fits
meanFits_long <- meanFits %>% 
  # Report only dynaViTE model  
  pivot_longer(cols = 2:(ncol(meanFits)), names_to = "Parameter") %>%
  pivot_wider(id_cols = "Parameter", names_from = c("model"), values_from = value) 
meanFits_long <- meanFits_long[,c(1, 5, 4, 3, 2)]
meanFits_long
Fits_textable <- xtable(meanFits_long, align = c("l", rep("c", ncol(meanFits_long))),
                        caption="\\raggedright Mean and standard deviation of parameter fits for experiment 4 ($t_0$ and $\\tau_0$ measured in seconds)")
addtorow <- list()
addtorow$pos <- list(-1, c(0:(nrow(meanFits_long)-1)))
addtorow$command <- c('\\toprule  & \\multicolumn{4}{c}{Model}\\\\ \\cmidrule{2-5} ','[1mm]')

print(Fits_textable, type="latex", 
      sanitize.text.function=function(x){x}, 
      include.rownames=FALSE, 
      add.to.row=addtorow, 
      hline.after = c(0, nrow(Fits_textable)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")
print(Fits_textable, type="latex", 
      file="../../../Supplement/figures/tablepars2_SAT.tex",
      sanitize.text.function=function(x){x}, 
      include.rownames=FALSE, 
      add.to.row=addtorow, 
      hline.after = c(0, nrow(Fits_textable)), booktabs = TRUE,
      caption.placement="top",
      table.placement="hp")

## Suppl Fig 9: Fitted response distributions            ----
model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")
p_ratingdist_SAT <- ggplot(data=Data_RatingDist_corr_cond, aes(x=rating, y=p))+
  geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
  geom_point(data=Preds_RatingDist_corr_cond,
             aes(shape=model, fill=model), 
             position=position_dodge(1))+
  facet_nested(condition~SAT+correct,
               labeller = labeller(correct=c("0"="Wrong", "1"="Correct"),
                                   condition=function(x) paste(x, "mm")))+
  scale_x_discrete(name="Confidence", breaks=1:6, 
                   labels = seq(50, 100, by=10))+
  scale_y_continuous(name = "Probability", breaks = c(0, 0.25, 0.5, 0.75))+
  scale_shape_manual(name = "Model prediction",values=c(21, 22,23,24)) +
  scale_fill_manual(name = "Model prediction",values=model_colors) +
  theme_bw() +
  theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        #axis.text.x = element_text(angle=90),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        legend.key = element_blank(),
        legend.position = "bottom", legend.text = element_text(size=12, family = "Times"),
        strip.text.y = element_text(angle = 0, size=9, family = "Times"),
        panel.spacing=unit(0, "lines"))+
  coord_cartesian(xlim = c(0.5,6.5),ylim=c(0, 0.82),  clip = 'off')
p_ratingdist_SAT

p_ratingdist_SAT+
  theme(plot.margin = margin(0, 0, 0, 0, "cm"),
        legend.box = "vertical", legend.spacing.y = unit(-.2, "cm"))
ggsave(paste0("figures/RespDist_SAT.tiff"),
       width=17.3, height=18, units="cm",dpi=1200)
ggsave(paste0("../../../Supplement/figures/RespDist_SAT.eps"),
       width=17.3, height=18, units="cm",dpi=1200, device = cairo_ps)

# The end                                                  ####