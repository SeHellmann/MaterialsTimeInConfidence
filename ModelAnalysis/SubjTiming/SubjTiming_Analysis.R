## Simulation of a dynaViTE model with subjective timing measure
#
# This simulation should demonstrate that the OU-model itself is not able to
# produce a double-increase pattern

# This code produces Supplementary Figures 10.

# Preamble and imports    
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())
Rcpp::sourceCpp("SimulatedynaViTESubjectiveTiming.cpp")
source("simulatedynaViTEsubjTming_part.R")
library(Rmisc)
library(tidyverse)
library(ggpubr)
library(ggh4x)
conflicted::conflict_prefer_all("dplyr", "plyr")
conflicted::conflicts_prefer(dplyr::mutate)
conflicted::conflicts_prefer(dplyr::filter)

windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
dir.create("figures",showWarnings = FALSE)

load("modelfits_threeexps.RData")
load("collected_Data_Fits_Predicts.RData")
Preds_MRating_corr_cond <- subset(Preds_MRating_corr_cond, model=="dynaViTE")
Preds_RatingDist_corr_cond <- Preds_RatingDist_corr_cond %>% filter(model=="dynaViTE")
Preds_RTdens_corr_cond_rating <- Preds_RTdens_corr_cond_rating %>% filter(model=="dynaViTE")
Preds_RTQuants_corr_rating <- Preds_RTQuants_corr_rating %>% filter(model=="dynaViTE")



dynaViTE_fits <- subset(fits, model=="dynaViTE")
# Use fitted parameters from Hawkins & Heathcote (2021, PsychRev)
# from the perceptual task and the neutral condition
dynaViTE_fits$Timerate <- 1.12
dynaViTE_fits$Timediff <- 0.51
if (!file.exists("simulation_subjtiming.RData")) {
  #nrow_parts <- table(Data$participant)
  set.seed(2904)
  simulations_subjtiming <- dynaViTE_fits %>% group_by(participant, experiment) %>%
    do(simulatedynaViTE_subjtimig(select(.,where(~!all(is.na(.x)))), n=1e+5, simult_conf=TRUE))
  
  recompute_rating <- function(conf, cur_participant) {
    cur_participant <- cur_participant$participant[1]
    print(cur_participant)
    cur_participant_Data <- Data %>% filter(participant==cur_participant) 
    proportions <- table(cur_participant_Data$rating)/nrow(cur_participant_Data)
    thresholds <- c(quantile(conf, probs = cumsum(proportions))[1:4], Inf)
    new_rating <- as.numeric(as.factor((cut(conf, c(-Inf, thresholds), include.lowest = TRUE))))
    return(new_rating)
  }
  
  simulations_subjtiming <- simulations_subjtiming %>% group_by(participant) %>%
    mutate(rating_old = rating,
           rating = recompute_rating(conf, cur_group()))
  
  #### Aggregate simulations     ####
  levels(Preds_MRating_corr_cond$condition)
  simulations_subjtiming <- simulations_subjtiming %>% 
    mutate(condition=ifelse(grepl("Experiment 1", experiment), condition, 
                            ifelse(grepl("Experiment 2", experiment), condition + 5, condition+10))) %>%
    mutate(condition = factor(condition, levels=1:13, labels=levels(Preds_MRating_corr_cond$condition))) %>%
    group_by(experiment, participant, condition) %>%
    mutate(nrows=n())
  Sim_RatingDist_part <- simulations_subjtiming %>% 
    group_by(experiment, participant, condition, rating, correct) %>% 
    summarise(p = n()/(mean(nrows)),.groups = "drop") %>%
    full_join(y = merge(distinct(simulations_subjtiming[,c("correct","rating", "experiment")]),
                        distinct(simulations_subjtiming[,c("condition", "participant", "experiment")]), by="experiment"),
              by = c("experiment", "participant", "condition", "rating", "correct")) %>%
    mutate(p = ifelse(is.na(p), 0, p)) 
  
  ### Sanity Checks:
  # sum(Sim_RatingDist_part$p)
  # table((Sim_RatingDist_part %>% group_by(condition, participant) %>% summarise(p = sum(p)))$p)
  
  # For the plots we won't differentiate between 
  # stimulus directions and participants
  Sim_RatingDist_corr_cond <- Sim_RatingDist_part %>% 
    group_by(experiment, correct, condition, rating) %>% 
    summarise(p = mean(p), .groups = "drop")
  #sum(Sim_RatingDist_corr_cond$p)
  
  ### b) simulations_subjtiming mean confidence rating                       ----
  Sim_MRating_corr_cond <- simulations_subjtiming %>% 
    group_by(participant, experiment, condition, correct) %>%
    summarise(MRating = mean(rating), .groups = "drop") 
  Sim_MRating_corr_cond_SE <- Sim_MRating_corr_cond %>%
    group_by(experiment) %>%
    reframe(summarySEwithin(pick(everything()),measurevar = "MRating", 
                            idvar = "participant", withinvars = c("condition", "correct"))) %>%
    select(experiment, condition, correct, sewithin = se) 
  Sim_MRating_corr_cond <- Sim_MRating_corr_cond %>%
    group_by(experiment, condition, correct) %>%
    summarise(MRating=mean(MRating, na.rm=TRUE), .groups="drop") %>%
    mutate(correct=as.factor(correct), condition=as.factor(condition)) %>%
    left_join(Sim_MRating_corr_cond_SE, by=c("condition", "correct", "experiment"))
  ### c) simulations_subjtiming response time quantiles                      ----  
  # Reaction Time Quantiles of the simulations_subjtiming grouped by rating and accuracy  
  Sim_RTQuants_corr_rating <- simulations_subjtiming %>%
    group_by(experiment, rating, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  # Reaction Time Quantiles of the simulations_subjtiming grouped by condition and accuracy  
  Sim_RTQuants_corr_cond <- simulations_subjtiming %>%
    group_by(experiment, condition, correct) %>%
    reframe(p=c(.1,.3,.5,.7,.9), q = quantile(rt, probs = c(.1,.3,.5,.7,.9))) 
  
  Sim_Acc <- simulations_subjtiming %>% group_by(participant, experiment, condition) %>%
    summarise(Acc = mean(correct), .groups="drop") %>%
    group_by(experiment) %>%
    reframe(summarySEwithin(pick(everything()),measurevar = "Acc", idvar = "participant", withinvars = c("condition"))) %>%
    rename(sewithin = se) 
  save(file="raw_simulation_subjtiming.RData", simulations_subjtiming)
  save(file="simulation_subjtiming.RData", 
       Sim_MRating_corr_cond, Sim_MRating_corr_cond_SE,
       Sim_RatingDist_corr_cond, Sim_RatingDist_part,Sim_Acc,
       Sim_RTQuants_corr_cond, Sim_RTQuants_corr_rating, dynaViTE_fits)
  rm(simulations_subjtiming)
} else {
  load("simulation_subjtiming.RData")
}






# B  Generate visualizations                               ----
## Include font style and define colors for plots          ----
# Use: extrafont::loadfonts() to permanently add Times 
#windowsFonts(Times=windowsFont("Times New Roman"))
two_colors_correct <- c("#1b9e77", "#fc8d62")
model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")

## Figure 4: Plot of Mean Ratings                          ----
pd <- position_dodge(0.02)
Preds_MRating_corr_cond_plot <- mutate(Preds_MRating_corr_cond,
                                       confidence ="Confidence measure") 
# Sim_MRating_corr_cond_plot <- merge(Sim_MRating_corr_cond, 
#                                      expand.grid(model=c("dynaViTE", "dynWEV", "2DSD+", "2DSD"), 
#                                                  confidence ="Confidence measure"))


p_MRating <- ggplot(Sim_MRating_corr_cond,
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
# # Only for manuscript generation
ggsave("../../../Supplement/figures/simul_subjtiming_meanconf.eps",
      width = 17.62, height=9/0.7, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/simul_subjtiming_meanconf.tiff",
       width = 18, height=13, units="cm",dpi=600)



## Figure 6: RTQuantiles accross correct X rating          ----
Sim_RTQuants_corr_rating_plot <- Sim_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))))%>% 
  filter(p %in% c(0.1, 0.5, 0.9))
Preds_RTQuants_corr_rating_plot <- Preds_RTQuants_corr_rating %>%
  mutate(X = factor(paste(rating, correct, sep="."),
                    levels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep=".")),
                    labels = c(paste(c(5:1, 1:5), rep(0:1, each=5), sep="."))),
         model = factor(model, levels=c("dynaViTE", "dynWEV", "2DSD+", "2DSD"),
                        labels=c("dynaViTE", "dynWEV", "2DSD+", "2DSD")))%>% 
  filter(p %in% c(0.1, 0.5, 0.9))
ggplot()+
  geom_line(data=mutate(Preds_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
                        rating = as.factor(rating)),
            aes(x=rating, y=log(q), group=as.factor(p),color=correct), linewidth=0.7)+
  geom_point(data=mutate(Sim_RTQuants_corr_rating_plot, correct=factor(correct, labels=c("Wrong", "Correct")),
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
#  # # Only for manuscript generation
# ggsave("../../../Supplement/figures/simul_subjtiming_RTQuantsConf.eps",
#        width = 15, height=13, units="cm",dpi=600, device = cairo_ps)
ggsave("figures/simul_subjtiming_RTQuantsConf.tiff",
       width = 15, height=13, units="cm",dpi=600)
# 
# ## Suppl Fig 6-8: Fitted response distributions            ----
# model_colors = c("#D81B60", "#1E88E5", "#FFC107","#004D40")
# for (i in 1:3) {
#   cur_experiment <- sort(unique(fits$experiment))[i]
#   p_ratingdist_KonfMask <- ggplot(data=filter(Sim_RatingDist_corr_cond, experiment==cur_experiment), aes(x=rating, y=p))+
#     geom_bar(stat = "identity", show.legend = FALSE, fill="white", col="black")+
#     geom_point(data=filter(Preds_RatingDist_corr_cond, experiment==cur_experiment),
#                aes(shape=model, fill=model), 
#                position=position_dodge(1))+
#     facet_nested(condition~experiment+correct,
#                  labeller = labeller(correct=c("0"="Wrong", "1"="Correct"),
#                                      condition=function(x) paste(x, "ms")))+
#     scale_x_continuous(name = "Confidence [% scale width]", breaks = 1:5,
#                        labels = c("0-20","20-40", "40-60", "60-80", "80-100")) +
#     scale_y_continuous(name = "Probability", breaks = c(0, 0.4, 0.8))+
#     scale_shape_manual(name = "Model prediction",values=c(21, 22,23,24)) +
#     scale_fill_manual(name = "Model prediction",values=model_colors) +
#     theme_bw() +
#     theme(plot.margin = margin(0, 0.2, 0, 0, "cm"),
#           text = element_text(size=9, family="Times"),
#           axis.text = element_text(size=9, family="Times", color="black"),
#           axis.text.x = element_text(angle=90),
#           axis.title.y = element_text( angle=90),
#           panel.grid.minor = element_blank(),  # switch off minor gridlines
#           panel.grid.major = element_blank(),
#           legend.key = element_blank(),
#           legend.position = "bottom", legend.text = element_text(size=12, family = "Times"),
#           strip.text.y = element_text(angle = 0, size=9, family = "Times"),
#           panel.spacing=unit(0, "lines"))+
#     coord_cartesian(xlim = c(0.5,5.5),ylim=c(0, 0.9),  clip = 'off')
#   show(p_ratingdist_KonfMask+
#     theme(plot.margin = margin(0, 0, 0, 0, "cm"),
#           legend.box = "vertical", legend.spacing.y = unit(-.2, "cm")))
#   
#   # ggsave(paste0("figures/RespDist", 
#   #               str_replace_all(cur_experiment, c("\n"="", " "="", "\\("="", "\\)"="", "\\."="")),
#   #               ".tiff"),
#   #        width=17.3, height=18, units="cm",dpi=1200)
#   # ggsave(paste0("C:/Users/PPA859/OneDrive - ku.de/Forschung/Manuskripte/TimeInDynWEV/Supplement/figures/RespDist", 
#   #               str_replace_all(cur_experiment, c("\n"="", " "="", "\\("="", "\\)"="", "\\."="")),
#   #               ".eps"),
#   #        width=17.3, height=18, units="cm",dpi=1200, device = cairo_ps)
#   
# }






# 
# Preds_Acc <- Preds_RatingDist_corr_cond %>% 
#   group_by(model, experiment, condition) %>%
#   summarise(Acc = sum(p*correct)/sum(p), .groups="drop") 
# p_Acc <- ggplot(Sim_Acc,
#                 aes(x=condition, y=Acc)) +
#   geom_line(data=Preds_Acc, aes(x = condition, linetype="Predicted", group=model), linewidth=1)+
#   #geom_line(linetype="dashed", alpha=0.5,position = pd)+
#   geom_errorbar(aes(ymin=Acc-sewithin, ymax=Acc+sewithin), colour="black", width=.2, linewidth=0.4) +
#   geom_point(fill="white", aes(shape="Observed"))+
#   facet_grid(rows=vars(model), 
#              cols=vars(experiment), scales = "free_x", drop=TRUE)+ #, dir="v"
#   ylab("Mean accuracy")+
#   scale_x_discrete(
#     name="Stimulus-onset-asynchrony [ms]                  Motion coherence [%]                                     Contrast [%]            ")+
#   scale_linetype_manual(name="", values=1) +
#   scale_shape_manual(values=c(21),name = "")  +
#   scale_y_continuous(breaks=c(0.5, 0.7, 0.9), name="Mean Accuracy")+
#   theme_bw() +
#   guides(shape=guide_legend(order=3), color=guide_legend(order=3))+
#   theme(plot.margin = margin(0, 0, 0, 0, "cm"),
#         axis.text.x = element_text(angle=90,hjust = 0.5, vjust = 0.5),
#         text = element_text(size=9, family="Times"),
#         axis.text = element_text(size=9, family="Times", color="black"),
#         axis.title.y = element_text( angle=90),
#         panel.grid.minor = element_blank(),  # switch off minor gridlines
#         panel.grid.major = element_blank(),
#         strip.text = element_text(size=9),
#         legend.position = "bottom",
#         # legend.position=c(.80, .005), legend.justification = c(0,0),
#         legend.margin = margin(-0.1, 0, 0, 0, "cm"),#legend.direction = "horizontal",
#         legend.box.margin = margin(-0.1, 0, 0, 0, "cm"),
#         legend.key = element_blank(),
#         legend.key.width=unit(2,"line"),
#         panel.spacing=unit(0, "lines"))
# p_Acc
# 







# 
# 
# load("raw_predictions.RData")
# rm("preds_RT")
# Sim_p_ratings <- simulations_subjtiming %>% 
#   group_by(participant, experiment) %>% 
#   mutate(nrow=n()) %>%
#   group_by(participant, experiment, rating)%>%
#   summarise(prop_rating = n()/nrow[1], .groups="drop") 
# sum(Sim_p_ratings$prop_rating)
# length(unique(preds_confidence$participant))
# for (i in unique(preds_confidence$participant)) {
#   temp <- subset(preds_confidence, model=="dynaViTE" & participant == i) 
#   temp <- temp %>% group_by(participant, experiment, rating) %>%
#     summarise(p = sum(p)/length(unique(condition)/2)) 
#   sum(temp$p)
# }
# Preds_p_ratings <- preds_confidence %>%
#   filter(model=="dynaViTE") %>%
#   group_by(participant, experiment, rating) %>%
#   summarise(prop_rating_pred = sum(p)/(2*length(unique(condition))), .groups="drop") 
# Preds_p_ratings <- left_join(Preds_p_ratings, Sim_p_ratings)
# 
# ggplot(Preds_p_ratings, aes(x=prop_rating, y=prop_rating_pred))+
#   geom_point()
