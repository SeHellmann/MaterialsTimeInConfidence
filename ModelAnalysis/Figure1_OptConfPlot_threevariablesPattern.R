## Visualization of the posterior probability of being correct
# as a function of decision time, decision evidence, and visibility evidence
# for different distribution for the levels of stimulus discriminability:
# Figure 1: Discrete uniformly distribution (most experimental manipulations)
# Suppl. Fig. 1: Continuously uniform distribution
# Suppl. Fig. 2: Folded normal distribution

# This code produces Figure 1 and
# Supplementary Figures 1 and 2 in the paper. 

# Preamble and imports                                     ----
rm(list = ls())
# use RStudio to find script file path
script_path <- dirname(rstudioapi::getSourceEditorContext()$path)
setwd(script_path)
print("Working directory set to:")
print(getwd())

library(tidyverse)
library(ggplot2)
library(ggpubr)

# Load fitted parameters                                   ----
load("../EmpiricalStudy/modelfits_threeexps.RData")
paramDf <- subset(fits, participant==215&model=="dynaViTE") %>%
  select(t0:lambda, v4, v5)
paramDf <- round(paramDf, 1)

text_size <- 17
color_legend_name = "Probability correct"
contour_breaks = 2:9/10
dir.create("figures",showWarnings = FALSE)

#########################################################################
#######   Figure 1: Discrete uniform distribution           #############
#########################################################################   
Ds <- as.numeric(paramDf[,paste("v", 1:5, sep="")])
odds_fct <- function(RT, X, V, Ds, a, z, sv, tau, svis, sigvis, ...) {
      SigV <- svis^2+(RT+tau)*sigvis^2
      Sc <- sv^2*(RT+tau)+1
      X <- -X
      TSSV <- (RT+tau)*(Sc+SigV)
      exp1 <- exp(-2*V*X/TSSV)
      v1<- rep(0, length(RT))
      v2 <- rep(0, length(RT))
      for (i in 1:length(Ds)) {
        v1 <- v1 + exp(-(TSSV*Ds[i]-V*Sc + X*SigV)/(2*Sc*SigV*TSSV))
        v2 <- v2 + exp(-(TSSV*Ds[i]-V*Sc - X*SigV)/(2*Sc*SigV*TSSV))
      }
      exp1*v1/v2
    }

input <- list(Decmin=-2, Decmax=6, Decsteps=30, Decfixed=3,
              Vismin=0, Vismax=6, Vissteps=30, Visfixed=1,
              RTmin=0.2, RTmax=3.5, RTsteps=30, RTfixed=2)
input <- c(input, paramDf)

predsRT <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                       vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                       rt = input$RTfixed)
predsDec <- expand.grid(dec=input$Decfixed, 
                        seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps))
predsVis <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=input$Visfixed,
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps))
predsRT$odds_correct <- odds_fct(predsRT$rt, predsRT$dec, predsRT$vis, 
                                 min = input$minD, max=input$maxD,
                                 Ds = Ds,
                                 a=input$a, sv=input$sv, 
                                 sigvis=input$sigvis, svis=input$svis, 
                                 tau=input$tau)
predsDec$odds_correct <- odds_fct(predsDec$rt, predsDec$dec, predsDec$vis, 
                                  min = input$minD, max=input$maxD,
                                  Ds = Ds,
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsVis$odds_correct <- odds_fct(predsVis$rt, predsVis$dec, predsVis$vis, 
                                  min = input$minD, max=input$maxD,
                                  Ds = Ds,
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsVis$fill <- predsVis$odds_correct/(1+predsVis$odds_correct)
predsDec$fill <- predsDec$odds_correct/(1+predsDec$odds_correct)
predsRT$fill <- predsRT$odds_correct/(1+predsRT$odds_correct)

minp <- min(predsRT$fill, predsDec$fill, predsVis$fill)
gRT <- ggplot(data = predsRT, aes(x = dec, y = vis, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final decision state", y = "Final visibility state",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.box.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
gVis <- ggplot(data = predsVis, aes(x = dec, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final decision state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  theme(plot.margin = margin(0.1, .2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.box.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
gDec <- ggplot(data = predsDec, aes(x = vis, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final visibility state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  # switch off minor gridlines
        panel.grid.major = element_blank(),
        strip.text = element_text(size=9),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.box.margin = margin(0.1, 0, 0, 0, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
ggarrange(gRT, gVis, gDec, nrow=1, common.legend = TRUE, legend="bottom")
ggsave("figures/posterior_prob_corrects.tiff",
       width = 13.62, height=6, units="cm",dpi=900)


#########################################################################
#######  Suppl. Fig. 2_ Coninuous uniform distribution          #########
#########################################################################        
Ds <- matrix(c(1.0, 1.2, 
              0.3, 2.4, 
              0.1, 3.7), nrow=3, byrow=TRUE)
odds_fct <- function(RT, X, V, min, max, a, z, sv, tau, svis, sigvis, ...) {
  SigV <- svis^2+(RT+tau)*sigvis^2
  Sc <- sv^2*(RT+tau)+1
  X <- -X
  denom <- sqrt(Sc*SigV*(RT+tau)*(Sc+SigV))
  exp(-2*V*X/((RT+tau)*(Sc+SigV))) * 
    (pnorm(((RT+tau)*(Sc+SigV)*max - V*Sc + X*SigV)/denom)-pnorm(((RT+tau)*(Sc+SigV)*min - V*Sc + X*SigV)/denom))/
    (pnorm(((RT+tau)*(Sc+SigV)*max - V*Sc - X*SigV)/denom)-pnorm(((RT+tau)*(Sc+SigV)*min - V*Sc - X*SigV)/denom))
}

input <- list(Decmin=-2, Decmax=4, Decsteps=30, Decfixed=1.5,
              Vismin=0, Vismax=8, Vissteps=30, Visfixed=1,
              RTmin=0.2, RTmax=3.5, RTsteps=30, RTfixed=2,
              plot_val="prob")
input <- c(input, paramDf)

predsRT <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                       vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                       rt = input$RTfixed,
                       range=1:3)
predsDec <- expand.grid(dec=input$Decfixed, 
                        seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps),
                        range=1:3)
predsVis <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=input$Visfixed,
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps),
                        range=1:3)
predsRT$odds_correct <- odds_fct(predsRT$rt, predsRT$dec, predsRT$vis, 
                                 min = Ds[predsRT$range, 1], max=Ds[predsRT$range, 2],
                                 a=input$a, sv=input$sv, 
                                 sigvis=input$sigvis, svis=input$svis, 
                                 tau=input$tau)
predsDec$odds_correct <- odds_fct(predsDec$rt, predsDec$dec, predsDec$vis, 
                                  min = Ds[predsDec$range, 1], max=Ds[predsDec$range, 2],
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsVis$odds_correct <- odds_fct(predsVis$rt, predsVis$dec, predsVis$vis, 
                                  min = Ds[predsVis$range, 1], max=Ds[predsVis$range, 2],
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsDec$Range <- paste0("[", Ds[predsDec$range,1], ", ", Ds[predsDec$range,2], "]")
predsVis$Range <- paste0("[", Ds[predsVis$range,1], ", ", Ds[predsVis$range,2], "]")
predsRT$Range <- paste0("[", Ds[predsRT$range,1], ", ", Ds[predsRT$range,2], "]")
predsVis$fill <- predsVis$odds_correct/(1+predsVis$odds_correct)
predsDec$fill <- predsDec$odds_correct/(1+predsDec$odds_correct)
predsRT$fill <- predsRT$odds_correct/(1+predsRT$odds_correct)

minp <- min(predsRT$fill, predsDec$fill, predsVis$fill)
gRT <- ggplot(data = predsRT, aes(x = dec, y = vis, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 2, 4, 6)) +
  labs(x = "Final decision state", y = "Final visibility state",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(Range), labeller=label_both)+
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        strip.text.y = element_blank() , 
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        legend.spacing.x = unit(1,"line"),
        panel.spacing=unit(0, "lines"))
gVis <- ggplot(data = predsVis, aes(x = dec, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final decision state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(Range), labeller=label_both)+
  theme(plot.margin = margin(0.1, .2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text.y = element_blank() , 
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        legend.spacing.x = unit(1,"line"),
        panel.spacing=unit(0, "lines"))
gDec <- ggplot(data = predsDec, aes(x = vis, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final visibility state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(Range), labeller=label_both)+
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text = element_text(size=12),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.key = element_blank(),
        legend.key.width=unit(2,"line"),
        legend.spacing.x = unit(1,"line"),
        panel.spacing=unit(0, "lines"))
ggarrange(gRT, gVis, gDec, nrow=1, common.legend = TRUE, legend="bottom",
                  widths=c(0.32, 0.32, 0.36))
ggsave("figures/posterior_prob_corrects_unifdist.tiff",
      width = 17.62, height=16, units="cm",dpi=900, bg="white")
# ## Only for manuscript generation
# ggsave("C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/posterior_prob_corrects_unifdist.eps",
#        width = 17.62, height=15, units="cm",dpi=1200, device = cairo_ps)



#########################################################################
#######  Suppl. Fig. 3: Absolute normal distribution         ############
#########################################################################        
SDs <- c(1,2,3)
odds_fct <- function(RT, X, V,sd, a, z, sv, tau, svis, sigvis, ...) {
  SigV <- svis^2+(RT+tau)*sigvis^2
  Sc <- sv^2*(RT+tau)+1
  X <- -X
  TSSV <- (RT+tau)*(Sc+SigV)
  denom <- sqrt(sd^2*Sc*SigV*(Sc*SigV+sd^2*TSSV))
  exp(-2*V*X*(sd^2/(sd^2*TSSV+Sc*SigV)))*
    pnorm(sd^2*(V*Sc-X*SigV)/denom) / pnorm(sd^2*(V*Sc+X*SigV)/denom)
}

input <- list(Decmin=-2, Decmax=4, Decsteps=30, Decfixed=1.5,
              Vismin=0, Vismax=8, Vissteps=30, Visfixed=1,
              RTmin=0.2, RTmax=3.5, RTsteps=30, RTfixed=2,
              plot_val="prob")
input <- c(input, paramDf)

predsRT <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                       vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                       rt = input$RTfixed,
                       SD=c(0.25, 1, 1.75))
predsDec <- expand.grid(dec=input$Decfixed, 
                        seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=seq(input$Vismin, input$Vismax, length.out=input$Vissteps),
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps),
                        SD=c(0.25, 1, 1.75))
predsVis <- expand.grid(dec=seq(input$Decmin, input$Decmax, length.out=input$Decsteps),
                        vis=input$Visfixed,
                        rt = seq(input$RTmin, input$RTmax, length.out=input$RTsteps),
                        SD=c(0.25, 1, 1.75))
predsRT$odds_correct <- odds_fct(predsRT$rt, predsRT$dec, predsRT$vis, 
                                 sd=predsRT$SD,
                                 a=input$a, sv=input$sv, 
                                 sigvis=input$sigvis, svis=input$svis, 
                                 tau=input$tau)
predsDec$odds_correct <- odds_fct(predsDec$rt, predsDec$dec, predsDec$vis, 
                                  sd=predsDec$SD,
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsVis$odds_correct <- odds_fct(predsVis$rt, predsVis$dec, predsVis$vis, 
                                  sd=predsVis$SD,
                                  a=input$a, sv=input$sv, 
                                  sigvis=input$sigvis, svis=input$svis, 
                                  tau=input$tau)
predsDec$SD <- paste0("SD=", predsDec$SD)
predsVis$SD <- paste0("SD=", predsVis$SD)
predsRT$SD <- paste0("SD=", predsRT$SD)
predsVis$fill <- predsVis$odds_correct/(1+predsVis$odds_correct)
predsDec$fill <- predsDec$odds_correct/(1+predsDec$odds_correct)
predsRT$fill <- predsRT$odds_correct/(1+predsRT$odds_correct)

minp <- min(predsRT$fill, predsDec$fill, predsVis$fill)
gRT <- ggplot(data = predsRT, aes(x = dec, y = vis, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), breaks=c(0, 2, 4, 6)) +
  labs(x = "Final decision state", y = "Final visibility state",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(SD))+
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(), 
        panel.grid.major = element_blank(),
        strip.text.y = element_blank() , 
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.spacing.x = unit(1,"line"),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
gVis <- ggplot(data = predsVis, aes(x = dec, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final decision state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  stat_contour(color = "red", breaks = 0.5,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(SD))+
  theme(plot.margin = margin(0.1, .2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        strip.text.y = element_blank() , 
        strip.background = element_blank(),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.spacing.x = unit(1,"line"),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
gDec <- ggplot(data = predsDec, aes(x = vis, y = rt, z = fill))  +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(x = "Final visibility state", y = "Decision time",
       fill = color_legend_name)+
  stat_contour()+
  geom_tile(aes(fill = fill)) +
  stat_contour(color = "white", breaks = contour_breaks,linewidth = .7) +#, bins = 15 
  scale_fill_viridis_c(breaks=contour_breaks, limits=c(minp, 1),
                       guide = guide_colourbar(ticks.linewidth = .7))+
  theme_bw() +
  facet_grid(rows=vars(SD))+
  theme(plot.margin = margin(0.1, 0.2, 0, 0, "cm"),
        text = element_text(size=9, family="Times"),
        axis.text = element_text(size=9, family="Times", color="black"),
        axis.title.y = element_text( angle=90),
        panel.grid.minor = element_blank(),  
        panel.grid.major = element_blank(),
        strip.text = element_text(size=12),
        legend.position = "bottom",
        legend.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.box.margin = margin(0.1, 0.1, 0, 0.1, "cm"),
        legend.spacing.x = unit(1,"line"),
        legend.key.width=unit(2,"line"),
        panel.spacing=unit(0, "lines"))
ggarrange(gRT, gVis, gDec, nrow=1, common.legend = TRUE, legend="bottom",
                  widths=c(0.32, 0.32, 0.36))
ggsave("figures/posterior_prob_corrects_absnorm.tiff",
       width = 17.62, height=16, units="cm",dpi=900, bg="white")
# ## Only for manuscript generation
# ggsave("C:/Users/PPA859/Documents/Manuskripte/TimeInDynWEV/Supplement/figures/posterior_prob_corrects_absnorm.eps",
#        width = 17.62, height=15, units="cm",dpi=1200, device = cairo_ps)

