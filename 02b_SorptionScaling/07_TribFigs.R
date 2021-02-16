library(mgcv)
library(tidyverse)
library(lubridate)
library(dataRetrieval)
library(readr)
library(readxl)
library(zoo)
library(openxlsx)
library(discharge)
library(ggpubr)
library(gridExtra)
library(egg)
library(grid)


#load functions
source("02b_SorptionScaling/zz_TravelTimeEstFuns.R")

# get P sorbed data
# get P sorbed/gDM above flow target
StartStumpf <- as.POSIXct("2019-03-01", format = "%Y-%m-%d")
EndStumpf <- as.POSIXct("2019-07-01", format = "%Y-%m-%d")

WC <- read.csv("04_generatedData/05d_WCboot.csv") 
UTLC <- read.csv("04_generatedData/05d_UTLCboot.csv")
SFT <- read.csv("04_generatedData/05d_SFTboot.csv") 


  

WC2 <- WC %>%   
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
  filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
  mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
         Svol_gPwindow.hf = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0),
         DRPloadDS.gPwin.wS.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadDS.gPwin.wS, 0)) %>% 
  select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, DRPloadDS.gPwin.wS, DRPloadDS.gPwin.wS.hf,
         Svol_gPwindow.hf, 
         SorbDistM, gTarPerFlow, "Qm3m" = "Qm3m.l1", bootN) %>% 
  mutate(stream = "WC") %>% 
  filter(!is.na(gTarPerFlow))



UTLC2 <- UTLC %>%  
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
    filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
    mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
           Svol_gPwindow.hf = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0),
           DRPloadDS.gPwin.wS.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadDS.gPwin.wS, 0)) %>% 
    select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, DRPloadDS.gPwin.wS, DRPloadDS.gPwin.wS.hf,
           Svol_gPwindow.hf, 
           SorbDistM, gTarPerFlow, Qm3m, bootN) %>% 
    mutate(stream = "UTLC")




SFT2 <-  SFT %>%  
  mutate(dateTime = as.POSIXct(dateTime, origin = "1970-01-01 00:00.00 UTC")) %>% 
  filter(dateTime >= StartStumpf & dateTime <= EndStumpf) %>% 
  mutate(DRPloadUS.gPwin.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadUS.gPwin, 0),
         Svol_gPwindow.hf = ifelse(gTarPerFlow == "Gtarget", Svol_gPwindow, 0),
         DRPloadDS.gPwin.wS.hf = ifelse(gTarPerFlow == "Gtarget", DRPloadDS.gPwin.wS, 0)) %>% 
  select(dateTime, DRPloadUS.gPwin, DRPloadUS.gPwin.hf, DRPloadDS.gPwin.wS, DRPloadDS.gPwin.wS.hf,
         Svol_gPwindow.hf, 
         SorbDistM, gTarPerFlow, Qm3m, bootN) %>% 
  mutate(stream = "SFT")

TribSorb <- rbind(WC2, UTLC2, SFT2) %>% 
  mutate(stream = as.factor(stream),
         bootN = as.factor(bootN),
         gTarPerFlow = as.factor(gTarPerFlow)) %>% 
  filter(!is.na(gTarPerFlow))

#####################
# Little functions for calculating probs
#####################
p <- c(0.025, 0.5, 0.975)

p_names <- map_chr(p, ~paste0(.x*100, "per"))

p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)


#####################
# only make calcuate for intervals with Q over target % flow
#####################
TribSorbS <- TribSorb %>% 
  group_by(dateTime, stream, gTarPerFlow) %>% 
  summarize_at(vars(DRPloadUS.gPwin:SorbDistM, Qm3m), funs(!!!p_funs)) %>% 
  mutate(DOY = as.numeric(as.character(strftime(dateTime, format = "%j")))) %>% 
  group_by(DOY, stream, gTarPerFlow) %>% 
  summarize_at(vars(DRPloadUS.gPwin_2.5per:Qm3m_97.5per), list("s" = sum, "m" = mean)) 

TribSorbS2 <- TribSorbS %>% 
  select(DOY:gTarPerFlow,
    DRPloadUS.gPwin_2.5per_s:Svol_gPwindow.hf_2.5per_s, 
    DRPloadUS.gPwin_50per_s:Svol_gPwindow.hf_50per_s, 
    DRPloadUS.gPwin_97.5per_s:Svol_gPwindow.hf_97.5per_s, 
    SorbDistM_2.5per_m, Qm3m_2.5per_m,
    SorbDistM_50per_m, Qm3m_50per_m,
    SorbDistM_97.5per_m, Qm3m_97.5per_m) %>% 
  mutate(SorbDistPerWat = ifelse(stream == "SFT", (SorbDistM_50per_m/38715.5)*100,
                            ifelse(stream == "WC", (SorbDistM_50per_m/56692.6)*100,
                              ifelse(stream == "UTLC", (SorbDistM_50per_m/114467.8)*100, as.numeric("NA")))),
         TribDist2wat = ifelse(stream == "SFT", (38715.5)/1000,
                                 ifelse(stream == "WC", (56692.6)/1000,
                                        ifelse(stream == "UTLC", (114467.8)/1000, as.numeric("NA")))),
         AtCapDist2Wat = TribDist2wat - SorbDistM_50per_m/1000,
         stream = fct_relevel(stream, c("UTLC", "WC", "SFT"))) %>% 
  filter(gTarPerFlow == "Gtarget") %>% 
  mutate(perDRPsorbed = Svol_gPwindow.hf_50per_s/DRPloadUS.gPwin.hf_50per_s*100,
         stream = fct_recode(stream, "STF" = "SFT"))



#####################
# Percent sorbed or something
#####################
TribSorbS2$s75per = TribSorbS2$DRPloadUS.gPwin.hf_50per_s*0.75
TribSorbS2$s50per = TribSorbS2$DRPloadUS.gPwin.hf_50per_s*0.5
TribSorbS2$s25per = TribSorbS2$DRPloadUS.gPwin.hf_50per_s*0.25
TribSorbS2$s10per = TribSorbS2$DRPloadUS.gPwin.hf_50per_s*0.1


#####################
# Percent sorbed at annual scale
#####################


# Portion sorbed during stumpf window - high flows

TribSorbSsw <- TribSorb %>% 
  group_by(stream, bootN) %>% #don't need gTarPerFlow because have that in windows
  summarize_at(vars(DRPloadUS.gPwin:Svol_gPwindow.hf), list("s" = sum)) %>% 
  mutate(DRPsorbed.hf = Svol_gPwindow.hf_s/DRPloadUS.gPwin.hf_s*100,
         DRPsorbed.allF = Svol_gPwindow.hf_s/DRPloadUS.gPwin_s*100,
         stream = fct_recode(stream, "STF" = "SFT")) 

TribSorbSsw %>% 
  group_by(stream) %>% 
  summarize_at(vars(DRPloadUS.gPwin_s:DRPsorbed.allF),  list(mean = mean), na.rm = TRUE)



Fig3cXaxis <- expression(atop("% DRP load sorbed","(Mar-June, >75% flows)"))

#####################
# FIGURES
#####################

# distance traveled
Fig3c <- ggplot() +
  geom_segment(data = TribSorbS2, aes(x = DOY, xend = DOY, y = TribDist2wat, yend = AtCapDist2Wat, color = stream), 
               arrow = arrow(length = unit(0.1,"cm")), size = 0.35) +
  stat_smooth(data = TribSorbS2, aes(x = DOY, y = AtCapDist2Wat, color = stream), span = 0.25, alpha = 0.5)  +
  coord_flip()  +
  scale_y_reverse(limits = c(120,0)) +
  theme_bw() +
  xlab("Day of year") +
  ylab("Distance to Waterville, OH (km)") +
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2")) +
  theme(legend.position = "none", #c(0.9,0.15)
        legend.text = element_text(size = 14),
        legend.title = element_text(size = 18),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank()) 
  
Fig3b <-  ggplot() +
  geom_line(data = TribSorbS2, aes(y = s50per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = TribSorbS2, aes(y = s25per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = TribSorbS2, aes(y = s10per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_point(data = TribSorbS2, aes(y = Svol_gPwindow.hf_50per_s/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, fill = stream, size = Qm3m_50per_m),
             shape = 21, alpha = 0.5)+
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.0001, 675),
                labels = function(x) sprintf("%g", x)) + #
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.0001, 675),
                labels = function(x) sprintf("%g", x)) +
  annotate("text", x = 600, y = 300*0.5, label = "50%", color = "grey40", size = 5) +
  annotate("text", x = 600, y = 300*0.25, label = "25%", color = "grey40", size = 5) +
  annotate("text", x = 600, y = 300*0.1, label = "10%", color = "grey40", size = 5) +
  scale_fill_manual(values = c("STF"="light salmon",
                               "UTLC"="medium sea green",
                               "WC" = "plum2"), name = "Stream") +
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2"), guide = FALSE) +
  scale_size(name = expression(paste("Discharge (", m^3," ", min^-1,")")))+
  xlab(expression(paste("DRP load (kg ",day^-1,")"))) +
  ylab(expression(paste("DRP sorption (kg ",day^-1,")"))) +
  theme_bw() +
  theme(legend.position = c(0.21,0.75),
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 12),
        legend.key.height = unit(0.5,"cm"),
        legend.spacing.y = unit(0.05,"cm"),
        legend.margin = margin(0.1,0,0,0, unit="cm"),
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  stat_smooth(data = TribSorbS2, aes(y = Svol_gPwindow.hf_50per_s/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, color = stream), method = "lm")

Fig3a <-  ggplot(TribSorbSsw)  +
  geom_density(aes(x = DRPsorbed.hf, fill = stream))+
  geom_rug(aes(x = DRPsorbed.hf, y = 0, color = stream), position = position_jitter(height = 0)) +
  xlab("% DRP sorbed") + #
  ylab("Density") +
  scale_fill_manual(values = c("STF"="light salmon",
                               "UTLC"="medium sea green",
                               "WC" = "plum2"))+
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2"))+
  xlim(25,50) +
  theme_bw() +
  theme(legend.position = "none",#c(0.9,0.9)
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank()) 


# https://cran.r-project.org/web/packages/egg/vignettes/Overview.html


Fig3a.g <- ggplotGrob(Fig3a)
Fig3b.g <- ggplotGrob(Fig3b)
Fig3c.g <- ggplotGrob(Fig3c)

Fig3a.gtf <- gtable_frame(Fig3a.g,  width = unit(0.3, "null"), height = unit(0.45, "null")) #draws cells: debug = TRUE, 
Fig3b.gtf <- gtable_frame(Fig3b.g, width = unit(0.7, "null"), height = unit(0.45, "null"))
Fig3c.gtf <- gtable_frame(Fig3c.g, width = unit(2, "null"), height = unit(0.2, "null"))


Fig3ab.gtf <- gtable_frame(gtable_cbind(Fig3a.gtf, Fig3b.gtf),
                           width = unit(0.45, "null"), height = unit(0.45, "null"))
Fig3.gtf <- gtable_rbind(Fig3ab.gtf, Fig3c.gtf)



pdf("05_Figures/07p_Fig3_take2.pdf", height = 7, width = 8)
grid.newpage()
grid.draw(Fig3.gtf)
grid.text("a", x = unit(0.02,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.41,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.02,"npc"), y = unit(0.35,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()



#####################
# SAVE/LOAD
#####################

save.image("03_Rdata/07_TribFigs_rdat")
# load("07_TribFigs_rdat")





