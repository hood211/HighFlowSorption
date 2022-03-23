# JMH, Feb 2022

# JMH Feb 2022
# Generates Fig. S1

# Libraries ----
library(tidyverse)
library(grid)
library(egg)


# Get data ----
wqMaum <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "02d_MaumeeWatervilleWaterQual.csv"), row.names = 1) %>% 
            mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d"),
                   HighFlow = ifelse(HighFlow == "HighFlow", "G75th", "L75th"),
                   # Y = as.character(strftime(dateTime, format = "%Y")),
                   Mc = as.factor(as.character(strftime(Date, format = "%m"))),
                   # create month groups: Mar-June, Jul, Aug-Feb
                   # TargetMonths = as.factor(ifelse(Mc %in% c("03", "04", "05", "06"), "MarJun",
                   #                                 ifelse(Mc == "07", "Jul","AugFeb")))) %>% 
                   TargetMonths = as.factor(ifelse(Mc %in% c("03", "04", "05", "06"), "MarJun",
                                                   ifelse(Mc %in% c("07", "08", "09"), "JulSep","OctFeb")))) %>% 
                  # remove 1978, 79, 80,81; could leave 78 for spring stuff
                  filter(Y != 1978 & Y != 1979 & Y != 1980 & Y != 1981) %>% 
                  select(Date, Qm3day:gTPday, HighFlow, Y, M, TargetMonths) %>% 
                  # sum to month group for low and high flow
                  group_by(Y, TargetMonths, HighFlow) %>% 
                  summarise(Qm3window = sum(Qm3day, na.rm = T),
                            gDRPwindow= sum(gSRPday, na.rm = TRUE),
                            gSSwindow = sum(gSSday, na.rm = TRUE),
                            gTPwindow = sum(gTPday, na.rm = TRUE))  %>% 
                  # easy codes for ploting
                  mutate(key = as.factor(paste0(TargetMonths,"_",HighFlow))) %>% 
                  # mutate(key2 = fct_rev(fct_relevel(key, c("AugFeb_L75th","MarJun_L75th","Jul_L75th","AugFeb_G75th", "MarJun_G75th", "Jul_G75th")))) %>% 
                  mutate(key2 = fct_rev(fct_relevel(key, c("OctFeb_L75th","MarJun_L75th","JulSep_L75th","OctFeb_G75th", "MarJun_G75th", "JulSep_G75th")))) %>% 
                  mutate(DRP_TP = gDRPwindow/gTPwindow)

# for graphical abstract
  summary(wqMaum[wqMaum$Y >= 2003 & wqMaum$TargetMonths == "MarJun" & wqMaum$HighFlow == "G75th",]$gTPwindow/1e6)
  
  # line 466
  # area of maumee watershed = 21,540 km2
  222*0.6/21540*1000
  255*0.6/21540*1000
  303*0.6/21540*1000

# FIG S1 ----

FigS1a <-   ggplot(wqMaum, 
                   aes(y = gDRPwindow/1e6, x = Y, fill = key2)) +
                  geom_bar(stat = "identity") +
                  scale_fill_manual(values = c("pink", "firebrick1", "firebrick4","lightblue", "steelblue3", "steelblue4"), name = NULL,
                                    labels = c(">75% flow: Jul-Sep",
                                               ">75% flow: Mar-Jun",
                                               ">75% flow: Oct-Feb",
                                               "<75% flow: Jul-Sep",
                                               "<75% flow: Mar-Jun",
                                               "<75% flow: Oct-Feb"))+
                  theme_bw() +
                  ylab(expression(atop("DRP load", "(Metric ton P)"))) +
                  xlab("Year") +
                  ylim(0,1000) +
                  theme(legend.position = "top", #c(0.88,0.83)
                        legend.text = element_text(size = 15),
                        legend.title = element_text(size = 12, face = "bold"),
                        legend.key.height = unit(0.45,"cm"),
                        legend.spacing.y = unit(0.05,"cm"),
                        legend.margin = margin(0.1,0,0,0, unit="cm"),
                        legend.background = element_rect(fill = "transparent"),
                        axis.title = element_text(size = 18),
                        axis.text = element_text(size = 12),
                        panel.background = element_rect(fill = "transparent"),
                        panel.grid.minor = element_blank()) 


FigS1b <-   ggplot(wqMaum, 
                   aes(y = gTPwindow/1e6, x = Y, fill = key2)) +
                  geom_bar(stat = "identity")+
                  scale_fill_manual(values = c("pink", "firebrick1", "firebrick4","lightblue", "steelblue3", "steelblue4"), name = NULL,
                                    labels = c(">75% flow: Jul-Sep",
                                               ">75% flow: Mar-Jun",
                                               ">75% flow: Oct-Feb",
                                               "<75% flow: Jul-Sep",
                                               "<75% flow: Mar-Jun",
                                               "<75% flow: Oct-Feb"))+
                  theme_bw() +
                  ylab(expression(atop("TP load", "(Metric ton P)"))) +
                  xlab("Year") +
                  ylim(0,5000) +
                  theme(legend.position = "none",
                        legend.text = element_text(size = 10),
                        legend.title = element_text(size = 12, face = "bold"),
                        legend.key.height = unit(0.5,"cm"),
                        legend.spacing.y = unit(0.05,"cm"),
                        legend.margin = margin(0.1,0,0,0, unit="cm"),
                        legend.background = element_rect(fill = "transparent"),
                        axis.title = element_text(size = 18),
                        axis.text = element_text(size = 12),
                        panel.background = element_rect(fill = "transparent"),
                        panel.grid.minor = element_blank()) 


FigS1c <-   ggplot(wqMaum, 
                   aes(y = gSSwindow/1e12, x = Y, fill = key2)) +
                    geom_bar(stat = "identity")+
                    scale_fill_manual(values = c("pink", "firebrick1", "firebrick4","lightblue", "steelblue3", "steelblue4"), name = NULL,
                                      labels = c(">75% flow: Jul-Sep",
                                                 ">75% flow: Mar-Jun",
                                                 ">75% flow: Oct-Feb",
                                                 "<75% flow: Jul-Sep",
                                                 "<75% flow: Mar-Jun",
                                                 "<75% flow: Oct-Feb"))+
                    theme_bw() +
                    ylab(expression(atop("Suspended sediment load", paste("(mT dry mass x ",10^6,")")))) +
                    xlab("Year") +
                    theme(legend.position = "none",
                          legend.text = element_text(size = 10),
                          legend.title = element_text(size = 12, face = "bold"),
                          legend.key.height = unit(0.5,"cm"),
                          legend.spacing.y = unit(0.05,"cm"),
                          legend.margin = margin(0.1,0,0,0, unit="cm"),
                          legend.background = element_rect(fill = "transparent"),
                          axis.title = element_text(size = 18),
                          axis.text = element_text(size = 12),
                          panel.background = element_rect(fill = "transparent"),
                          panel.grid.minor = element_blank()) 

FigS1a.g <- ggplotGrob(FigS1a)
FigS1b.g <- ggplotGrob(FigS1b)
FigS1c.g <- ggplotGrob(FigS1c)

FigS1a.gtf <- gtable_frame(FigS1a.g, width = unit(1, "null"), height = unit(0.75, "null"))
FigS1b.gtf <- gtable_frame(FigS1b.g, width = unit(1, "null"), height = unit(0.75, "null"))
FigS1c.gtf <- gtable_frame(FigS1c.g, width = unit(1, "null"), height = unit(0.75, "null"))

FigS1.gtf <- gtable_frame(gtable_rbind(FigS1a.gtf, FigS1b.gtf, FigS1c.gtf),
                          width = unit(1, "null"), height = unit(1.5, "null"))

png("05_Figures/08_FigS1.png", units="in", width=8, height=11, res=300)
grid.newpage()
grid.draw(FigS1.gtf)
grid.text("a", x = unit(0.02,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.02,"npc"), y = unit(0.66,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.02,"npc"), y = unit(0.33, "npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()



                  

# save ----
# save.image(file.path(here::here("03_Rdata"),"08_MaumeeWatervillePloadsFigS1_Rdata"))
# load(file.path(here::here("03_Rdata"),"08_MaumeeWatervillePloadsFigS1_Rdata"))
