# JMH, Mar 2022
# Generates Fig 2

# Libraries ----
library(tidyverse)
library(egg)
library(grid)

# data ----

## for fig 2a ----
Fig2aDat_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "05_TribDat4Fig2a_0col.csv"), row.names = 1) %>% 
                  mutate(colloidRatio = "per0")
Fig2aDat_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "05b_TribDat4Fig2a_50col.csv"), row.names = 1) %>% 
                  mutate(colloidRatio = "per50")

Fig2aDat <- rbind(Fig2aDat_0col, Fig2aDat_50col)

## for other figs 2b ----
Fig2bcDat_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "05_TribDat4Fig2bc_0col.csv"), row.names = 1)%>% 
                  mutate(colloidRatio = "per0")
Fig2bcDat_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "05b_TribDat4Fig2bc_50col.csv"), row.names = 1)%>% 
                  mutate(colloidRatio = "per50")

Fig2bcDat <- rbind(Fig2bcDat_0col, Fig2bcDat_50col)

# Figures ----

## Fig 2a % DRP sorbed ----
Fig2a <- ggplot(Fig2aDat)  +
  geom_density(aes(x = DRPsorbed.hf, fill = stream, linetype = colloidRatio)) +
  # geom_rug(aes(x = DRPsorbed.hf, y = 0, color = stream), position = position_jitter(height = 0)) +
  xlab("% DRP sorbed") + #
  ylab("Density") +
  scale_fill_manual(values = c("STF"="light salmon",
                               "UTLC"="medium sea green",
                               "WC" = "plum2"))+
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2"))+
  xlim(0,60) +
  theme_bw() +
  theme(legend.position = "none",#c(0.9,0.9)
        legend.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 18),
        axis.text = element_text(size = 12),
        plot.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank()) 


## Fig 2b P sorp v. DRP load ----
Fig2bMain <- ggplot() +
  geom_line(data = Fig2bcDat_0col, aes(y = s50per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = Fig2bcDat_0col, aes(y = s25per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = Fig2bcDat_0col, aes(y = s10per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_point(data = Fig2bcDat_0col, aes(y = (Svol_gPwindow.hf_50per_s+0.075)/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, fill = stream, size = Qm3m_50per_m),
             shape = 21, alpha = 0.5) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.0005, 250),
                labels = function(x) sprintf("%g", x)) + #
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.00001, 250),
                labels = function(x) sprintf("%g", x)) +
  annotate("text", x = 0.0008, y = 0.002*0.6, label = "50%", color = "grey40", size = 5) +
  annotate("text", x = 0.0008, y = 0.002*0.25, label = "25%", color = "grey40", size = 5) +
  annotate("text", x = 0.0008, y = 0.002*0.1, label = "10%", color = "grey40", size = 5) +
  scale_fill_manual(values = c("STF"="light salmon",
                               "UTLC"="medium sea green",
                               "WC" = "plum2"), name = "Stream") +
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2"), guide = FALSE) +
  scale_size(name = expression(paste("Discharge (", m^3," ", min^-1,")")))+
  xlab(expression(paste("DRP load (kg ",day^-1,")"))) +
  ylab(expression(paste("P sorption (kg ",day^-1,")"))) +
  theme_bw() +
  theme(legend.position = c(0.21,0.70),
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
  stat_smooth(data = Fig2bcDat_0col, aes(y = (Svol_gPwindow.hf_50per_s+0.075)/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, color = stream), method = "lm")+
  annotate(geom = "text", x = 10, y = 0.055, label = "colloidal-P:DRP = 0.5", fontface = "bold")


Fig2bInset <- ggplot() +
  geom_line(data = Fig2bcDat_50col, aes(y = s50per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = Fig2bcDat_50col, aes(y = s25per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_line(data = Fig2bcDat_50col, aes(y = s10per/1000, x = DRPloadUS.gPwin.hf_50per_s/1000), color = "grey", linetype = "dashed")+
  geom_point(data = Fig2bcDat_50col, aes(y = (Svol_gPwindow.hf_50per_s+0.075)/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, fill = stream, size = Qm3m_50per_m),
             shape = 21, alpha = 0.5) +
  scale_x_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.0005, 250),
                labels = function(x) sprintf("%g", x)) + #
  scale_y_log10(breaks = c(0.0001, 0.001, 0.01, 0.1,1, 10,100), limits = c(0.00001, 250),
                labels = function(x) sprintf("%g", x)) +
  scale_fill_manual(values = c("STF"="light salmon",
                               "UTLC"="medium sea green",
                               "WC" = "plum2"), name = "Stream") +
  scale_color_manual(values = c("STF"="light salmon",
                                "UTLC"="medium sea green",
                                "WC" = "plum2"), guide = FALSE) +
  scale_size(name = expression(paste("Discharge (", m^3," ", min^-1,")")))+
  xlab(NULL) +
  ylab(NULL) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title = element_text(size = 18),
        # axis.text = element_text(size = 12),
        axis.text = element_blank(),
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "transparent", color = "transparent"),
        panel.grid.minor = element_blank()) +
  guides(fill = guide_legend(override.aes = list(size=5))) +
  stat_smooth(data = Fig2bcDat_50col, aes(y = (Svol_gPwindow.hf_50per_s+0.075)/1000, x = DRPloadUS.gPwin.hf_50per_s/1000, color = stream), method = "lm")


Fig2b <- Fig2bMain + annotation_custom(
  ggplotGrob(Fig2bInset),
  xmin = log10(0.19), xmax  = log10(500), ymin = log10(0.000004), ymax = log10(0.04))

## Fig 2c distance traveled ----
Fig2c <- ggplot() +
  geom_segment(data = Fig2bcDat_0col, aes(x = DOY, xend = DOY, y = TribDist2wat_km, yend = AtCapDist2Wat, color = stream), 
               arrow = arrow(length = unit(0.1,"cm")), size = 0.35) +
  geom_segment(data = Fig2bcDat_50col, aes(x = DOY, xend = DOY, y = TribDist2wat_km, yend = AtCapDist2Wat, color = stream), 
               size = 1) +
  stat_smooth(data = Fig2bcDat_0col, aes(x = DOY, y = AtCapDist2Wat, color = stream), span = 0.25, alpha = 0.5)  +
  coord_flip()  +
  scale_y_reverse(limits = c(120,0)) +
  theme_bw() +
  xlim(65,205)+
  ylab("Distance to Waterville, OH (km)") +
  xlab("Day of year") +
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
        panel.grid.minor = element_blank()) +
  annotate("text", x = 205, y = 120, label = "km traveled before P sorption saturates",  size = 5, hjust = 0, fontface = "bold") +#
  
  annotate("text", y = 108, x = 190, label = "UTLC", color = "medium sea green",  size = 5, hjust = 0) +
  geom_segment(aes(y = 120, yend = 110, x = 190, xend = 190), color = "medium sea green", arrow = arrow(length = unit(0.02, "npc")), size = 1) +
  
  annotate("text", y = 83, x = 190, label = "WC", color = "plum2",  size = 5, hjust = 0) +
  geom_segment(aes(y = 95, yend = 85, x = 190, xend = 190), color = "plum2", arrow = arrow(length = unit(0.02, "npc")), size = 1) +
  
  annotate("text", y = 60, x = 190, label = "STF", color = "light salmon",  size = 5, hjust = 0) +
  geom_segment(aes(y = 72, yend = 62, x = 190, xend = 190), color = "light salmon", arrow = arrow(length = unit(0.02, "npc")), size = 1) 





Fig2a.g <- ggplotGrob(Fig2a)
Fig2b.g <- ggplotGrob(Fig2b)
Fig2c.g <- ggplotGrob(Fig2c)


Fig2a.gtf <- gtable_frame(Fig2a.g,  width = unit(0.3, "null"), height = unit(0.35, "null")) #draws cells: debug = TRUE, 
Fig2b.gtf <- gtable_frame(Fig2b.g, width = unit(0.7, "null"), height = unit(0.35, "null"))
Fig2c.gtf <- gtable_frame(Fig2c.g, width = unit(2, "null"), height = unit(0.3, "null"))


Fig2ab.gtf <- gtable_frame(gtable_cbind(Fig2a.gtf, Fig2b.gtf),
                           width = unit(0.45, "null"), height = unit(0.45, "null"))
Fig2.gtf <- gtable_rbind(Fig2ab.gtf, Fig2c.gtf)


## Print fig ----
png("05_Figures/05_Fig2_combined.png", units="in", width=8, height=7, res=300)
grid.newpage()
grid.draw(Fig2.gtf)
grid.text("a", x = unit(0.02,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("b", x = unit(0.41,"npc"), y = unit(0.97,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
grid.text("c", x = unit(0.02,"npc"), y = unit(0.39,"npc"), gp=gpar(fontsize = 25, fontface = "bold"))
dev.off()
