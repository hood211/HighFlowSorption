# JMH Mar 2020
# Generates Fig 3

library(tidyverse)
library(grid)
library(egg)

# get data ----
maum_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maum_50col.csv"), row.names = 1)
maum_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maum_0col.csv"), row.names = 1) %>% 
              select(Y, gPsorbedgDM_50per:perSorbP_50per)
maum <- maum_0col %>% 
        full_join(maum_50col, by = "Y", suffix = c("_0col", "_50col")) %>% 
        mutate(Yg = fct_relevel(Yg, c("L2003", "G2003")))


maumCyanos_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maumCyano_50col.csv"), row.names = 1) %>% 
                    mutate(col = "50per")

maumCyanos_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "07_Fig3Dat_maumCyano_0col.csv"), row.names = 1)%>% 
                    mutate(col = "0per")

maumCyanos <- rbind(maumCyanos_0col, maumCyanos_50col) %>% 
              filter(!(var == "Observed" & col == "50per")) %>% 
              mutate(var2 = paste0(var, "_", col),
                     var2 = case_when(var2 == "Observed_0per" ~ "Observed",
                                      var2 == "Without P sorption_0per" ~ "0",
                                      var2 == "Without P sorption_50per" ~ "0.5"),
                     var2 = fct_relevel(var2, "Observed", "0", "0.5"))

# for results ----
      summary(maum$perSorbP_50per_0col)
      summary(maum$perSorbP_50per_50col)
      summary(maum[maum$Yg == "G2003",]$gPsorbP_MjHf_50per_0col/1e6)
      summary(maum[maum$Yg == "G2003",]$gPsorbP_MjHf_50per_50col/1e6)
      
      # used in graphical table of contents
      summary(maum[maum$Yg == "G2003",]$gDRP_MjHf_50per/1e6)
      summary(maum[maum$Yg == "G2003",]$gPsorbP_MjHf_50per_0col/1e6)
      summary(maum[maum$Yg == "G2003",]$gPsorbP_MjHf_50per_50col/1e6)
      
# Figures ----
## Fig 3a ----
Fig3a <- ggplot() +
  geom_ribbon(data = maum, aes(ymin = gDRP_MjHf_50per/1e6  + gPsorbP_MjHf_50per_0col/1e6, 
                               ymax = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per_50col/1e6, x = Y), 
              fill = "gold2", alpha = 40/100) +
  # geom_segment(data = maum, aes(y = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per/1e6, yend =  gDRP_MjHf_50per/1e6+5, x = Y, xend=Y), 
  #              arrow = arrow(length = unit(0.1, "inches")), color = "black", size = 1.15)+
  geom_point(data = maum, aes(y = gDRP_MjHf_50per/1e6, x = Y), fill = "black", shape = 21, size = 3) +
  geom_line(data = maum, aes(y = gDRP_MjHf_50per/1e6, x = Y), fill = "grey70", shape = 21, size = 0.5) +
  geom_point(data = maum, aes(y = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per_50col/1e6, x = Y), fill = "gold2", shape = 22, size = 3) +
  geom_point(data = maum, aes(y = gDRP_MjHf_50per/1e6 + gPsorbP_MjHf_50per_0col/1e6, x = Y), fill = "gold2", shape = 21, size = 3) +
  ylab("DRP load (tons P)") +
  xlab("Year") +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black"),
        panel.background = element_rect(fill = "transparent")) +
  geom_point(aes(y =620, x = 1975), fill = "gold2", size = 5, shape = 21) +
  annotate("text", label = "Without P sorption (colloidal-P:DRP = 0)", y = 620, x = 1976, size = 7, hjust = 0) +
  geom_point(aes(y =580, x = 1975), fill = "gold2", size = 5, shape = 22) +
  annotate("text", label = "Without P sorption (colloidal-P:DRP = 0.5)", y = 580, x = 1976, size = 7, hjust = 0) +
  geom_point(aes(y =540, x = 1975), color = "black", size = 5) +
  annotate("text", label = "Observed", y = 540, x = 1976, size = 7, hjust = 0) 

## Fig 3b ----
Fig3b <-ggplot() +
  geom_segment(data = maum, aes(y = gPsorbP_MjHf_50per_50col/1e6,
                                yend = gPsorbP_MjHf_50per_0col/1e6,
                                x = gDRP_MjHf_50per/1e6, xend = gDRP_MjHf_50per/1e6,  fill = Yg)) +
  geom_point(data = maum, aes(y = gPsorbP_MjHf_50per_0col/1e6, x = gDRP_MjHf_50per/1e6,  fill = Yg), shape = 21, size = 4) +
  geom_point(data = maum, aes(y = gPsorbP_MjHf_50per_50col/1e6, x = gDRP_MjHf_50per/1e6,  color = Yg), shape = 22, size = 4, fill = "white") +
  geom_point(size = 4, shape = 21) +
  geom_segment(aes(x = 0, xend = 290, y = 0, yend = 290), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 560, y = 0, yend = 560*0.5), color = "grey", linetype = "dashed") +
  geom_segment(aes(x = 0, xend = 560, y = 0, yend = 560*0.25), color = "grey", linetype = "dashed") +
  annotate("text", x = 280, y = 300, label = "100%", color = "grey40", size = 6) +
  annotate("text", x = 600, y = 290, label = "50%", color = "grey40", size = 6) +
  annotate("text", x = 600, y = 145, label = "25%", color = "grey40", size = 6) +
  scale_fill_manual(values = c("dodgerblue", "firebrick"), name = NULL, labels = c("< 2003","≥ 2003")) +
  scale_color_manual(values = c("dodgerblue", "firebrick"), name = NULL, labels = c("< 2003","≥ 2003")) +
  theme_bw() +
  xlab("DRP load (tons P)") +
  ylab("P sorption (tons P)") +
  ylim(0,300) +
  xlim(0,600) +
  theme(
    legend.position = "none",
    axis.title = element_text(size = 24),
    axis.text = element_text(size = 16),
    panel.border = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_blank(),
    axis.line.x = element_line(size = 1, color = "black"),
    axis.line.y = element_line(size = 1, color = "black")) +
  annotate("text", label = "Colloidal-P:DRP = 0", y = 280, x = 10, size = 8, hjust = 0, fontface = "bold") +
  geom_point(aes(y =260, x = 20), fill = "dodgerblue", size = 5, shape = 21) +
  annotate("text", label = "≤ 2003", y = 260, x = 35, size = 7, hjust = 0) +
  geom_point(aes(y =240, x = 20), fill = "firebrick", size = 5, shape = 21) +
  annotate("text", label = "≥ 2003", y = 240, x = 35, size = 7, hjust = 0) +
  annotate("text", label = "Colloidal-P:DRP = 0.5", y = 220, x = 10, size = 8, hjust = 0, fontface = "bold") +
  geom_point(aes(y =200, x = 20), color = "dodgerblue", size = 5, shape = 22) +
  annotate("text", label = "≤ 2003", y = 200, x = 35, size = 7, hjust = 0) +
  geom_point(aes(y =180, x = 20), color = "firebrick", size = 5, shape = 22) +
  annotate("text", label = "≥ 2003", y = 180, x = 35, size = 7, hjust = 0) 


## Fig 3c ----
# estimate of cyanos with and without sorption

Fig3c <-  ggplot(maumCyanos,
                 aes(y = values, x = var2, label = Yn2))+
  geom_boxplot(fill = "grey90", alpha = 0.7, outlier.color = "transparent") + #"#7fbf7b"
  geom_text(position = "jitter", size = 6) +
  # facet_grid(var~.) +
  # ylim(0,30) +
  theme_bw() +
        ylab(expression(atop("Cyanobacteria bloom", paste("(kiloton dry mass)")))) +
  xlab(NULL) +
  theme(axis.title = element_text(size = 24),
        axis.text.y = element_text(size = 16),
        axis.text.x = element_text(size = 20),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  annotate("text", x = 1, y = 125, label = "Observed", size = 7, fontface = "bold")+
  annotate("text", x = 2.5, y = 125, 
           label = "Without P sorption", size = 7, fontface = "bold") +
  annotate("text", x = 2.5, y = 115, 
           label = "(Colloidal-P:DRP)", size = 7, fontface = "bold") +
  geom_vline(xintercept = 1.5, color = "grey", size = 2)


# Print ----
Fig3a.g <- ggplotGrob(Fig3a)
Fig3b.g <- ggplotGrob(Fig3b)
Fig3c.g <- ggplotGrob(Fig3c)



Fig3a.gtf <- gtable_frame(Fig3a.g,  width = unit(1.8, "null"), height = unit(0.9, "null"))
Fig3b.gtf <- gtable_frame(Fig3b.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))
Fig3c.gtf <- gtable_frame(Fig3c.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))

Fig3bC.gtf <- gtable_frame(gtable_cbind(Fig3b.gtf, Fig3c.gtf),  width = unit(1.8, "null"), height = unit(0.9, "null"))

Fig3.gtf <- gtable_rbind(Fig3a.gtf, Fig3bC.gtf)




png("05_Figures/07_Fig3.png", units = "in", height = 12, width = 14, res = 300)
grid.newpage()
grid.draw(Fig3.gtf)
grid.text("a", x = unit(0.01,"npc"), y = unit(0.98,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
grid.text("b", x = unit(0.01,"npc"), y = unit(0.50,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
grid.text("c", x = unit(0.51,"npc"), y = unit(0.50,"npc"), gp=gpar(fontsize = 30, fontface = "bold"))
dev.off()


## Graphical abstract ----

FigGraphAb <- ggplot(maumCyanos,
       aes(y = values, x = var2, label = Yn2, fill = var2))+
  geom_boxplot(outlier.color = "transparent") + #"#7fbf7b"
  scale_fill_manual(values = c("dodgerblue1", "#7fbf7b", "green")) +
  ylim(0,45) +
  theme_bw() +
  ylab(NULL) +
  xlab(NULL) +
  theme(legend.position = "none",
        axis.title = element_text(size = 24, face = "bold"),
        axis.text.y = element_text(size = 24),
        axis.text.x = element_blank(),
        panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  annotate("text", x = 1, y = 40, label = "Observed", size = 12, fontface = "bold")+
  annotate("text", x = 2.5, y = 40, 
           label = "Without P sorption", size = 12, fontface = "bold") +
  geom_vline(xintercept = 1.5, color = "grey", size = 2)

png("05_Figures/07_GraphicalAbstractCyanoFig.png", units = "in", height =4, width = 8.5, res = 300)
FigGraphAb
dev.off()


# save ----
# save.image(file.path(here::here("03_Rdata"), "07c_Fig3_comp_Rdat"))
load(file.path(here::here("03_Rdata"), "07c_Fig3_comp_Rdat"))
