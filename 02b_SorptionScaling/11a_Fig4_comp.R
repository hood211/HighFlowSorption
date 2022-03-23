# JMH Feb 2022
# Generates Fig 4

# Libraries ----
library(tidyverse)
library(grid)
library(egg)

# Data ---- 
G03_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_G03.csv"))
L03_0col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_L03.csv"))

G03_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_50perDRP_G03.csv"))
L03_50col <- read.csv(file.path(here::here("04a_generatedDataOnGit"), "11a_Fig4dat_50perDRP_L03.csv"))

G03_0col_2 <-  G03_0col %>% 
          left_join(G03_50col %>% 
                      select(Y, mtDRP_MjHf_pred75ss, mtDRP_MjHf_pred75ss_se), by = "Y", suffix = c("_0col", "_50col"))

# Fig 4 ----
Fig4 <- ggplot() +
  geom_line(data = L03_0col, aes(y = mtDRP_MjHf_pred, x = Qm3e6_MjHf), color = "dodgerblue", size = 1.3) +
  geom_ribbon(data = L03_0col, aes(ymin = mtDRP_MjHf_pred - mtDRP_MjHf_pred_se,
                                   ymax = mtDRP_MjHf_pred + mtDRP_MjHf_pred_se, 
                                   x = Qm3e6_MjHf), fill = "dodgerblue", alpha = 0.25) +
  geom_line(data = G03_0col_2, aes(y = mtDRP_MjHf_pred, x = Qm3e6_MjHf), color = "firebrick", size = 1.3) +
  geom_ribbon(data = G03_0col_2, aes(ymin = mtDRP_MjHf_pred - mtDRP_MjHf_pred_se,
                                     ymax = mtDRP_MjHf_pred + mtDRP_MjHf_pred_se, 
                                     x = Qm3e6_MjHf), fill = "firebrick", alpha = 0.25) +
  geom_point(data = L03_0col, aes(y = mtDRP_MjHf, x = Qm3e6_MjHf), fill = "dodgerblue", shape = 21, size = 4) +
  geom_point(data = G03_0col, aes(y = mtDRP_MjHf, x = Qm3e6_MjHf), fill = "firebrick", shape = 21, size = 4)+
  geom_ribbon(data = G03_0col_2, aes(ymin = mtDRP_MjHf_pred75ss_0col,
                                     ymax = mtDRP_MjHf_pred75ss_50col, x = Qm3e6_MjHf), 
              fill = "salmon", size = 1, alpha = 0.75, color = "black", linetype = "dashed") +
  scale_y_continuous("DRP load (tons P)") +
  xlab(expression(paste("Discharge (",m^3, "x ", 10^6,")"))) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_line(color = "grey90"),
        panel.background = element_rect(fill = "transparent"),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 16),
        panel.border = element_blank(),
        axis.line.x = element_line(size = 1, color = "black"),
        axis.line.y = element_line(size = 1, color = "black")) +
  geom_point(aes(x = 400, y = 400), shape = 21, fill = "firebrick", size = 3) +
  geom_point(aes(x = 400, y = 375), shape = 21, fill = "salmon", size = 3) +
  geom_point(aes(x = 400, y = 350), shape = 21, fill = "dodgerblue", size = 3) +
  annotate("text", x = 500, y = 400, label = "≥ 2003 observed", size = 7, color = "black", hjust = 0) +
  annotate("text", x = 500, y = 375, label = "≥ 2003 historic SS loads", size = 7, color = "black", hjust = 0) +
  annotate("text", x = 500, y = 350, label = "< 2003 observed", size = 7, color = "black", hjust = 0) 

Fig4.g <- ggplotGrob(Fig4)
Fig4.gtf <- gtable_frame(Fig4.g,  width = unit(0.8, "null"), height = unit(0.9, "null"))

png("05_Figures/11a_Fig4.png", units = "in", height = 6, width = 7, res = 300)
grid.newpage()
grid.draw(Fig4.gtf)
dev.off()


# save ----
# save.image(file.path(here::here("03_Rdata"), "11a_Fig4_comp_Rdat"))
