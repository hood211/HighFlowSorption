# Fall 2020 JMH
# loads and cleans slope data which was measured by WK in QGIS

library(tidyverse)

# load data ----


slopes <-  read.csv("01_RawData/SlopesClean20201015jmh.csv") %>% 
  mutate(Gauge_start = as.factor(Gauge_start),
         Gauge_end = as.factor(Gauge_end),
         Gauge_start = fct_recode(Gauge_start, AGDF = "Auglaize_defiance",
                                  MMDF = "Maumee_defiance",
                                  BRFL = "Blanchard_above_Findlay",
                                  BRDP = "Blanchard_at_Dupont",
                                  BROT = "Blanchard_at_Ottawa",
                                  BRHB = "Blanchard_heidelberg",
                                  BRBR = "Blanchard_near_Mt_Blanchard",
                                  BRGB = "Blancharg_at_Gilboa"),
         Gauge_end = fct_recode(Gauge_end, MMWV = "Maumee_waterville",
                                AGDF = "Auglaize_defiance",
                                  MMDF = "Maumee_defiance ",
                                  BRFL = "Blanchard_above_Findlay",
                                  BRDP = "Blanchard_at_Dupont",
                                  BROT = "Blanchard_at_Ottawa",
                                  BRHB = "Blanchard_heidelberg",
                                  BRBR = "Blanchard_near_Mt_Blanchard",
                                  BRGB = "Blanchard_at_Gilboa"),
         Seg = paste0(Gauge_start,"___", Gauge_end)) %>% 
  select(-Stream)

# Export data ----
write.csv(slopes, "04a_generatedDataOnGit/03d_slopes.csv")
