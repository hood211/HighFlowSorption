# JMH Jan 2022
# Estimates cyanoHABs with and without P sorption. 
# Assumes colloid-P:DRP = 0.5

# Libraries ----
library(tidyverse)

# data ----
## this is the raw bootstrap output for Maumee
maum <- read.csv(file.path(here::here("04a_generatedDataTooBigForGit"), "04b_RawBootstrapResults_Maumee_50perDRP.csv"))[-1,-1]

## Need this for TP data
wqMaum <- read.csv("04a_generatedDataOnGit/02d_MaumeeWatervilleWaterQual.csv", row.names = 1) %>% 
  mutate(Date = as.POSIXct(Date, format = "%Y-%m-%d", tz = "America/New_York"),
         SmpTimeWindowDay = 1)   #daily data so window in days -=1

wqMaum_MjHf <- wqMaum %>% 
  # filter Mar-July and high flow
  filter(M >= 3 & M < 7) %>% 
  filter(HighFlow == "HighFlow") %>% 
  # remove concentrations & SS/SRP loads (exactly the same)
  select(-c(Date,SS_gm3:TP_gm3, gSRPday, gSSday, HighFlow, M:SmpTimeWindowDay)) %>% 
  group_by(Y) %>% 
  summarize_at(vars(gTPday), sum) %>% 
  mutate(Y = as.factor(Y)) %>% 
  rename(gTP_MjHf = gTPday)

wqMaum_MjLowf <- wqMaum %>% 
  # filter Mar-July and high flow
  filter(M >= 3 & M < 7) %>% 
  filter(HighFlow == "LowFlow") %>% 
  # remove concentrations & SS/SRP loads (exactly the same)
  select(-c(Date,SS_gm3:TP_gm3, gSSday, HighFlow, M:SmpTimeWindowDay)) %>% 
  group_by(Y) %>% 
  summarize_at(vars(gTPday, gSRPday), sum) %>% 
  mutate(Y = as.factor(Y)) %>% 
  rename(gTP_MjLowf = gTPday,
         gSRP_MjLowf = gSRPday)



## Little functions for calculating probs ----
p <- c(0.025, 0.5, 0.975)
p_names <- map_chr(p, ~paste0(.x*100, "per"))
p_funs <- map(p, ~partial(quantile, probs = .x, na.rm = TRUE)) %>% 
  set_names(nm = p_names)

# Probs for Mar-June high flows ----
maumMjHf <- maum %>% 
  mutate(Date = as.POSIXct(Date, origin = "1970-01-01 00:00.00 UTC"),
         M = as.numeric(as.character(strftime(Date, format = "%m"))),
         Y = as.factor(as.character(strftime(Date, format = "%Y")))) %>% 
  # filter Mar-July and high flow
  filter(M >= 3 & M < 7) %>% 
  filter(HighFlow == "HighFlow") %>% 
  # remove concentrations
  select(-c(SS_gm3, SRP_gm3, HighFlow)) %>% 
  # sum for year
  group_by(Y, bootN) %>% 
  summarize_at(vars(Qm3day:gDRPwindowWOsorp), list(sum), na.rm = TRUE) %>% 
  #MjHf = Mar-June for high flow
  rename(Qm3_MjHf = Qm3day,
         gDRP_MjHf = gSRPday,
         gSS_MjHf = gSSday,
         gPsorbP_MjHf = gPsorbWindowP,
         gDRPwoSorp_MjHf = gDRPwindowWOsorpP) %>% 
  mutate(perSorbP = gPsorbP_MjHf/gDRP_MjHf) %>%  # potential doesn't use constraint on DRP, use this
  group_by(Y) %>% 
  summarize_at(vars(Qm3_MjHf:perSorbP), funs(!!!p_funs)) %>% 
  mutate(Yn = as.numeric(as.character(Y)))

maumMjHf2 <- maumMjHf %>% 
  left_join(wqMaum_MjHf, by = "Y") %>% 
  left_join(wqMaum_MjLowf, by = "Y") 

# Calc cyano index ----
# only doing this for 2002 and later
# this will calc. cyanos with all loads not just Highflow
mtDRP_MjHf_50per

wqMaumCyanos <- maumMjHf2 %>%
  # convert g to metric tons
  mutate(mtDRP_MjAllF = (gDRP_MjHf_50per + gSRP_MjLowf)/1e6, # DRP for all flows
         mtTP_MjAllF = (gTP_MjHf + gTP_MjLowf)/1e6, # TP for all flows
         mtPsorbP_MjHf_50per  = gPsorbP_MjHf_50per/1e6, # Sorption, only high flows
         # DRP without sorption - sorbed P goes to DRP
         mtDRPwoSorpP_MjAllf = mtDRP_MjAllF + mtPsorbP_MjHf_50per,
         # PP OBSERVED
         mtPP_MjAllF = mtTP_MjAllF - mtDRP_MjAllF,
         # PP WITHOUT SORPTION - SORBED P TAKEN FROM PP
         mtPPwoSorpP_MjAllF = mtPP_MjAllF - mtPsorbP_MjHf_50per,
         # bioavailable P - using each papers eq s
         # Their bioavailable P models produce very different estimates! Bertani is much higher
         # bioavailable P, from Stumpf et al. 2016
         mtBioavaP_MjAllf_S = mtDRP_MjAllF + (1-0.7)*(0.26 * mtPP_MjAllF),
         mtBioavaP_woSorp_MjAllf_S = mtDRPwoSorpP_MjAllf + (1-0.7)*(0.26 * (mtPPwoSorpP_MjAllF)),
         # Bertani
         mtBioavaP_MjAllf_B = mtDRP_MjAllF + 0.63 * mtPP_MjAllF, #0.63 is the portion of PP that is bioavailable
         mtBioavaP_woSorp_MjAllf_B = mtDRPwoSorpP_MjAllf + 0.63 * mtPPwoSorpP_MjAllF,
         # from Bertani et al. 2016 (Table 2)
         # not bioP is a monthly "weighted" average, but works out to dividing by months
         # Bb + B0 + Bw +Bt
         Bert_Year = Yn - 2007, # see eq. 4 and 5 in Bertani et al. 2016
         cyanos_Bert_obs = ifelse((7.91 + -18.9 + (mtBioavaP_MjAllf_B/4)* 111.8 + 3.57 * Bert_Year) < 0,
                                  7.91,
                                  7.91 + -18.9 + (mtBioavaP_MjAllf_B/4)* 111.8 + 3.57 * Bert_Year)/1000, # bloom size (1000 MT)
         cyanos_Bert_woSorp = ifelse((7.91 + -18.9 + (mtBioavaP_woSorp_MjAllf_B/4)* 111.8 + 3.57 * Bert_Year) < 0,
                                     7.91,
                                     7.91 + -18.9 + (mtBioavaP_woSorp_MjAllf_B/4)* 111.8 + 3.57 * Bert_Year)/1000, # bloom size (1000 MT)
         # from Stumpf et al. 2016
         # Table 2, Mar-Jun column
         cyano_Stump_obs = 0.47 * 10^(mtBioavaP_MjAllf_S * 0.00306) *4, # 4 is a conversion from cells to biomass based on Obenour
         cyano_Stump_woSorp = 0.47 * 10^(mtBioavaP_woSorp_MjAllf_S * 0.00306)*4,
         cyano_modelMean_obs = (cyanos_Bert_obs + cyano_Stump_obs)/2,
         cyano_modelMean_woSorp = (cyanos_Bert_woSorp + cyano_Stump_woSorp)/2) %>%  # without constraint on DRP, use this
  # after 2002
  filter(Yn >= 2002) %>%
  select(Yn, mtDRP_MjAllF:cyano_modelMean_woSorp) %>%
  pivot_longer(cols = c(mtDRP_MjAllF:cyano_modelMean_woSorp), names_to = "var", values_to = "values")

ggplot(wqMaumCyanos %>% 
         filter(var %in% c("cyano_modelMean_obs", "cyano_modelMean_woSorp")), aes(y = values, x = Yn, color = var)) +
  geom_boxplot()

wqMaumCyanosW <- wqMaumCyanos %>% 
  filter(var %in% c("cyano_modelMean_obs", "cyano_modelMean_woSorp")) %>% 
  pivot_wider(id_col = Yn, names_from = var, values_from = values) %>% 
  mutate(CyanosPerHigher = cyano_modelMean_woSorp/cyano_modelMean_obs)

summary(wqMaumCyanosW)

  

# Export data ----
# write.csv(maumMjHf, file.path(here::here("04a_generatedDataOnGit"), "06b_Scale2WVall_MjHf_50perDRP.csv"))
# write.csv(wqMaumCyanos, file.path(here::here("04a_generatedDataOnGit"), "06b_Scale2WVall_MjHf_cyanos_50perDRP.csv"))

# Save image
# save.image(file.path(here::here("03_Rdata"), "06b_MaumeeBootMungePcyanos_50perDRP_Rdat"))
