# JMH Feb 2022
# Fits Langmuir model to isotherm data
# merges with other isotherm data

# Libaries ----
library(tidyverse)
library(nlstools) #for 95% conf int 
# # this package is better at choosing initials
library(minpack.lm)
library(broom)

# get data ----
iso <- read.csv("01_RawData/SedExp_R_WK.csv") %>% 
                # get isotherm exp data
                filter(Experiment == "Psorp") %>% 
                select(Date, Stream, Added_ugP_L:MassSed_g) %>% 
                mutate(Date = ifelse(Date == "Jan-24-2019", "Jan-24-19", Date)) %>% 
                mutate(Date = as.POSIXct(Date, format = "%b-%d-%y"),
                       # 0.04 L is the volume of water)
                       P_change_ugP_g = ((Added_ugP_L - Final_ugP_L)*0.04)/MassSed_g,
                       Stream = ifelse(Stream == "Platter ", "Platter",Stream),
                       # partitioning coefficient PU et al. 2021 Chemosphere, 263 128334 - L/mg
                       # THERE’S A MISTAKE in Lu final_ugP needs to be divided not multiplie
                       Kd_Pu = ((Added_ugP_L/1000) - (Final_ugP_L/1000))*0.04/(MassSed_g/1000) /(Final_ugP_L/1000),
                       # Kd_Pu2 = P_change_ugP_g/Final_ugP_L
                       ) %>% 
                mutate(ID = paste0(Stream,"_", Date))


AmbDat <- read.csv("01_RawData/RawSorptionData.csv")[,-c(1,2)] %>% 
              mutate(Date = as.POSIXct(SampleDate, format = "%m/%d/%y")) %>% 
              mutate(Stream = case_when(Stream == "Little Flat Rock" ~ "LFR",
                                        Stream == "Platter Creek" ~ "Platter",
                                        Stream == "Unnamed Trib to Lost Creek" ~ "Lost_Trib",
                                        Stream == "South Turkey Foot" ~ "S_Turkey",
                                        Stream == "West Creek" ~ "West",
                                        Stream == "Potato Run" ~ "Potato_Run")) %>% 
              mutate(ID = paste0(Stream,"_", Date))


iso2 <- iso %>% 
        full_join(AmbDat %>% 
                    select(-Stream, -Date), by = "ID")

# looking for weird mass values
hist(iso2$MassSed_g)

# Look at the isotherm relationships
    ggplot(iso2, aes(y = P_change_ugP_g, x = Final_ugP_L)) +
      geom_point()+
      facet_grid(Stream ~ Date) +
      geom_hline(yintercept = 0) +
      stat_smooth(se = F, color = "blue", size = 0.75, formula = y ~ x + I(x^2), method = "lm")
    
  # Fig for response to reviewr ----
    ggplot(iso2  %>% 
             filter(Added_ugP_L != 0) %>% 
             filter(Date != "2019-01-24 EST") %>% 
             filter(Date != "2019-02-08 EST"), aes(y = ((Added_ugP_L-Final_ugP_L)), 
                                                   x = Added_ugP_L, color = ID)) +
      geom_point()+
      # geom_line()+
      # facet_wrap(vars(ID)) +
      stat_smooth(method = "lm", size = 0.75, se = F) +
      geom_hline(yintercept = 0) +
      xlim(0,650) +
      theme(legend.position = "none") +
      ylab("Added µgP/L - Final µgP/L") +
      xlab("Added P (µg P/L)") +
      geom_abline(intercept = 0, slope = 1)
 
    
# Fit isotherm ---
    # langmuir parameters
    par0.Langmuirv2 <-c(
      K = 0.01, #bonding energy L/g
      Qmax = 300, #P sorption maximum ug P/g
      NAP = 0 # native P ug P/g
    )
    
    
    # Run all of them
    isoModParam <- iso2 %>% 
      mutate(ID = paste0(Stream,"_", Date)) %>% 
      filter(!(Stream == "Platter" & Date == "2019-06-30" & P_change_ugP_g < 5)) %>%
      filter(!(Stream == "Platter" & Date == "2019-01-24" & P_change_ugP_g > 200)) %>%
      filter(!(Stream == "West" & Date == "2019-01-24" & P_change_ugP_g > 200)) %>%
      filter(!(Stream == "LFR" & Date == "2019-02-08" & P_change_ugP_g >150)) %>%
      filter(!(Stream == "Platter" & Date == "2019-02-08" & P_change_ugP_g > 130)) %>%
      filter(!(Stream == "S_Turkey" & Date == "2019-02-08" & P_change_ugP_g >160)) %>%
      filter(!(Stream == "West" & Date == "2019-02-08" & P_change_ugP_g > 100)) %>%
      filter(!(Stream == "Potato_Run" & Date == "2019-03-10" & P_change_ugP_g > 100)) %>%
      # something's wrong with Whitney's parameters
      # filter(!(Stream == "LFR" & Date == "2019-04-12" & P_change_ugP_g > 90 & P_change_ugP_g < 0)) %>%
      filter(!(Stream == "S_Turkey" & Date == "2019-04-12" & P_change_ugP_g < 0)) %>%
      filter(!(Stream == "West" & Date == "2019-04-12" & P_change_ugP_g > 260)) %>%
      filter(!(Stream == "S_Turkey" & Date == "2019-04-20" & P_change_ugP_g > 100)) %>%
      filter(!(Stream == "Platter" & Date == "2019-05-01" & P_change_ugP_g > 90)) %>%
      filter(!(Stream == "S_Turkey" & Date == "2019-05-01" & P_change_ugP_g > 150)) %>%
      filter(!(Stream == "Potato_Run" & Date == "2019-05-20" & P_change_ugP_g > 100)) %>%
      filter(!(Stream == "Potato_Run" & Date == "2019-06-20" & P_change_ugP_g > 100)) %>%
      filter(!(Stream == "LFR" & Date == "2019-06-30" & P_change_ugP_g > 190)) %>%
      filter(!(Stream == "Platter" & Date == "2019-06-30" & P_change_ugP_g < 5)) %>%
      # sites that won't fit
      filter(!(ID %in% c("Lost_Trib_2019-04-12"))) %>%
      #LFR_2019-04-12
      nest(data = c(-ID)) %>% 
      mutate(
        fit = map(data, ~ nlsLM(P_change_ugP_g ~ ((Qmax * K * Final_ugP_L )/ (1 + K * Final_ugP_L)) - NAP, 
                                data = .x, start = par0.Langmuirv2, control = list(maxiter = 300), trace = TRUE)),
        tidied = map(fit, tidy)) %>% 
      unnest(tidied) %>% 
      select(-data, -fit)
    
    isoModParam2 <- isoModParam %>% 
      rename("est" = "estimate", "se" = "std.error", "t" = "statistic", "P" = "p.value") %>% 
      pivot_wider(id_cols = ID, names_from = term, values_from = est:P)
    
    
    # join with raw data ----
    iso3 <- iso2 %>% 
      mutate(ID = paste0(Stream,"_", Date)) %>% 
      full_join(isoModParam2, by = "ID") %>% 
      mutate(P_change_ugP_gMOD = ((est_Qmax * est_K * Final_ugP_L )/ (1 + est_K * Final_ugP_L)) - est_NAP,
             EPC0 = est_NAP/((est_Qmax * est_K) - (est_NAP * est_K)),
             EPCdif = EPC0-EPC,
              # but this needs to be done in the scaling
             P_change_ugP_gMODatAmb = ((est_Qmax * est_K * AmbP )/ (1 + est_K * AmbP)) - est_NAP,
             sorp_capacityJMH = est_Qmax - NAP) %>% 
      filter(!(SampleDate %in% c("1/24/19", "2/8/19"))) %>% 
      filter(ID != "LFR_2019-04-12")
    
    #Export isotherm dat----
    # write.csv(iso3, "01_RawData/00_WK_IsoAmb_JMH.csv")
 
    # Save ----
    # save.image(file.path(here::here("03_Rdata"),"00_EPCredo_Rdat"))
    
    