### 1. Data preparation 

library(tidyverse)
library(lme4)
library(openxlsx)

### Load data ------
data_serpentine <- readRDS("data_serpentine.RDS")
soil <- data_serpentine[[1]]
traits.itv <- data_serpentine[[2]]
abundance <- data_serpentine[[3]]

# Number of species with traits in each plot
sp_rich_tr <- apply(abundance > 0, 1, sum) 

# average biomass of species with traits 
mean(soil$PercBiomObs)
min(soil$PercBiomObs)
max(soil$PercBiomObs)

# ## relationship between species richness, biomass and soil Ni
summary(lm(soil$richness ~ log(soil$Ni)))
summary(lm(soil$Greenbiom ~ log(soil$Ni)))

## Descriptive table of plot characteristics
soil |> 
  mutate(plot_no = 1:26) |>
  mutate(sp_traits = sp_rich_tr) |>
  select(plot_no, Site, richness, sp_traits, biomass = Greenbiom, PercBiomObs, Ni) |>
  write.csv("output/site_summary.csv")


## proportion of species without N measured
a <- as.vector(rowSums(abundance))
prop_nN <- abundance |> 
  mutate(ab_tot = as.vector(rowSums(abundance))) |>
  rowwise() |> 
  mutate(tot_nN = sum(BROMCOM,PHLEPRAT, TRIFCAMP)/ab_tot)
mean(prop_nN$tot_nN)
min(prop_nN$tot_nN)
max(prop_nN$tot_nN)


## Data preparation ####
# Clean and transform soil
soil.b <- soil %>% 
  mutate(MgCaRatio = log(Mg/Ca),
         Ni = log(Ni), 
         Cr = log(Cr), 
         Co = log(Co), 
         Mn = log(Mn), 
         P = log(P),  
         K = log(K)) %>%
  select(site = Site, Plot, Ni, Cr, Co, MgCaRatio, Mn, Zn, P, K)

## Transform abundance 
abun <- sqrt(abundance) |> rownames_to_column("site")

### Traits 
## Add leaf shape and transform trait values if necessary to improve distribution
traits.itv <- traits.itv %>% 
  filter(Sp %in% colnames(abun)) |>
  mutate(LS = LL/as.numeric(LW)) %>%
  mutate(LeafArea = LeafArea*100) %>% #transform in mm2
  mutate(LDMC = LDMC,
         LeafArea = sqrt(LeafArea),
         SLA = SLA,
         LT = sqrt(LT),
         LNC = sqrt(N), 
         LS = sqrt(LS)) |> 
  # remove impossible values
  filter(is.na(LDMC)| LDMC < 1000)

# Species level trait 
traits_sp = traits.itv %>%
  select(-species) %>%
  summarise(.by = "Sp",
            LS = median(LS, na.rm=TRUE),
            LeafArea = median(LeafArea, na.rm=TRUE),
            LDMC = median(LDMC, na.rm=TRUE),
            SLA = median(SLA, na.rm=TRUE),
            LT = median(LT, na.rm=TRUE),
            LNC = median(N,  na.rm=TRUE)) %>%
  ungroup() 


## average trait of each species in each site
traits_sp_site <- traits.itv |>
  select(Sp, site, LeafArea, SLA, LDMC, LT, LNC = N, LS) |>
  unite(sp_site, c(Sp, site), sep = "_", remove = F) |>
  group_by(sp_site, Sp, site) |>
  summarise_all(median, na.rm = T)


## Calculate missing values for species:site combinations missing
traits_sp_site <- abun |> 
  pivot_longer(where(is.numeric), names_to = "Sp", values_to = "abundance") |>
  # keep sp_site with abundance > 0
  filter(abundance > 0) |>
  # add the sp_site identifier
  separate(site, into = c("Site"), remove = T) |>
  unite(sp_site, c("Sp", "Site"), sep = "_", remove = T) |>
  group_by(sp_site) |>
  summarise(abundance = sum(abundance)) |>
  left_join(traits_sp_site) |> 
  separate(sp_site, into = c("Sp", "site")) |>
  left_join(soil |>
              select(site = Site, Ni)|>
              group_by(site) |>
              summarise(Ni = log(mean(Ni))))

# run models of traits along the gradient with intercept and slope per species
mod_LA_md <- lmer(LeafArea ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_SLA_md <- lmer(SLA ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LDMC_md <- lmer(LDMC ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LT_md <- lmer(LT ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LNC_md <- lmer(LNC ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LS_md <- lmer(LS ~ 0 + (1 + Ni| Sp), data = traits_sp_site)


mod_pred = data.frame(Sp = traits_sp_site$Sp,
                      site = traits_sp_site$site,
                      pred_LA = predict(mod_LA_md, newdata = traits_sp_site),
                      pred_SLA = predict(mod_SLA_md, newdata = traits_sp_site),
                      pred_LDMC = predict(mod_LDMC_md, newdata = traits_sp_site),
                      pred_LT = predict(mod_LT_md, newdata = traits_sp_site),
                      pred_LNC = predict(mod_LNC_md, newdata = traits_sp_site, allow.new.levels = TRUE),
                      pred_LS = predict(mod_LS_md, newdata = traits_sp_site)
)

# replace 0 by NA because is due to new levels in the prediction
mod_pred[mod_pred == 0] <- NA 

## replace in the trait_sp_site the non measured values by
## the predicted values
traits_sp_site[c(13, 16, 17, 20, 25, 52), 4:9] <- mod_pred[c(13, 16, 17, 20, 25, 52), 3:8]
traits_sp_site[9, 8] <- mod_pred[9, 7]

traits_sp_site <- select(traits_sp_site, - abundance)

## Write resulting data
# dir.create("./data")

write.csv(abun, "./data/abundance.csv")
write.csv(soil.b, "./data/soil.csv")
write.csv(traits_sp, "./data/traits_sp.csv")
write.csv(traits_sp_site,"./data/traits_sp_site.csv")
write.csv(traits.itv, "./data/traits_itv.csv")

