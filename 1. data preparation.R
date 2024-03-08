### Data preparation 
library(tidyverse)
library(lme4)
library(openxlsx)

### Data --------------------------------------------------------
soil <- read.xlsx("./data/soil.xlsx")
traits.itv <- read.xlsx("./data/traitsITV.xlsx")
abundance <- read.xlsx("./data/abundance.xlsx", rowNames = T) 
# setwd("../")

#species richness in each plot
apply(abundance>0, 1, sum) 

# ## relationship between species richness, biomass and soil Ni
summary(lm(soil$richness ~ log(soil$Ni)))
summary(lm(soil$Allbiom ~ log(soil$Ni)))
summary(lm(soil$Greenbiom ~ log(soil$Ni)))


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
  group_by(Site) %>% 
  mutate(id = row_number()) |>
  ungroup() |>
  unite(id, c(Site, id), sep= "_", remove = F) |>
  select(site = Site, id, Ni, Cr, Co, MgCaRatio, Mn, Zn, P, K)

# Transform abundance 
abundance_b <- abundance |>
  rownames_to_column("site") |>
  separate(site, into = c("site")) |>
  group_by(site) %>% 
  mutate(id = row_number()) |>
  ungroup() |>
  unite(id, c(site, id), sep= "_", remove = F) |>
  select(-site) |>
  column_to_rownames("id")

## Transform abundance (either square root or Hellinger)
abun <- sqrt(abundance_b) |> rownames_to_column("site")
# abun <- hellinger(abundance_b)

## Add leaf shape 
traits.itv <- traits.itv %>%
  ## transform trait values if necessary to improve distribution
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

# trait_sp
traits_sp = traits.itv %>%
  select(-species) %>%
  summarise(.by = "Sp",
            LS = median(LS, na.rm=TRUE),
            LeafArea = median(LeafArea, na.rm=TRUE),
            LDMC = median(LDMC, na.rm=TRUE),
            SLA = median(SLA, na.rm=TRUE),
            LT = median(LT, na.rm=TRUE),
            LNC = median(N,  na.rm=TRUE)) %>%
  ungroup %>%
  #column_to_rownames("Sp") %>%
  filter(!row_number() %in% c(16)) # Remove a species present in none of the sites 

## trait sp_site
traits_sp_site <- traits.itv |>
  select(Sp, site, LeafArea, SLA, LDMC, LT, LNC = N, LS) |>
  unite(sp_site, c(Sp, site), sep = "_", remove = F) |>
  group_by(sp_site, Sp, site) |>
  summarise_all(median, na.rm = T)

## join with abundance > 0 (using square root of abundance)
traits_sp_site <- abun |> 
  #rownames_to_column("site") |>
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

mod_LA  <- lmer(LeafArea ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_SLA <- lmer(SLA ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LDMC <- lmer(LDMC ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LT <- lmer(LT ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LNC <- lmer(LNC ~ 0 + (1 + Ni| Sp), data = traits_sp_site)
mod_LS <- lmer(LS ~ 0 + (1 + Ni| Sp), data = traits_sp_site)

mod_pred = data.frame(Sp = traits_sp_site$Sp,
                      site = traits_sp_site$site,
                      pred_LA = predict(mod_LA, newdata = traits_sp_site),
                      pred_SLA = predict(mod_SLA, newdata = traits_sp_site),
                      pred_LDMC = predict(mod_LDMC, newdata = traits_sp_site),
                      pred_LT = predict(mod_LT, newdata = traits_sp_site),
                      pred_LNC = predict(mod_LNC, newdata = traits_sp_site, allow.new.levels = TRUE),
                      pred_LS = predict(mod_LS, newdata = traits_sp_site)
)
# replace 0 by NA because is due to new levels in the prediction
mod_pred[mod_pred == 0] <- NA 

# Test model fit 
# ggplot(traits_sp_site, aes(x = Ni, y = SLA, colour = Sp)) +
#   geom_point() +
#   geom_line(aes(y = mod_pred$pred_SLA)) +
#   facet_wrap(~Sp)

## replace in the trait_sp_site the non measured values by
## the predicted values
traits_sp_site[c(13, 16, 17, 25, 52, 58, 59, 61), 4:9] <- mod_pred[c(13, 16, 17, 25, 52, 58, 59, 61), 3:8]
traits_sp_site <- select(traits_sp_site, - abundance)

# # Transform similarly to the median trait values => no, done before !
# traits_sp_site <- traits_sp_site |>
#   mutate(LeafArea = log(LeafArea), 
#          LT = log(LT), 
#          LNC = log(LNC),
#          LS = log(LS))

write.csv(abun, "./data/abundance_clean.csv")
write.csv(soil.b, "./data/soil_clean.csv")
write.csv(traits_sp, "./data/traits_sp.csv")
write.csv(traits_sp_site,"./data/traits_sp_site.csv")
write.csv(traits.itv, "./data/traits_itv_clean.csv")
