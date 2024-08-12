# 4. Intraspecific trait variation 
library(tidyverse)
library(brms)
library(tidybayes)
library(lmerTest)
library(FD)
library(rstan)
library(labdsv)
library(marginaleffects)
library(modelr)
source("Functions.R")

# load data
soil.b <- read.csv("./data/soil.csv") |> select(-X)
traits_sp_site <- read.csv("./data/traits_sp_site.csv") |> select(-X)
traits_sp <- read.csv("./data/traits_sp.csv") |> select(-X)
abun <- read.csv("./data/abundance.csv") |> select(-X)
traits.itv <- read.csv("./data/traits_itv.csv")|> select(-X)
CWM_sp_site <- read.csv("./output/CWM/CWM_sp_site.csv")|> select(-X)
CWM_sp  <- read.csv("./output/CWM/CWM_sp.csv")|> select(-X)

# dir.create("./output/varpart")

## 1 Variation partitionning # ------
# Data preparation
traits.norm <- traits.itv  |>
  select(site, Sp, LeafArea, LDMC, SLA, LT, LNC, LS )|>
  left_join(soil.b |> 
              select(site, Ni) |> 
              group_by(site) |> 
              summarize(Ni = mean(Ni))) |>
  mutate_at(c("LS","LeafArea", "LDMC", "SLA", "LT", "LNC", "Ni"), scale2)# standardise traits

# Test variance partition between intra and inter
fit_N <- brm(formula = SLA ~ 1 + Ni + (1 + Ni | Sp),
             data = traits.norm,
             control = list(adapt_delta = 0.95, 
                            max_treedepth = 13),
             iter = 3000, warmup = 1000, 
             file = "./output/varpart/N"
)

fit_SLA <- update(fit_N, SLA ~ ., newdata = traits.norm, file = "./output/varpart/SLA")
fit_LA <- update(fit_N, LeafArea ~ ., newdata = traits.norm, file = "./output/varpart/LA")
fit_LT <- update(fit_N, LT ~ ., newdata = traits.norm, file = "./output/varpart/LT")
fit_LDMC <- update(fit_N, LDMC ~ ., newdata = traits.norm, file = "./output/varpart/LDMC")
fit_LS <- update(fit_N, LS ~ ., newdata = traits.norm, file = "./output/varpart/LS")

fit_N <- readRDS("./output/varpart/N.rds")
fit_SLA <- readRDS("./output/varpart/SLA.rds")
fit_LA <- readRDS("./output/varpart/LA.rds")
fit_LT <- readRDS("./output/varpart/LT.rds")
fit_LDMC <- readRDS("./output/varpart/LDMC.rds")
fit_LS <- readRDS("./output/varpart/LS.rds")

varpart_result <-rbind(
  varpart_fit(fit_SLA, "SLA"),
  varpart_fit(fit_N, "N"), 
  varpart_fit(fit_LA, "LA"), 
  varpart_fit(fit_LT, "LT"), 
  varpart_fit(fit_LDMC, "LDMC"), 
  varpart_fit(fit_LS, "LS")
)

# write.csv(varpart_result, "./output/varpart/varpart.csv")


#### Figure varpart
var <- varpart_result %>%
  filter(str_detect(component, "pct")) %>%
  select(Component = component, Trait=trait, Var = y, Varmin = ymin, Varmax = ymax) %>%
  mutate(intrainter = fct_relevel(Component, 
                                  "Sp_pct", "Sp__Ni_sq_pct", "residual_pct"), 
         Trait = fct_relevel(Trait, "LA", "SLA", "LT", "LDMC", "N", "LS"))

var |> select(Component, Trait, Var) |>
  pivot_wider(names_from = Trait, values_from = Var)

ggplot(var, aes(fill=intrainter, y=Var, x=Trait)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = Varmin, ymax = Varmax), width=.2,
                position=position_dodge(.9)) +
  labs(x = "", y = "Variance", fill = "") +
  theme_minimal()+
  scale_fill_grey(labels = c("Interspecific", "Intraspecific", "Residual"))

# ggsave("./output/figure/varpart.pdf", height =3, width =  5)

## 3.2 species level variation variation along the gradient
var_Leaf_area = mod_intra(fit_LA)
var_SLA = mod_intra(fit_SLA)
var_LT = mod_intra(fit_LT)
var_LDMC = mod_intra(fit_LDMC)
var_LNC = mod_intra(fit_N)
var_LS = mod_intra(fit_LS)


var_intra <- data.frame(
  Species = colnames(abun[2:21]), var_Leaf_area, var_SLA, var_LT,
  var_LDMC, var_LNC, var_LS)

# write.csv(var_intra, "./output/varpart/var_intra_table.csv")

ggplot(traits.norm,aes(y = LT, x = Ni, colour = Sp)) +
  geom_point() +
  geom_smooth(method = "lm") + 
  facet_wrap(~Sp)

# 3.3 Use modified trait.flex function from Leps 2011 to calculate 
# contribution of ITV to CWM
CWM_site <- CWM_sp_site |> 
  rename(CWM_LA_intra = CWM_LA, CWM_SLA_intra = CWM_SLA, CWM_LT_intra = CWM_LT,
         CWM_LDMC_intra = CWM_LDMC, CWM_LNC_intra = CWM_LNC, CWM_LS_intra = CWM_LS)
CWM_med <- CWM_sp |> 
  rename(CWM_LA_inter = CWM_LA, CWM_SLA_inter = CWM_SLA, CWM_LT_inter = CWM_LT,
         CWM_LDMC_inter = CWM_LDMC, CWM_LNC_inter = CWM_LNC, CWM_LS_inter = CWM_LS)  

CWM <- cbind(CWM_site, CWM_med)

t_fl_LA <- trait.flex(specific = "CWM_LA_intra", nonspecific = "CWM_LA_inter",
                      envi = "Ni", data = CWM)
t_fl_SLA <- trait.flex(specific = "CWM_SLA_intra", nonspecific = "CWM_SLA_inter",
                       envi = "Ni", data = CWM)
t_fl_LT <- trait.flex(specific = "CWM_LT_intra", nonspecific = "CWM_LT_inter",
                      envi = "Ni", data = CWM)
t_fl_LDMC <- trait.flex(specific = "CWM_LDMC_intra", nonspecific = "CWM_LDMC_inter",
                        envi = "Ni", data = CWM)
t_fl_LNC <- trait.flex(specific = "CWM_LNC_intra", nonspecific = "CWM_LNC_inter",
                       envi = "Ni", data = CWM)
t_fl_LS <- trait.flex(specific = "CWM_LS_intra", nonspecific = "CWM_LS_inter",
                      envi = "Ni", data = CWM)
t_fl <- rbind(LA = t_fl_LA$RelSumSq[1,], 
              SLA = t_fl_SLA$RelSumSq[1,], 
              LDMC = t_fl_LDMC$RelSumSq[1,], 
              LT = t_fl_LT$RelSumSq[1,], 
              LNC = t_fl_LNC$RelSumSq[1,], 
              LS = t_fl_LS$RelSumSq[1,]) |>
  rownames_to_column("trait") |>
  pivot_longer(where(is.numeric), names_to = "variable") 
t_fl$variable <-  factor(t_fl$variable, level = c("Turnover", "Intraspec.", "Covariation", "Total"))
t_fl$trait <-  factor(t_fl$trait, level = c("LA", "SLA", "LDMC", "LT", "LNC","LS"))

ggplot(t_fl, aes(fill=variable, y=value, x=trait)) + 
  geom_bar(position="dodge", stat="identity") +
  labs(x = "", y = "Variance", fill = "") +
  theme_minimal()+
  scale_fill_grey()

# ggsave("./output/figure/varpart_cwm.jpeg", height = 4, width = 7)

## 3.3 intraspecific trait covariation and diversity ----------
traits_integration <- traits.itv |> 
  select(Sp, site, LeafArea, LDMC, SLA, LT) |>
  mutate_at(.vars = c("LeafArea", "LDMC", "SLA", "LT"), scale2)|>
  unite(sp_site, c(Sp, site), sep = "_") |>
  distinct()|>
  na.omit()

trait_int <- data.frame(population = unique(traits_integration$sp_site), 
                        trait_integration = NA, 
                        FRic = NA, 
                        FDis = NA)

for(i in 1:length(trait_int$population)) {
  pop <- unique(traits_integration$sp_site)[i]
  tr <- traits_integration |> filter(sp_site == pop)
  FD <- FD::dbFD(x = tr[,2:5])
  
  trait_int[i,2] <- var(eigen(cor(tr[,2:5], use = "complete.obs"))$values)
  trait_int[i,3] <- FD$FRic
  trait_int[i,4] <- FD$FDis
}

soil_pop <- soil.b |> 
  select(site, Ni) |>
  group_by(site) |>
  summarise(Ni = mean(Ni)) |>
  mutate(Ni = scale(Ni))

data_int <- trait_int |> 
  separate(population, into = c("species", "site")) |>
  left_join(soil_pop) |>
  filter(species != "TRIFCAMP" & species != "BROMCOM")

mod_coordination <- brm(trait_integration ~ 1 + Ni + (1 + Ni | species), 
                        data = data_int, 
                        prior = c(prior(normal(0, 2), class = b, coef = Ni),
                                  prior(normal(1, 2), class = Intercept)))

mod_FRic <- brm(FRic ~ 1 + Ni + (1 + Ni | species), 
                          data = data_int,
                prior = c(prior(normal(0, 2), class = b, coef = Ni),
                          prior(normal(1, 2), class = Intercept)))

mod_FDis <- brm(FDis ~ 1 + Ni + (1 + Ni | species), 
                          data = data_int, 
                prior = c(prior(normal(0, 2), class = b, coef = Ni),
                          prior(normal(1, 2), class = Intercept)))

