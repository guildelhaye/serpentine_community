### 2. Exploratory analyses 
library(tidyverse)
library(vegan)
library(ggfortify)

## load data



soil <- read.csv("./data/soil.csv") |> select(-c(X))
traits_sp_site <- read.csv("./data/traits_sp_site.csv") |> select(-X)
abun <- read.csv("./data/abundance.csv") |> select(-X)
traits_sp <- read.csv("./data/traits_sp.csv") |> select(-X)

# dir.create("./output")
# dir.create("./output/figure")

###### 1. Correlation structures in soil and trait data -----
## Soil analyses ----------------------------------------------
# correlation between variables
plot(soil[,3:10])
round(cor(soil[,3:10]), 2)

## PCA
PC2 <- prcomp(soil[,3:10], center = TRUE, scale = TRUE)

#Plot
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
autoplot(PC2, data = soil, colour = "site", type = "site", loadings = TRUE,
         loadings.label = TRUE, loadings.label.size = 4, loadings.colour = 'black', 
         loadings.label.colour = 'black', frame = TRUE) +
  theme_minimal() +
  theme(legend.position = "none")

# ggsave("./output/figure/soil_pca.pdf")

##Relationship between traits ----------------------------
traits_pca <- traits_sp_site |> 
  na.omit()|>
  left_join(soil |> 
              select(site, Ni) |> 
              summarise(.by = site, Ni = log(mean(Ni)))) |>
  mutate(Sp = recode(Sp,
                     "AEGBIUNC"= "A. biuncialis",
                     "ALYLESB" = "O. lesbiaca", 
                     "ANAGARVE" = "A. arvensis",
                     "AVENA" = "A. barbata",
                     "CREPCOM" = "C. commutata",
                     "CYNOECHI" = "C. echinatus",
                     "DACTGLOM" = "D. glomerata",
                     "FILAERIO" = "F. eriocephala",
                     "HIRSINCA" = "H. incana",
                     "HORDBULB" = "H. bulbosum",
                     "LAGOCUMI" = "L. cuminoides",
                     "LOLIRIGI" = "L. rigidum",
                     "PLANLAGO" = "P. lagopus",
                     "TORINODO" = "T. nodosa",
                     "TRACHDIST" = "T. distachya",
                     "TRIFANGU" = "T. angustifolium",
                     "TRIFARVE" = "T. arvense"))

# PCA
PC2 <- prcomp(na.omit(traits_pca[,c(3:8)]), center = TRUE, scale = TRUE) #check

#Plot
par(mar = c(4.1, 4.4, 4.1, 1.9), xaxs="i", yaxs="i")
autoplot(PC2, data = traits_pca, colour = "Sp", 
         loadings = TRUE,
         loadings.label = TRUE, 
         loadings.label.size = 4, 
         loadings.colour = 'black', 
         loadings.label.colour = 'black', 
         frame = TRUE) +
  theme_minimal() 


# ggsave("./output/figure/trait_pca.pdf")

## Correlation between traits at sp level
round(cor(traits_sp[,c(2:6)]), 2)

## NMDS 
nmds = metaMDS(abun[,2:21], distance = "bray")
en = envfit(nmds, soil[,c(3:10)], permutations = 999, na.rm = TRUE)


ordiplot(nmds,type="points", xlim = c(-2.5, 2))
orditorp(nmds,display="species",col="red",air=0.01, cex=0.5)
plot(en)

# ggsave("./output/figure/NMDS.pdf", height = 4, width = 5)

#without A. lesb
nmds.nal = metaMDS(abun[,c(2, 4:21)], distance = "bray")
ordiplot(nmds.nal,type="points", xlim = c(-2.5, 2))
orditorp(nmds.nal,display="species",col="red",air=0.01, cex=0.5)
plot(en)

## Proportion of annual species
# 3 species are perennial, DACTGLOM, HORDBULB, PHLEPRAT
abund_annual <- abun |> 
  mutate(abun_tot = rowSums(abun[,-3] |> # remove O. lesbiaca
                              column_to_rownames("site"))) |>
           rowwise() |> 
  mutate(prop_annual = 1 - (sum(DACTGLOM, HORDBULB, PHLEPRAT)/abun_tot)) |> # annual = 1 - perennial
  left_join(soil |> select(site = Plot, Ni)) 

summary(lm(prop_annual ~ Ni, data = abund_annual))

ggplot(abund_annual, aes(y = prop_annual, x = Ni)) +
  geom_point() +
  geom_smooth(method = "lm") +
  ylab("Proportion annual species") +
  xlab("Soil Ni content (log)")

# ggsave("./output/figure/proportion_annual.pdf", height = 4, width = 4)
