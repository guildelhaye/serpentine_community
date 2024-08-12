## Helper functions

# function to scale the data
scale2 <- function(x, na.rm = T) (x - mean(x, na.rm = na.rm)) / sd(x, na.rm)

## Helper function to extract the proportion of variance
varpart_fit <- function(model, trait){
  model %>% 
    tidy_draws() %>% 
    select(starts_with("sd_"), sigma) %>% 
    transmute_all(.funs = list(sq = ~(. ^ 2))) %>% 
    mutate(total_var = rowSums(.)) %>% 
    mutate_at(.vars = vars(-total_var), 
              .funs = list(pct = ~(. / total_var))) %>% 
    map_df(.f = ~ median_hdci(., .width = 0.95), .id = "component") %>% 
    mutate(
      component = str_remove(component, "sd_"),
      component = str_remove(component, "__Intercept_sq"),
      component = str_replace(component, "sigma_sq", "residual")) %>%
    mutate(trait = paste0(trait))
}

## Function to extract species specific slopes from a 
## multilevel model with species as random variables
mod_intra <- function(mod){
  fit <- mod %>% 
    spread_draws(b_Ni, r_Sp[species,Ni]) %>%
    filter(Ni == "Ni") |>
    mutate(slope = r_Sp) |>
    ungroup() |>
    select(species, slope) |>
    group_by(species) |>
    median_qi(slope, .width = c(.95)) %>%
    mutate(lower = round(.lower, 2),
           upper = round(.upper, 2),
           slope = round(slope, 2)) |>
    unite("CI", c(lower, upper), sep = " : ") |>
    mutate(slopeCI = paste0(slope, " (", CI, ")"))
  
  pl <- fit|> 
    ggplot(aes(y = species, x = slope, xmin = .lower, xmax = .upper)) +
    geom_pointinterval()  +
    geom_vline(xintercept = 0, linetype = "dashed")
  print(pl)
  return(fit$slopeCI)
}

## function to extract draws
m_ex <- function(model, trait) {
  ax <- model %>% tidy_draws() %>% 
    select(var = starts_with("b_Ni:")) |> 
    mutate(trait = trait)
}

m_ex2 <- function(model, trait) {
  ax <- model %>% tidy_draws() %>% 
    select(var = starts_with("b_Ni")) |> 
    mutate(trait = trait)
}


## Function extract quantiles
extract_quantile <- function(model) {
  post <- model %>% 
    tidy_draws() %>% 
    select(var = starts_with("b_Ni:"))
  quant <- quantile(post$var, probs = c(0.05, 0.1, 0.15, 0.5, 0.85,0.9, 0.95))
  p_above_zero <- sum(post$var > 0)/length(post$var)
  res <- list(quantiles = quant, 'prob > 0' = p_above_zero)
  return(res)}

## Function trait flex modified from Leps et al 2011 Ecography
trait.flex = function(specific, nonspecific, envi, data){
  
  data = CWM
  envi = as.character(envi)
  specific = as.character(specific)
  nonspecific = as.character(nonspecific)
  
  trait <- data.frame(
    envi = data |> select(all_of(envi)),
    specific = as.vector(data |> select(all_of(specific))),
    nonspec = as.vector(data |> select(all_of(nonspecific))))
  
  colnames(trait) <- c("envi", "specific", "nonspec")
  
  trait <- trait |> mutate(diff = specific - nonspec)
  
  res.1 <- summary(aov(lm(specific ~ envi, data = trait)))[[1]]
  res.2 <- summary(aov(lm(nonspec ~ envi, data = trait)))[[1]]
  res.3 <- summary(aov(lm(diff ~ envi, data = trait)))[[1]]
  
  nrows <- dim(res.1)[1]
  ss.turn <- res.2[,2]
  ss.var  <- res.3[,2]
  ss.tot  <- res.1[,2]
  ss.covar <- ss.tot - ss.turn - ss.var
  ss.row.names <- dimnames(res.1)[[1]]
  
  if(nrows > 1) {
    ss.turn <- c(ss.turn, sum(ss.turn))
    ss.var  <- c(ss.var,  sum(ss.var))
    ss.tot  <- c(ss.tot,  sum(ss.tot))
    ss.covar<- c(ss.covar,sum(ss.covar))
    ss.row.names <- c(ss.row.names, "Total")
    nrows   <- nrows + 1
  } else {
    # replace row title
    ss.row.names[1] <- "Total"
  }
  
  SS.tab <- data.frame( Turnover = ss.turn, 
                        Intraspec. = ss.var,
                        Covariation = ss.covar, 
                        Total = ss.tot,
                        row.names = ss.row.names)
  
  # Calculate relative fractions
  TotalSS <- SS.tab[nrows, 4] # lower right corner
  SS.tab.rel <- SS.tab / TotalSS
  
  # Collect significances
  if(nrows > 1)  # get rid of the "Total" label again
    ss.row.names <- ss.row.names[-nrows]
  P.tab <- data.frame( Turnover = res.2[,5], Intraspec. = res.3[,5],
                       Total = res.1[,5], row.names = ss.row.names)
  
  res <- list( #SumSq=SS.tab, 
    RelSumSq=SS.tab.rel , Pvals=P.tab #, 
    #anova.turnover=res.2, anova.total=res.1, anova.diff=res.3
  )
  return(res)
}

#### CWM using long format tables of traits and abundance to include easily intraspecific trait variation
cwm <- function(abundance, traits) {
  com_mean  <-  abundance |>
    left_join(traits) |> 
    mutate(s_LA = abundance * LeafArea, 
           s_SLA = abundance * SLA,
           s_LT = abundance * LT,
           s_LDMC = abundance * LDMC,
           s_LNC = abundance * LNC,
           s_LS = abundance * LS) |>
    group_by(id) |>
    summarise(CWM_LA = sum(s_LA, na.rm = T)/sum(abundance, na.rm = T), 
              CWM_SLA = sum(s_SLA, na.rm = T)/sum(abundance, na.rm = T), 
              CWM_LT = sum(s_LT, na.rm = T)/sum(abundance, na.rm = T), 
              CWM_LDMC = sum(s_LDMC, na.rm = T)/sum(abundance, na.rm = T), 
              CWM_LNC = sum(s_LNC, na.rm = T)/sum(abundance, na.rm = T), 
              CWM_LS = sum(s_LS, na.rm = T)/sum(abundance, na.rm = T))|> 
    select(id, starts_with("CWM_"))
  return(com_mean)
}

#### Univariate FDis
sesFDisu <- function(traits, comm, i, si, numberReps){ # i is the trait name in ' ', ex: 'SLA'
# traits = traits_sp
#   comm = abun
#   i = "SLA"
#   si = "all"
#   numberReps = 10

  if(si != "all") {
    comm = comm |> filter(str_detect(site, si)) |> select(-site)
    traits <- traits |> filter(site == si) |> column_to_rownames("Sp") |> select(-site)
  } else {
    traits <- traits |> column_to_rownames("Sp")
    }
  
  traits <- traits[rownames(traits) %in% colnames(comm),]
  Tr <- traits %>% select(all_of(i)) %>% filter(is.na(.) == F)
  comm <- comm[,colnames(comm) %in% rownames(Tr)]
  
  FDobs <- dbFD(Tr, comm, calc.FGR = F, 
                calc.CWM = F, calc.FDiv = F, calc.FRic = F)
  #Lets create a matrix to store results from each iteration (one column per iteration) 
  resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                          dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  for(rep in 1:numberReps){ 
    traitsRand <- data.frame(Tr[sample(1:nrow(Tr)),])
    rownames(traitsRand) <- rownames(Tr) 
    FDnull <- dbFD(traitsRand, comm, m = ndim, calc.FGR = F, 
                   calc.CWM = F, calc.FDiv = F, calc.FRic = F, message = F)
    resultsRandom[,rep] <- FDnull$FDis
  }
  obsFDis <- FDobs$FDis
  meanNullFDis <- rowMeans(resultsRandom) 
  ES <- obsFDis - meanNullFDis
  sdNull <- apply(resultsRandom, 1, sd) 
  SESFDis <- ES / sdNull
  return(SESFDis)
}

### Multivariate FDis ####
sesFDism <- function(traits, comm, si, ndim, numberReps){#ndim is the number of axes to keep in the PCoA

  if(si != "all") {
    comm = comm |> filter(str_detect(site, si)) |> select(-site)
    traits <- traits |> filter(site == si) |> column_to_rownames("Sp") |> select(-site)
  } else {
    traits <- traits |> column_to_rownames("Sp")
  }
  
  comm <- comm[,colnames(comm) %in% rownames(traits)]
  traits <- traits[rownames(traits) %in% colnames(comm),]
  
  # rownames(traits) %in% colnames(comm)
  # colnames(comm) %in% rownames(traits)
   
  FDobs <- dbFD(dist(traits), comm, m = ndim, calc.FGR = F, 
                calc.CWM = F, calc.FDiv = F, 
                calc.FRic = F)#, corr = c("lingoes"))
  #Lets create a matrix to store results from each iteration (one column per iteration) 
  resultsRandom <- matrix(NA, nrow = nrow(comm), ncol = numberReps, 
                          dimnames = list(rownames(comm), paste0("Sim.", 1:numberReps)))
  for(rep in 1:numberReps){ 
    traitsRand <- traits[sample(1:nrow(traits)),] 
    rownames(traitsRand) <- rownames(traits) 
    FDnull <- dbFD(dist(traitsRand), comm, m = ndim, calc.FGR = F, 
                   calc.CWM = F, calc.FDiv = F, calc.FRic = F, message = F, corr = c("lingoes"))
    resultsRandom[,rep] <- FDnull$FDis
  }
  
  obsFDis <- FDobs$FDis
  meanNullFDis <- rowMeans(resultsRandom) 
  ES <- obsFDis - meanNullFDis
  sdNull <- apply(resultsRandom, 1, sd) 
  SESFDis <- ES / sdNull
  return(SESFDis)
}

#### Plot FDis figure
plot_fdis <- function(mod_1, mod_2, trait) {
  fdis <- data_fdis %>%
    data_grid(Ni = seq_range(Ni, n = 10), site) %>%
    add_epred_draws(mod_1) |>
    mutate(model = rep("total"))
  
  fdis_OL <- data_fdis_OL %>%
    data_grid(Ni = seq_range(Ni, n = 10), site) %>%
    add_epred_draws(mod_2) |>
    mutate(model = rep("no O.lesb"))
  
  rbind(fdis, fdis_OL) |>
    ggplot(aes(x = Ni, y = .epred, colour = model)) +
    stat_lineribbon(alpha = 0.5) +
    scale_fill_brewer(palette = "Greys") +
    scale_color_brewer(palette = "Set2") + 
    theme_classic() +
    theme(legend.position = "none") +
    ylim(-4, 4)+
    ylab(paste0(trait))+
    xlab("Soil Ni content (z-score)")
}
