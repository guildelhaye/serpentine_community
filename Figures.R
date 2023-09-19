### Figures ###
library(tidyverse)
library(openxlsx)
library(ggpubr)

setwd("../output")

### Variance partitionning ---------------------------------------------------
var <- read.csv("varpart.csv") %>%
  filter(str_detect(component, "pct")) %>%
  select(Component = component, Trait=trait, Var = y, Varmin = ymin, Varmax = ymax) %>%
  mutate(intrainter = fct_relevel(Component, 
                                  "Sp_pct", "Sp:Serp_pct", "residual_pct"))

ggplot(var, aes(fill=intrainter, y=Var, x=Trait)) + 
  geom_bar(position="dodge", stat="identity") +
  geom_errorbar(aes(ymin = Varmin, ymax = Varmax), width=.2,
                position=position_dodge(.9)) +
  labs(x = "", y = "Variance", fill = "") +
  theme_minimal()+
  scale_fill_grey(labels = c("Interspecific", "Intraspecific", "Residual"))

### Functional diversity -----------------------------------------------------
# Four files. Results_FD_Final_Aver, Results_FD_Final_ITV, 

FD_results <- read.csv("Results_FD_Final_ITV_hel.csv") %>% 
   select(-X)  %>%
   mutate(CWM.LT = CWM.LT/1000, # transform lt in mm
          CWM.LDMC=CWM.LDMC/10) # transform mg/g in percent 

## Function to plot the graph with line type as function of significance test
# CWM ----------------
graph_cwm <- function(var,# "variable name within the dataset" 
                     name_var # "Name we want for the variable in the graph
                     ) {
 # var = "CWM.LeafArea"
  ## Calculate the regression coefficients
formula = paste0(var, " ~ log(Ni) ")
(a <- summary(lm(formula, data = FD_results %>% filter(Com == "Tot"))))
if(a$coefficients[2, 4] < 0.05){x = rep(TRUE, 26)} else {x = rep(FALSE, 26)}
(b <- summary(lm(formula, data = FD_results %>% filter(Com == "No_Alesb"))))
if(b$coefficients[2, 4] < 0.05){y = rep(TRUE, 26)} else {y = rep(FALSE, 26)}
# vector assigning statistical significance for both curves
(sig <- c(x, y))
#local dataset of trait + significance
data  = FD_results %>% select(paste0(var), Ni, Com) %>% cbind(sig)
 
plot <- ggplot(data, aes(x = Ni, y = data[,1], colour = Com))+
  geom_jitter(height = 0.05, alpha = 0.3)+
  geom_smooth(aes(linetype = sig), method = "lm", se = F) +
  scale_linetype_manual(values=c("TRUE"= "solid", "FALSE"="dashed")) +
  theme_classic() +
  xlab(label = "Soil Ni content (mg kg-1)")+
  ylab(label = paste0(name_var)) +
  theme(legend.position = "none") 
}

graph_LA_CWM <- graph_cwm("CWM.LeafArea", name_var = paste("Leaf area (mm2)"))
graph_SLA_CWM <- graph_cwm("CWM.SLA", name_var = paste("Specific leaf area (cm2 g-1)"))
graph_LDMC_CWM <- graph_cwm("CWM.LDMC", name_var = paste("Leaf dry matter content (%)"))
graph_LT_CWM <- graph_cwm("CWM.LT", name_var = paste("Leaf thickness (mm)"))
graph_LNC_CWM <- graph_cwm("CWM.LNC", name_var = paste("Leaf nitrogen content (%)"))
graph_LS_CWM <- graph_cwm("CWM.LS", name_var = paste("Leaf shape"))


ggarrange(graph_LA_CWM, graph_SLA_CWM, graph_LDMC_CWM, graph_LT_CWM, 
          graph_LNC_CWM , graph_LS_CWM, ncol = 3, nrow = 2)
ggsave("CWM_Aver_hel.pdf")

#FDis --------------------
graph_fdis<- function(var,# "variable name within the dataset" 
                     name_var # "Name we want for the variable in the graph
) {
  ## Calculate the regression coefficients
  formula = paste0(var, " ~ log(Ni)")
  (a <- summary(lm(formula, data = FD_results %>% filter(Com == "Tot"))))
  if(a$coefficients[2, 4] < 0.05){x = rep(TRUE, 26)} else {x = rep(FALSE, 26)}
  (b <- summary(lm(formula, data = FD_results %>% filter(Com == "No_Alesb"))))
  if(b$coefficients[2, 4] < 0.05){y = rep(TRUE, 26)} else {y = rep(FALSE, 26)}
  # vector assigning statistical significance for both curves
  (sig <- c(x, y))
  #local dataset of trait + significance
  data  = FD_results %>% select(paste0(var), Ni, Com) %>% cbind(sig)
  
  (plot <- ggplot(data, aes(x =Ni, y = data[,1], colour = Com))+
      geom_jitter(alpha = 0.5)+
      geom_smooth(aes(linetype = sig), method = "lm", se = FALSE) +
      scale_linetype_manual(values=c("TRUE"= "solid", "FALSE"="dashed")) +
      theme_classic() +
      xlab(label = "Soil Ni content (mg kg-1)"))+
    ylab(label = paste0(name_var)) +
   # ylim(-4, 6) +
    theme(legend.position = "none")+
    geom_hline(yintercept= 0, lty = 1, lwd = 1, alpha = 0.3)
}

graph_LA_fdis <- graph_fdis("FDis.LeafArea", name_var = paste("Leaf area"))
graph_SLA_fdis <- graph_fdis("FDis.SLA", name_var = paste("Specific leaf area"))
graph_LDMC_fdis <- graph_fdis("FDis.LDMC", name_var = paste("Leaf dry matter content"))
graph_LT_fdis <- graph_fdis("FDis.LT", name_var = paste("Leaf thickness"))
graph_LNC_fdis <- graph_fdis("FDis.LNC", name_var = paste("Leaf nitrogen content"))
graph_LS_fdis <- graph_fdis("FDis.LS", name_var = paste("Leaf shape"))

ggarrange(graph_LA_fdis, graph_SLA_fdis, graph_LDMC_fdis, graph_LT_fdis, 
          graph_LNC_fdis , graph_LS_fdis, ncol = 3, nrow = 2)
ggsave("FDis_uni_Aver_hel.pdf")

# FDis multivarate ---------------
formula = paste0("FDis.Tot ~ log(Ni)")
a <- summary(lm(formula, data = FD_results %>% filter(Com == "Tot")))
if(a$coefficients[2, 4] < 0.05){x = rep(TRUE, 26)} else {x = rep(FALSE, 26)}
b <- summary(lm(formula, data = FD_results %>% filter(Com == "No_Alesb")))
if(b$coefficients[2, 4] < 0.05){y = rep(TRUE, 26)} else {y = rep(FALSE, 26)}
# vector assigning statistical significance for both curves
(sig <- c(x, y))
#local dataset of trait + significance
data  = FD_results %>% select(FDis.Tot, Ni, Com) %>% cbind(sig)

ggplot(data, aes(x =Ni, y = FDis.Tot, colour = Com))+
  geom_jitter(height = 0.05, alpha = 0.5)+
  geom_smooth(aes(linetype = sig), method = "lm", se = FALSE) +
  scale_linetype_manual(values=c("TRUE"= "solid", "FALSE"="dashed")) +
  theme_classic() +
  xlab(label = "Soil Ni content (mg kg-1)")+
  ylab(label = paste0("sesFDis")) +
  ylim(-3, 4.3) +
  theme(legend.position = "none")+
   geom_hline(yintercept= 0, lty = 1, lwd = 1, alpha = 0.3)

ggsave("FDis_multi_ITV_hel.pdf")
