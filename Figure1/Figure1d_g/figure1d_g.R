library(tidyverse)
library(ggpubr)
library(ggrepel)
library(ggthemes)
source("themes.R")
library(extrafont)
library(extrafontdb)
library(Rttf2pt1)
library(ggpattern)

#Figure1D
coding_props_plot <- 
  ggplot(coding_props) +
  geom_bar(aes(Species, prop, fill = Type), stat = "identity", position = "dodge") +
  theme_pubr()

#Figure1e
busted_pvals_plot <- 
  ggplot() +
  geom_point(data = busted_pvals %>% filter(class == "non-significant"), 
             aes(-log10(Ht_pval2), -log10(He_pval2), shape = type), 
             size = 2.5, stroke = 0.1, color = "lightgray", alpha = 0.25) +
  geom_point(data = busted_pvals %>% filter(class == "He_selection" & type == " SC"), 
             aes(-log10(Ht_pval2), -log10(He_pval2)), 
             size = 2.5, stroke = 0.25, color = "#F68D3F", alpha = 0.25, shape = 16) +
  geom_point(data = busted_pvals %>% filter(class == "He_selection" & type == "GRN"), 
             aes(-log10(Ht_pval2), -log10(He_pval2)), 
             size = 4, stroke = 0.25, fill = "#F68D3F", alpha = 1, shape = 24) +
  geom_point(data = busted_pvals %>% filter(class == "Ht_selection" & type == " SC"), 
             aes(-log10(Ht_pval2), -log10(He_pval2)), 
             size = 2.5, stroke = 0.25, color = "#74C476", alpha = 0.25, shape = 16) +
  geom_point(data = busted_pvals %>% filter(class == "Ht_selection" & type == "GRN"), 
             aes(-log10(Ht_pval2), -log10(He_pval2)), 
             size = 4, stroke = 0.25, fill = "#74C476", alpha = 1, shape = 24) +
  geom_point(data = busted_pvals %>% filter(class == "both_selection"), 
             aes(-log10(Ht_pval2), -log10(He_pval2)), 
             size = 2.5, stroke = 0.25, color = "#7373FF", alpha = 0.25, shape = 16) +
  scale_shape_manual(values = c(16,24)) +
  theme_pubr() +
  theme(legend.position = "none")

#Figure1f
noncoding_props_plot <- 
  ggplot(noncoding_props) +
  geom_bar(aes(Species, prop, fill = Type), stat = "identity", position = "dodge") +
  theme_pubr()

#Figure1g
noncoding_selection_plot <- 
  ggplot() + 
  geom_jitter(data=figure1g_data %>% filter((fdr.Ht_corr >= 0.1 & fdr.He_corr >= 0.1) | 
                                           (fdr.He_corr < 0.1 & zeta.He < 1.5) |
                                           (fdr.Ht_corr < 0.1 & zeta.Ht < 1.5)) %>%  #nonsig non GRN
                filter(GRN == FALSE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 2.5, stroke = 0.25, color = "lightgray", alpha = 0.1, shape = 16, width = 0.03, height = 0.03) +
  geom_point(data=figure1g_data %>% filter((fdr.Ht_corr >= 0.1 & fdr.He_corr >= 0.1) | 
                                          (fdr.He_corr < 0.1 & zeta.He < 1.5) |
                                          (fdr.Ht_corr < 0.1 & zeta.Ht < 1.5)) %>%   # nonsig GRN
               filter(GRN == TRUE), 
             aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
             size = 2.5, stroke = 0.25, color = "lightgray", alpha = 0.25, shape = 17) +  
  geom_jitter(data=figure1g_data %>% filter(zeta.He >= 1.5 & fdr.He_corr <= 0.1) %>%  #sigboth nonGRN (no GRN results)
                filter(zeta.Ht >= 1.5 & fdr.Ht_corr <= 0.1) %>% 
                filter(GRN == FALSE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 2.5, shape = 16, color = "#7373FF", alpha = 0.25, stroke = 0.25, width = 0.03) +
  geom_jitter(data=figure1g_data %>% filter(zeta.He >= 1.5 & fdr.He_corr <= 0.1) %>% #sigHe nonGRN
                filter(fdr.Ht_corr >= 0.1) %>% 
                filter(GRN == FALSE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 2.5, shape = 16, color = "#FD8D3C", alpha = 0.25, stroke = 0.25, width = 0.03) +
  geom_jitter(data=figure1g_data %>% filter(zeta.He >= 1.5 & fdr.He_corr <= 0.1) %>%#sigHe GRN
                filter(fdr.Ht_corr >= 0.1) %>% 
                filter(GRN == TRUE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 4, shape = 24, fill = "#FD8D3C", color = "black", stroke = 0.25, width = 0.03) +
  geom_jitter(data=figure1g_data %>% filter(zeta.Ht >= 1.5 & fdr.Ht_corr <= 0.1) %>% #sigHt nonGRN
                filter(fdr.He_corr >= 0.1) %>% 
                filter(GRN == FALSE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 2.5, shape = 16, color = "#74C476", alpha = 0.25, stroke = 0.25, width = 0.03) +
  geom_jitter(data=figure1g_data %>% filter(zeta.Ht >= 1.5 & fdr.Ht_corr <= 0.1) %>% #sigHt GRN
                filter(fdr.He_corr >= 0.1) %>% 
                filter(GRN == TRUE), 
              aes(-log10(fdr.Ht_corr), -log10(fdr.He_corr)), 
              size = 4, shape = 24, fill = "#74C476", stroke = 0.25, color = "black", width = 0.03) +
  theme_pubr()

  