library(tidyverse)

#Figure 4d
pmar_similarity_pairwise <- 
  ggplot(figure4d_data, aes(Species, Similarity, fill = Sequence)) +
  geom_boxplot(outlier.size = .75) +
  labs(y = "Ortholog Pairwise Similarity (%)",
       title = "Pmar1") +
  theme_pubclean() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(size = 12))