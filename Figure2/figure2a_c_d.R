library(tidyverse)

#Figure2a_top

ggplot() +
  geom_point(data = B_LH_res_anno_final %>% filter(padj > 0.05 | 
                                                     abs(log2FoldChange) < 1), 
             aes(log2FoldChange, -log10(padj)), size = 0.1, color = "gray") +
  geom_point(data = B_LH_res_anno_final %>% filter(padj <= 0.05 & 
                                                     log2FoldChange > 1), 
             aes(log2FoldChange, -log10(padj)), size = 0.5, color = "#75C376") +
  geom_point(data = B_LH_res_anno_final %>% filter(padj <= 0.05 & 
                                                     log2FoldChange < -1), 
             aes(log2FoldChange, -log10(padj)), size = 0.5, color = "darkgoldenrod2") +
  scale_x_continuous(limits = c(-10.5,10.5)) +
  theme_tufte() +
  theme(legend.position = "none")

#Figure2a_bottom

plot_filt <- B_LH_res_anno_final %>% filter(padj <= 0.05 & 
                                              abs(log2FoldChange) > 1)
y <- density(plot_filt$log2FoldChange, n = 2^12)

figure3e_b <- 
  ggplot(data.frame(x = y$x, y = y$y), aes(x, y)) + geom_line(size = 1) + 
  geom_segment(aes(xend = x, yend = 0, colour = x)) +
  scale_x_continuous(limits = c(-10.5,10.5)) +
  theme_tufte() +
  scale_color_gradient2(low = 'darkgoldenrod2', mid = "white", high = '#75C376') +
  theme(legend.position = "none",
        axis.title.x = element_blank())

#Figure2c_top
  ggplot(fig2c_top) +
  geom_bar(aes(category, Percentage, fill = GRN), 
           stat = "identity", position = "dodge") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

#Figure2c_bottom
  ggplot(fig2c_bottom) +
  geom_bar(aes(category, Percentage, fill = GRN), 
           stat = "identity", position = "dodge") +
  theme_bw(base_size = 12, base_family = "Helvetica") +
  theme(legend.position = "none",
        panel.grid.minor = element_blank())

#Figure2d
ggplot() +
  geom_jitter(data = figure2d_data %>% 
                filter(GRN == FALSE), 
              aes(sig_peak_count, sig_peak_count_ns), color = "gray", 
              size = 0.1, position = pos) +
  geom_jitter(data = figure2d_data%>% 
                filter(GRN == TRUE), 
              aes(sig_peak_count, sig_peak_count_ns), color = "maroon",
              size = 1, position = pos) +
  labs(x = "Number of Differentially Accessible OCRs",
     y = "Number of He-Accelerated OCRs") +
  scale_x_continuous(breaks=seq(0,11,1)) +
  scale_y_continuous(breaks=seq(0,4,1)) +
  theme_bw()+
  theme(legend.position = "none")

