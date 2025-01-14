## ----VENN DIAGRAM-------------------------------------------------------------

library(ggvenn)
library(viridis) 

#Import the files
heat <- read.csv("./filtered_DEGs2/Filtered_DEG2_Heat_vs_Control.csv", row.names = 1)
salinity <- read.csv("./filtered_DEGs2/Filtered_DEG2_Salinity_vs_Control.csv", row.names = 1)
salinity_heat <- read.csv("./filtered_DEGs2/Filtered_DEG2_S+H_vs_Control.csv", row.names = 1)

 
# Upregulated DEGs
heat_up<- rownames(subset(heat, log2FoldChange > 1))   
salinity_up<- rownames(subset(salinity, log2FoldChange > 1))           
salinity_heat_up  <- rownames(subset(salinity_heat, log2FoldChange > 1) )              

# Create the Venn diagram of the upregulated genes
venn_up <- list(Heat = heat_up, Salinity = salinity_up, "Salinity+Heat" = salinity_heat_up )

ggvenn(
  venn_up, 
  fill_color = c("#F89C74", "#66C5CC", "#F6CF71"),
  stroke_size = 0.5, set_name_size = 5, text_size = 4 ) +
  labs(title = "Up-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))



# Downregulated DEGs
heat_down<- rownames(subset(heat, log2FoldChange < 1))   
salinity_down<- rownames(subset(salinity, log2FoldChange < 1))           
salinity_heat_down  <- rownames(subset(salinity_heat, log2FoldChange < 1) )              

# Create the Venn diagram of the upregulated genes
venn_down <- list(Heat = heat_down, Salinity = salinity_down, "Salinity+Heat" = salinity_heat_down )

ggvenn(
  venn_down, 
  fill_color = c("#F89C74", "#66C5CC", "#F6CF71"),
  stroke_size = 0.5, set_name_size = 5, text_size = 4) +
  labs(title = "Down-regulated genes") +
  theme(plot.title = element_text(hjust = 0.5, size = 18, face = "bold"))




