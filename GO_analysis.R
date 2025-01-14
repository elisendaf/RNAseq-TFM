## ----GENE ONTOLOGY------------------------------------------------------------

library(data.table)
library(clusterProfiler)
library(GO.db)
library(ggplot2)
library(dplyr)



# 1. Load the annotation file and select only the important columns
annotation <- fread("./gene_association.itag4.0.nr", sep = "\t", header = F,
                     select = c(2, 5, 9))

colnames(annotation) <- c("GeneID", "GO", "Aspect")  # Name the columns


# 2. Load the DEGs files for all conditions
# HEAT 
deg_heat <- fread("./filtered_DEGs2/Filtered_DEG2_Heat_vs_Control.csv")
deg_heat <- sub("(\\.[^.]*)\\..*$", "\\1", deg_heat[[1]])  # Clean Solycs to be in the same format as in the annotation file

# SALINITY 
deg_salinity <- fread("./filtered_DEGs/Filtered_DEG_Salinity_vs_Control.csv")
deg_salinity <- sub("(\\.[^.]*)\\..*$", "\\1", deg_salinity[[1]])  # Clean Solycs

# HEAT + SALINITY
deg_sh <- fread("./filtered_DEGs/Filtered_DEG_S+H_vs_Control.csv")
deg_sh <- sub("(\\.[^.]*)\\..*$", "\\1", deg_sh[[1]])  # Clean Solycs



##### MOLECULAR FUNCTION--------------------------------------------------------

# 3. Filter the GO names by molecular function (F)
gene2go_mf <- annotation[annotation$Aspect == "F", .(TERM = GO, GENE = GeneID)]


## HEAT
# 4. Gene ontology enrichment using enricher function
go_heat <- enricher(
   gene = deg_heat,     #List of DEG gene IDs
   TERM2GENE = gene2go_mf, #Mapping between GO terms and genes
   pvalueCutoff = 0.05,     #Filter by padj <= 0.05
   pAdjustMethod = "fdr"    #Adjusts p-values using FDR (False Discovery Rate)
  )
  
# Map the GO IDs to the term descriptions
go_heat@result$Description <- mapIds(GO.db, 
                                     keys = go_heat@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

  
# Dotplot
dotplot(go_heat, showCategory = 10) +
  ggtitle("Heat") +
  ylab("Molecular Function")


### SALINITY
# Gene ontology enrichment using enricher function
go_salinity <- enricher(
  gene = deg_salinity, 
  TERM2GENE = gene2go_mf, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_salinity@result$Description <- mapIds(GO.db, 
                                     keys = go_salinity@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

# Dotplot
dotplot(go_salinity, showCategory = 20) +
  ggtitle("Salinity") +
  ylab("Molecular Function")



### SALINITY + HEAT
# Gene ontology enrichment using enricher function
go_sh <- enricher(
  gene = deg_sh, 
  TERM2GENE = gene2go_mf, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_sh@result$Description <- mapIds(GO.db, 
                                     keys = go_sh@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

# Dotplot
dotplot(go_sh, showCategory = 20) +
  ggtitle("Salinity + Heat") +
  ylab("Molecular Function")



##### BIOLOGICAL PROCESS--------------------------------------------------------

# Filter the GO names by biological process (P)
gene2go_bp <- annotation[annotation$Aspect == "P", .(TERM = GO, GENE = GeneID)]


### HEAT
# Gene ontology enrichment using enricher function
go_bp_heat <- enricher(
  gene = deg_heat,     #List of DEG gene IDs
  TERM2GENE = gene2go_bp, #Mapping between GO terms and genes
  pvalueCutoff = 0.05,     #Filter by padj <= 0.05
  pAdjustMethod = "fdr"    #Adjusts p-values using FDR (False Discovery Rate)
)

# Map the GO IDs to the term descriptions
go_bp_heat@result$Description <- mapIds(GO.db, 
                                     keys = go_bp_heat@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

# Dotplot
dotplot(go_bp_heat, showCategory = 10) +
  ggtitle("Heat") +
  ylab("Biological process")


barplot(go_bp_heat, drop = T, showCategory = 12)

### SALINITY
# Gene ontology enrichment using enricher function
go_bp_salinity <- enricher(
  gene = deg_salinity, 
  TERM2GENE = gene2go_bp, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_bp_salinity@result$Description <- mapIds(GO.db, 
                                     keys = go_bp_salinity@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

# Dotplot
dotplot(go_bp_salinity, showCategory = 20) +
  ggtitle("Salinity") +
  ylab("Biological process")



### SALINITY + HEAT
# Gene ontology enrichment using enricher function
go_bp_sh <- enricher(
  gene = deg_sh, 
  TERM2GENE = gene2go_bp, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_bp_sh@result$Description <- mapIds(GO.db, 
                                     keys = go_bp_sh@result$ID,
                                     column = "TERM",
                                     keytype = "GOID",
                                     multiVals = "first")

# Dotplot
dotplot(go_bp_sh, showCategory = 20) +
  ggtitle("Salinity + Heat") +
  ylab("Biological process")






### CELLULAR COMPARTMENT--------------------------------------------------------

# Filter the GO names by cellular compartment (C)
gene2go_cc <- annotation[annotation$Aspect == "C", .(TERM = GO, GENE = GeneID)]


### HEAT
# Gene ontology enrichment using enricher function
go_cc_heat <- enricher(
  gene = deg_heat,     #List of DEG gene IDs
  TERM2GENE = gene2go_cc, #Mapping between GO terms and genes
  pvalueCutoff = 0.05,     #Filter by padj <= 0.05
  pAdjustMethod = "fdr"    #Adjusts p-values using FDR (False Discovery Rate)
)

# Map the GO IDs to the term descriptions
go_cc_heat@result$Description <- mapIds(GO.db, 
                                        keys = go_cc_heat@result$ID,
                                        column = "TERM",
                                        keytype = "GOID",
                                        multiVals = "first")

# Dotplot
dotplot(go_cc_heat, showCategory = 10) +
  ggtitle("Heat") +
  ylab("Biological process")


### SALINITY
# Gene ontology enrichment using enricher function
go_cc_salinity <- enricher(
  gene = deg_salinity, 
  TERM2GENE = gene2go_cc, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_cc_salinity@result$Description <- mapIds(GO.db, 
                                            keys = go_cc_salinity@result$ID,
                                            column = "TERM",
                                            keytype = "GOID",
                                            multiVals = "first")

# Dotplot
dotplot(go_cc_salinity, showCategory = 20) +
  ggtitle("Salinity") +
  ylab("Biological process")

emapplot(go_cc_salinity)

### SALINITY + HEAT
# Gene ontology enrichment using enricher function
go_cc_sh <- enricher(
  gene = deg_sh, 
  TERM2GENE = gene2go_cc, 
  pvalueCutoff = 0.05, 
  pAdjustMethod = "fdr"
)

# Map the GO IDs to the term descriptions
go_cc_sh@result$Description <- mapIds(GO.db, 
                                      keys = go_cc_sh@result$ID,
                                      column = "TERM",
                                      keytype = "GOID",
                                      multiVals = "first")

# Dotplot
dotplot(go_cc_sh, showCategory = 20) +
  ggtitle("Salinity + Heat") +
  ylab("Biological process")




### COMBINED PLOTS -------------------------------------------------------------

## HEAT

# Extract results from the enrichments and add category labels
mf_results_h <- go_heat@result %>%
  mutate(Category = "MF") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

bp_results_h <- go_bp_heat@result %>%
  mutate(Category = "BP") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

# Combine MF and BP results
combined_results_h <- rbind(mf_results_h, bp_results_h)

# Filter for top results (e.g., top 10 from each category by p-value)
combined_results_h <- combined_results_h %>%
  group_by(Category) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()

# Plot the combined data
ggplot(combined_results_h, aes(x = Category, y = reorder(Description, -p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for p.adjust
  labs(
    title = " Heat GO Enrichment Analysis",
    x = NULL,
    y = NULL,
    size = "Gene Count",
    color = "p.adjust"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


## SALINITY
# Extract results from the enrichments and add category labels
mf_results_s <- go_salinity@result %>%
  mutate(Category = "MF") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

bp_results_s <- go_bp_salinity@result %>%
  mutate(Category = "BP") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

# Combine MF and BP results
combined_results_s <- rbind(mf_results_s, bp_results_s)

# Filter for top results (e.g., top 10 from each category by p-value)
combined_results_s <- combined_results_s %>%
  group_by(Category) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()

# Plot the combined data
ggplot(combined_results_s, aes(x = Category, y = reorder(Description, -p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for p.adjust
  labs(
    title = "Salinity GO Enrichment Analysis",
    x = NULL,
    y = NULL,
    size = "Gene Count",
    color = "p.adjust"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )


## SALINITY + HEAT
# Extract results from the enrichments and add category labels
mf_results_sh <- go_sh@result %>%
  mutate(Category = "MF") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

bp_results_sh <- go_bp_sh@result %>%
  mutate(Category = "BP") %>%
  dplyr::select(ID, Description, p.adjust, Count, Category)

# Combine MF and BP results
combined_results_sh <- rbind(mf_results_sh, bp_results_sh)

# Filter for top results (e.g., top 10 from each category by p-value)
combined_results_sh <- combined_results_sh %>%
  group_by(Category) %>%
  slice_min(order_by = p.adjust, n = 10) %>%
  ungroup()

# Plot the combined data
ggplot(combined_results_sh, aes(x = Category, y = reorder(Description, -p.adjust), size = Count, color = p.adjust)) +
  geom_point(alpha = 0.8) +
  scale_color_gradient(low = "red", high = "blue") +  # Color gradient for p.adjust
  labs(
    title = "Salinity + heat GO Enrichment Analysis",
    x = NULL,
    y = NULL,
    size = "Gene Count",
    color = "p.adjust"
  ) +
  theme_minimal() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5)
  )













