#UMAP for NK cells
# 1. Subset NK cells
nk_data <- subset(data, functional_group_final == "NK_cell")
nk_data <- as.data.frame(nk_data)

# 2. Choose NK phenotype markers
nk_features <- c(
  "CD56",    # NK identity / CD56bright vs dim
  "CD16",    # cytotoxic NK
  "NKG2A",   # inhibitory / immature-like
  "CD57",    # mature / senescent NK
  "CD38",    # activation
  "HLA-DR",  # activation / antigen presentation-like
  "CD45RA",
  "CD45RO"
)

nk_features <- intersect(nk_features, colnames(nk_data))

# 3. Scale
nk_scaled <- scale(nk_data[, nk_features])
# Elbow plot for NK cells
wss <- sapply(2:10, function(k) {
  kmeans(nk_scaled, centers = k, nstart = 10, iter.max = 50)$tot.withinss
})

plot(2:10, wss, type = "b", pch = 19,
     xlab = "k", ylab = "Within-cluster sum of squares",
     main = "NK cells: elbow plot")
set.seed(123)
nk_k <- 5
nk_data$cluster <- kmeans(nk_scaled, centers = nk_k, nstart = 25, iter.max = 50)$cluster
nk_data$cluster <- as.factor(nk_data$cluster)

table(nk_data$cluster)
set.seed(123)
nk_km1 <- kmeans(nk_scaled, centers = nk_k, nstart = 25, iter.max = 50)

set.seed(456)
nk_km2 <- kmeans(nk_scaled, centers = nk_k, nstart = 25, iter.max = 50)

table(nk_km1$cluster, nk_km2$cluster)
library(uwot)

set.seed(123)
nk_umap <- umap(nk_scaled)

nk_data$UMAP1 <- as.numeric(nk_umap[,1])
nk_data$UMAP2 <- as.numeric(nk_umap[,2])
library(ggplot2)

ggplot(nk_data, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(alpha = 0.6, size = 0.7) +
  theme_minimal() +
  labs(title = "NK cell clusters (UMAP)")
ggplot(nk_data, aes(UMAP1, UMAP2, color = yH2AX)) +
  geom_point(alpha = 0.6, size = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "DDR landscape in NK cells")

library(dplyr)

nk_data$high_DDR <- nk_data$yH2AX > quantile(nk_data$yH2AX, 0.9)

nk_cluster_enrichment <- nk_data %>%
  group_by(cluster) %>%
  summarise(
    total_cells = n(),
    high_ddr_cells = sum(high_DDR),
    fraction_high = high_ddr_cells / total_cells,
    median_gamma = median(yH2AX),
    .groups = "drop"
  )

nk_cluster_enrichment
aggregate(nk_data[, nk_features],
          by = list(cluster = nk_data$cluster),
          median)

cor(nk_data$yH2AX, nk_data$CD56, method = "spearman")
cor(nk_data$yH2AX, nk_data$CD16, method = "spearman")
cor(nk_data$yH2AX, nk_data$CD57, method = "spearman")
cor(nk_data$yH2AX, nk_data$NKG2A, method = "spearman")
cor(nk_data$yH2AX, nk_data$HLA-DR, method = "spearman")

