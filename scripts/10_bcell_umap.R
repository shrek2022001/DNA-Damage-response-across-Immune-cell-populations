#------------UMAPS------------------------------------------
b_data <- subset(data, functional_group_final == "B_cell")
b_data <- as.data.frame(b_data)

b_features <- c(
  "IgD",     # naive
  "CD27",    # memory
  "CD38",    # activation/plasmablast
  "CD20",    # mature B cells
  "HLA-DR"   # antigen presentation
)

b_features <- intersect(b_features, colnames(b_data))
b_scaled <- scale(b_data[, b_features])

set.seed(123)
b_data$cluster <- kmeans(b_scaled, centers = 4, nstart = 10)$cluster
b_data$cluster <- as.factor(b_data$cluster)
table(b_data$cluster)

library(uwot)

set.seed(123)
b_umap <- umap(b_scaled)

b_data$UMAP1 <- b_umap[,1]
b_data$UMAP2 <- b_umap[,2]

library(ggplot2)

ggplot(b_data, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(alpha = 0.6, size = 0.6) +
  theme_minimal() +
  labs(title = "B cell clusters (UMAP)")
ggplot(b_data, aes(UMAP1, UMAP2, color = yH2AX)) +
  geom_point(alpha = 0.6, size = 0.6) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "DDR landscape in B cells")
b_data$high_DDR <- b_data$yH2AX > quantile(b_data$yH2AX, 0.9)

ggplot(b_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(data = subset(b_data, high_DDR),
             color = "red", size = 0.6) +
  theme_minimal() +
  labs(title = "High DDR B cells")
cor(b_data$yH2AX, b_data$IgD, method = "spearman")
cor(b_data$yH2AX, b_data$CD27, method = "spearman")
cor(b_data$yH2AX, b_data$CD38, method = "spearman")

