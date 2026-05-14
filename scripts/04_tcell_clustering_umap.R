#For T cells
t_data <- subset(data, functional_group_final == "T_cell")
t_features <- c(
  "CD4", "CD8a", "CD45RA", "CD45RO",
  "CCR7", "CD27", "CD28",
  "CD161", "TCRgd"
)

t_features <- intersect(t_features, colnames(t_data))
set.seed(123)

t_scaled <- scale(t_data[, t_features])

wss <- sapply(2:10, function(k){
  kmeans(t_scaled, centers = k, nstart = 10)$tot.withinss
})

plot(2:10, wss, type = "b", pch = 19,
     xlab = "k", ylab = "Within-cluster sum of squares")

# Run clustering with k=5
set.seed(123)
t_data$cluster <- kmeans(t_scaled, centers = 5)$cluster

table(t_data$cluster)
cluster_summary <- aggregate(
  t_data[, t_features],
  by = list(cluster = t_data$cluster),
  median
)

cluster_summary

gamma_cluster <- aggregate(
  t_data$yH2AX,
  by = list(cluster = t_data$cluster,
            timepoint = t_data$timepoint),
  median
)

colnames(gamma_cluster)[3] <- "median_gamma"
data$timepoint <- factor(data$timepoint,
                         levels = c(0, 30, 60, 360))
ggplot(gamma_cluster,
       aes(x = timepoint, y = median_gamma,
           color = factor(cluster),
           group = cluster)) +
  geom_line() +
  geom_point() +
  theme_minimal()

library(uwot)

t_data <- subset(data, functional_group_final == "T_cell")

t_features <- c(
  "CD4", "CD8a", "CD45RA", "CD45RO",
  "CCR7", "CD27", "CD28",
  "CD161", "TCRgd"
)

t_features <- intersect(t_features, colnames(t_data))

t_scaled <- scale(t_data[, t_features])

set.seed(123)
t_data$cluster <- kmeans(t_scaled, centers = 5, nstart = 10)$cluster
# Use same features as clustering
library(uwot)

t_scaled <- scale(t_data[, t_features])

set.seed(123)
t_umap <- umap(t_scaled)

# IMPORTANT: force numeric vectors
t_data$UMAP1 <- as.numeric(t_umap[,1])
t_data$UMAP2 <- as.numeric(t_umap[,2])
str(t_data$cluster)
table(t_data$cluster)
library(uwot)

set.seed(123)
t_umap <- umap(t_scaled)

t_data$UMAP1 <- as.numeric(t_umap[,1])
t_data$UMAP2 <- as.numeric(t_umap[,2])
t_data$cluster <- as.factor(t_data$cluster)

library(ggplot2)

ggplot(t_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(alpha = 0.5, size = 0.6) +
  theme_minimal() +
  labs(
    color = "Cluster",
    title = "T cell clusters (UMAP)"
  )
ggplot(t_data, aes(x = UMAP1, y = UMAP2, color = yH2AX)) +
  geom_point(alpha = 0.6, size = 0.6) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "DDR landscape (γH2AX) in T cells")
t_data$high_DDR <- t_data$yH2AX > quantile(t_data$yH2AX, 0.9)

ggplot(t_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey80", size = 0.4) +
  geom_point(data = subset(t_data, high_DDR),
             color = "red", size = 0.6) +
  theme_minimal() +
  labs(title = "High DDR cells (top 10%)")

#Cluster level DDR instead of cells:
library(dplyr)

cluster_ddr <- t_data %>%
  group_by(cluster) %>%
  summarise(mean_gamma = mean(yH2AX))

t_data <- merge(t_data, cluster_ddr, by = "cluster")

ggplot(t_data, aes(UMAP1, UMAP2, color = mean_gamma)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "Cluster-level DDR")

#identify which cluster is high DDR
aggregate(t_data$yH2AX,
          by = list(cluster = t_data$cluster),
          median)
cluster_enrichment
aggregate(t_data[, t_features],
          by = list(cluster = t_data$cluster),
          median)
annotate("text", x = 4, y = 5, label = "High DDR", size = 5)

