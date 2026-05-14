#----Myeloid cells------------
my_data <- subset(data, functional_group_final == "Myeloid")
my_data <- as.data.frame(my_data)

my_features <- c(
  "CD14", "CD11c", "CD123", "CD16",
  "HLA-DR", "CD38", "CD45RA", "CD45RO",
  "CCR7", "CD27"
)

my_features <- intersect(my_features, colnames(my_data))

my_scaled <- scale(my_data[, my_features])

library(uwot)

set.seed(123)
my_umap <- umap(my_scaled)

my_data$UMAP1 <- as.numeric(my_umap[,1])
my_data$UMAP2 <- as.numeric(my_umap[,2])
wss <- sapply(2:10, function(k) {
  kmeans(my_scaled, centers = k, nstart = 25, iter.max = 50)$tot.withinss
})

plot(2:10, wss, type = "b", pch = 19,
     xlab = "k", ylab = "Within-cluster sum of squares",
     main = "Myeloid cells: elbow plot")
my_k <- 5

set.seed(123)
my_data$cluster <- kmeans(my_scaled, centers = my_k, nstart = 25, iter.max = 50)$cluster
my_data$cluster <- as.factor(my_data$cluster)
library(ggplot2)

ggplot(my_data, aes(UMAP1, UMAP2, color = cluster)) +
  geom_point(alpha = 0.6, size = 0.7) +
  theme_minimal() +
  labs(title = "Myeloid clusters (UMAP)")
ggplot(my_data, aes(UMAP1, UMAP2, color = yH2AX)) +
  geom_point(alpha = 0.6, size = 0.7) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "DDR landscape in Myeloid cells")

ggplot(t_data, aes(UMAP1, UMAP2, color = yH2AX)) +
  geom_point(size = 0.5) +
  facet_wrap(~timepoint) +
  scale_color_viridis_c() +
  theme_minimal() +
  labs(title = "DDR over time in T cells")

ggplot(t_data, aes(UMAP1, UMAP2, color = yH2AX)) +
  geom_point(size = 0.5) +
  facet_wrap(~timepoint) +
  scale_color_gradient(
    low = "grey85",
    high = "red"
  ) +
  theme_minimal() +
  labs(title = "DDR dynamics over time (T cells)")

