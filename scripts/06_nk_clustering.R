#------NK cells -------------------------

nk_data <- subset(data, functional_group_final == "NK_cell")

dim(nk_data)

nk_features <- c(
  "CD56", "CD16", "NKG2A", "CD57",
  "CD38", "HLA-DR", "CD45RA", "CD45RO"
)

nk_features <- intersect(nk_features, colnames(nk_data))
nk_features

set.seed(123)

nk_scaled <- scale(nk_data[, nk_features])

k <- 5
nk_data$cluster <- kmeans(nk_scaled, centers = k, nstart = 10)$cluster

table(nk_data$cluster)

nk_cluster_summary <- aggregate(
  nk_data[, nk_features],
  by = list(cluster = nk_data$cluster),
  median
)

nk_cluster_summary

gamma_nk <- aggregate(
  nk_data$yH2AX,
  by = list(
    cluster = nk_data$cluster,
    timepoint = nk_data$timepoint
  ),
  median
)

colnames(gamma_nk)[3] <- "median_gamma"
gamma_nk$timepoint <- factor(gamma_nk$timepoint,
                             levels = c(0, 30, 60, 360))
ggplot(gamma_nk,
       aes(x = timepoint, y = median_gamma,
           color = factor(cluster),
           group = cluster)) +
  geom_line() +
  geom_point() +
  theme_minimal()

