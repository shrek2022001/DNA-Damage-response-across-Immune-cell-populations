#Myeloid cells
my_data <- subset(data, functional_group_final == "Myeloid")

dim(my_data)
my_features <- c(
  "CD14", "CD11c", "CD123", "CD16",
  "HLA-DR", "CD38", "CD45RA", "CD45RO",
  "CCR7", "CD27"  # optional context markers
)

my_features <- intersect(my_features, colnames(my_data))
my_features

set.seed(123)

my_scaled <- scale(my_data[, my_features])

k <- 5  
my_data$cluster <- kmeans(my_scaled, centers = k, nstart = 10)$cluster

table(my_data$cluster)
my_cluster_summary <- aggregate(
  my_data[, my_features],
  by = list(cluster = my_data$cluster),
  median
)

my_cluster_summary

gamma_my <- aggregate(
  my_data$yH2AX,
  by = list(cluster = my_data$cluster,
            timepoint = my_data$timepoint),
  median
)

colnames(gamma_my)[3] <- "median_gamma"
gamma_my$timepoint <- factor(gamma_my$timepoint,
                             levels = c(0, 30, 60, 360))
ggplot(gamma_my,
       aes(x = timepoint, y = median_gamma,
           color = factor(cluster),
           group = cluster)) +
  geom_line() +
  geom_point() +
  theme_minimal()

