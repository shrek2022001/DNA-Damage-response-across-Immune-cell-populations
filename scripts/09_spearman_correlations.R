#Correlation for B cells
b_data <- subset(data, functional_group_final == "B_cell")

b_markers <- c("CD19", "CD20", "IgD", "CD27", "CD38", "HLA-DR")
b_markers <- intersect(b_markers, colnames(b_data))

cor_results <- sapply(b_markers, function(m) {
  cor(b_data[[m]], b_data$yH2AX, use = "complete.obs", method = "spearman")
})

cor_results

cor_df <- data.frame(
  marker = names(cor_results),
  correlation = cor_results
)

cor_df <- cor_df[order(-abs(cor_df$correlation)), ]
cor_df

library(ggplot2)

ggplot(cor_df, aes(x = reorder(marker, correlation), y = correlation)) +
  geom_bar(stat = "identity") +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Correlation of markers with Y-H2AX (B cells)",
    x = "Marker",
    y = "Spearman correlation"
  )

# T cells 
t_data <- subset(data, functional_group_final == "T_cell")
t_markers <- c(
  "CD3", "CD4", "CD8a", 
  "CD45RA", "CD45RO",
  "CCR7", "CD27", "CD28",
  "CD161", "TCRgd"
)

# keep only existing ones
t_markers <- intersect(t_markers, colnames(t_data))
cor_results_t <- sapply(t_markers, function(m) {
  cor(t_data[[m]], t_data$yH2AX,
      method = "spearman",
      use = "complete.obs")
})
library(dplyr)

t_cor_df <- data.frame(
  marker = names(cor_results_t),
  correlation = cor_results_t
)

# sort for nicer plot
t_cor_df <- t_cor_df %>%
  arrange(correlation)
library(ggplot2)

ggplot(t_cor_df,
       aes(x = correlation,
           y = reorder(marker, correlation))) +
  geom_bar(stat = "identity", fill = "grey40") +
  theme_minimal() +
  labs(
    title = "Correlation of markers with Y-H2AX (T cells)",
    x = "Spearman correlation",
    y = "Marker"
  )

#Myeloid cells
my_data <- subset(data, functional_group_final == "Myeloid")
my_markers <- c(
  "CD14", "CD11c", "CD123",
  "CD16", "HLA-DR",
  "CD38", "CD45RA", "CD45RO",
  "CCR7", "CD27"
)

# keep only existing
my_markers <- intersect(my_markers, colnames(my_data))
cor_results_my <- sapply(my_markers, function(m) {
  cor(my_data[[m]], my_data$yH2AX,
      method = "spearman",
      use = "complete.obs")
})
library(dplyr)

my_cor_df <- data.frame(
  marker = names(cor_results_my),
  correlation = cor_results_my
)

my_cor_df <- my_cor_df %>%
  arrange(correlation)
library(ggplot2)

ggplot(my_cor_df,
       aes(x = correlation,
           y = reorder(marker, correlation))) +
  geom_bar(stat = "identity", fill = "grey40") +
  theme_minimal() +
  labs(
    title = "Correlation of markers with Y-H2AX (Myeloid cells)",
    x = "Spearman correlation",
    y = "Marker"
  )

#NK cells
nk_data <- subset(data, functional_group_final == "NK_cell")
nk_markers <- c(
  "CD56", "CD16", "NKG2A", "CD57",
  "CD38", "HLA-DR",
  "CD45RA", "CD45RO"
)

nk_markers <- intersect(nk_markers, colnames(nk_data))
cor_results_nk <- sapply(nk_markers, function(m) {
  cor(nk_data[[m]], nk_data$yH2AX,
      method = "spearman",
      use = "complete.obs")
})
library(dplyr)

nk_cor_df <- data.frame(
  marker = names(cor_results_nk),
  correlation = cor_results_nk
)

nk_cor_df <- nk_cor_df %>%
  arrange(correlation)
library(ggplot2)

ggplot(nk_cor_df,
       aes(x = correlation,
           y = reorder(marker, correlation))) +
  geom_bar(stat = "identity", fill = "grey40") +
  theme_minimal() +
  labs(
    title = "Correlation of markers with Y-H2AX (NK cells)",
    x = "Spearman correlation",
    y = "Marker"
  )


