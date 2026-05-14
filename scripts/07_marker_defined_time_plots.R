#Something new 
b_data <- subset(data, functional_group_final == "B_cell")
b_markers <- c("CD19", "CD20", "IgD", "CD27", "CD38", "HLA-DR")
b_markers <- intersect(b_markers, colnames(b_data))

library(dplyr)

result_list <- list()

for (m in b_markers) {
  
  # define high-expression cells (top 50%)
  threshold <- median(b_data[[m]], na.rm = TRUE)
  
  subset_cells <- b_data[b_data[[m]] > threshold, ]
  
  summary <- subset_cells %>%
    group_by(timepoint) %>%
    summarise(median_gamma = median(yH2AX), .groups = "drop")
  
  summary$marker <- m
  
  result_list[[m]] <- summary
}

b_plot <- bind_rows(result_list)
b_plot$timepoint <- factor(b_plot$timepoint, levels = c(0, 30, 60, 360))

library(ggplot2)

ggplot(b_plot,
       aes(x = timepoint,
           y = median_gamma,
           color = marker,
           group = marker)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 2.5) +
  theme_minimal() +
  labs(
    title = "B cells: Y-H2AX response across marker-defined populations",
    x = "Timepoint",
    y = "Median γH2AX",
    color = "Marker"
  )

#marker based splits : T cells:
t_data <- subset(data, functional_group_final == "T_cell")

#t_markers <- c("CD3", "CD4", "CD8a", "CD45RA", "CD45RO", "CCR7", "CD27", "CD28", "CD161", "TCRgd")
t_markers <- c("CD4", "CD8a", "CD45RA", "CD45RO", "CCR7", "CD161")
t_markers <- intersect(t_markers, colnames(t_data))

library(dplyr)

result_list <- list()

for (m in t_markers) {
  
  # define high-expression cells (top 50%)
  threshold <- median(t_data[[m]], na.rm = TRUE)
  
  subset_cells <- t_data[t_data[[m]] > threshold, ]
  
  summary <- subset_cells %>%
    group_by(timepoint) %>%
    summarise(median_gamma = median(yH2AX), .groups = "drop")
  
  summary$marker <- m
  
  result_list[[m]] <- summary
}

t_plot <- bind_rows(result_list)
t_plot$timepoint <- factor(t_plot$timepoint, levels = c(0, 30, 60, 360))

library(ggplot2)

ggplot(t_plot,
       aes(x = timepoint,
           y = median_gamma,
           color = marker,
           group = marker)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 1.5) +
  theme_minimal() +
  labs(
    title = "T cells: Y-H2AX response across marker-defined populations",
    x = "Timepoint",
    y = "Median γH2AX",
    color = "Marker"
  )

#Myeloid cells
# 1. Subset Myeloid cells
my_data <- subset(data, functional_group_final == "Myeloid")

# 2. Define Myeloid markers
my_markers <- c("CD14", "CD11c", "CD123", "CD16", "HLA-DR", "CD38", "CCR7")
my_markers <- intersect(my_markers, colnames(my_data))

library(dplyr)

result_list <- list()

# 3. Loop through markers
for (m in my_markers) {
  
  # Define high-expression cells (top 50%)
  threshold <- median(my_data[[m]], na.rm = TRUE)
  
  subset_cells <- my_data[my_data[[m]] > threshold, ]
  
  summary <- subset_cells %>%
    group_by(timepoint) %>%
    summarise(median_gamma = median(yH2AX), .groups = "drop")
  
  summary$marker <- m
  
  result_list[[m]] <- summary
}

# 4. Combine results
my_plot <- bind_rows(result_list)

# 5. Fix timepoint order
my_plot$timepoint <- factor(my_plot$timepoint, levels = c(0, 30, 60, 360))

# 6. Plot
library(ggplot2)

ggplot(my_plot,
       aes(x = timepoint,
           y = median_gamma,
           color = marker,
           group = marker)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  theme_minimal() +
  labs(
    title = "Myeloid cells: Y-H2AX response across marker-defined populations",
    x = "Timepoint",
    y = "Median γH2AX",
    color = "Marker"
  )

#NK cells
# 1. Subset NK cells
nk_data <- subset(data, functional_group_final == "NK_cell")

# 2. Define NK markers
#nk_markers <- c("CD56", "CD16", "NKG2A", "CD57", "CD38", "CD45RA", "CD45RO")
nk_markers <- c("CD56", "CD16", "CD57")
nk_markers <- intersect(nk_markers, colnames(nk_data))

library(dplyr)

result_list <- list()

# 3. Loop through markers
for (m in nk_markers) {
  
  # Define high-expression cells (top 50%)
  threshold <- median(nk_data[[m]], na.rm = TRUE)
  
  subset_cells <- nk_data[nk_data[[m]] > threshold, ]
  
  summary <- subset_cells %>%
    group_by(timepoint) %>%
    summarise(median_gamma = median(yH2AX), .groups = "drop")
  
  summary$marker <- m
  
  result_list[[m]] <- summary
}

# 4. Combine results
nk_plot <- bind_rows(result_list)

# 5. Fix timepoint order
nk_plot$timepoint <- factor(nk_plot$timepoint, levels = c(0, 30, 60, 360))

# 6. Plot
library(ggplot2)

ggplot(nk_plot,
       aes(x = timepoint,
           y = median_gamma,
           color = marker,
           group = marker)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2.0) +
  theme_minimal() +
  labs(
    title = "NK cells: Y-H2AX response across marker-defined populations",
    x = "Timepoint",
    y = "Median γH2AX",
    color = "Marker"
  )


