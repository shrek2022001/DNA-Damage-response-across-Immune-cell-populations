#---------Let's work on just B-cells-----------
b_data <- subset(data, functional_group_final == "B_cell")

dim(b_data)

#Choosing markers related to B-cells
b_markers <- c("CD19", "CD20", "IgD", "CD38", "HLA-DR", "CD27", "CD45RA", "CD45RO", "Ki67", "yH2AX")
b_markers <- intersect(b_markers, colnames(b_data))
b_markers

# compare across time
aggregate(
  b_data[, b_markers],
  by = list(timepoint = b_data$timepoint),
  FUN = median
)

#Compare across subsets
aggregate(
  b_data[, b_markers],
  by = list(subset = b_data$subset),
  FUN = median
)

#Graph 
ggplot(b_data, aes(x = IgD, y = yH2AX)) +
  geom_point(alpha = 0.25, size = 0.5) +
  facet_wrap(~timepoint) +
  theme_minimal()
ggplot(b_data, aes(x = IgD, y = yH2AX)) +
  geom_point(alpha = 0.2, size = 0.4) +
  geom_smooth(method = "lm", color = "red") +
  facet_wrap(~timepoint) +
  theme_minimal()

#Let's test the trands with CD27
ggplot(b_data, aes(x = CD27, y = yH2AX)) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = "lm") +
  facet_wrap(~timepoint) +
  theme_minimal()

ggplot(b_data, aes(x = IgD, y = CD27, color = yH2AX)) +
  geom_point(alpha = 0.4, size = 0.7) +
  facet_wrap(~timepoint) +
  theme_minimal()

