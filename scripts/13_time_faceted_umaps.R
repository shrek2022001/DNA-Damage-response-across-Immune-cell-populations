#----Faceted UMAP over time---------
ggplot(t_data, aes(x = UMAP1, y = UMAP2, color = yH2AX)) +
  geom_point(size = 0.5) +
  facet_wrap(~timepoint) +
  scale_color_gradient(
    low = "grey90",
    high = "red"
  ) +
  theme_minimal() +
  labs(
    title = "DDR dynamics across time (T cells)",
    color = "γH2AX"
  )

t_data$high_DDR <- t_data$yH2AX > quantile(t_data$yH2AX, 0.9)

ggplot(t_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(
    data = subset(t_data, high_DDR),
    color = "red",
    size = 0.6
  ) +
  facet_wrap(~timepoint) +
  theme_minimal() +
  labs(title = "High DDR cells over time (T cells)")

#B cells over time DDR response
b_data$high_DDR <- b_data$yH2AX > quantile(b_data$yH2AX, 0.9)

ggplot(b_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(
    data = subset(b_data, high_DDR),
    color = "red",
    size = 0.6
  ) +
  facet_wrap(~timepoint) +
  theme_minimal() +
  labs(title = "High DDR cells over time (B cells)")

library(dplyr)

b_time_summary <- b_data %>%
  group_by(timepoint) %>%
  summarise(
    mean_gamma = mean(yH2AX),
    mean_IgD = mean(IgD),
    mean_CD27 = mean(CD27)
  )

b_time_summary
ggplot(b_data, aes(x = IgD, y = yH2AX, color = factor(timepoint))) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "lm") +
  theme_minimal()

#---DDR for T cells-----
nk_data$high_DDR <- nk_data$yH2AX > quantile(nk_data$yH2AX, 0.9)

ggplot(nk_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(
    data = subset(nk_data, high_DDR),
    color = "red",
    size = 0.6
  ) +
  facet_wrap(~timepoint) +
  theme_minimal() +
  labs(title = "High DDR cells over time (NK cells)")

#--Myeloid cells
my_data$high_DDR <- my_data$yH2AX > quantile(my_data$yH2AX, 0.9)

ggplot(my_data, aes(UMAP1, UMAP2)) +
  geom_point(color = "grey85", size = 0.4) +
  geom_point(
    data = subset(my_data, high_DDR),
    color = "red",
    size = 0.6
  ) +
  facet_wrap(~timepoint) +
  theme_minimal() +
  labs(title = "High DDR cells over time (Myeloid)")
