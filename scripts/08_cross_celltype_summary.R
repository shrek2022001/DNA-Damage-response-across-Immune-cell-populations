#DDR 
peak_gamma <- data %>%
  group_by(functional_group_final, timepoint) %>%
  summarise(median_gamma = median(yH2AX), .groups = "drop") %>%
  group_by(functional_group_final) %>%
  summarise(peak_gamma = max(median_gamma))

peak_gamma

#Fold change from baseline
baseline <- data %>%
  filter(timepoint == 0) %>%
  group_by(functional_group_final) %>%
  summarise(baseline_gamma = median(yH2AX))

response <- data %>%
  group_by(functional_group_final, timepoint) %>%
  summarise(median_gamma = median(yH2AX), .groups = "drop")

fold_change <- merge(response, baseline, by = "functional_group_final")

fold_change$fc <- fold_change$median_gamma / fold_change$baseline_gamma

fold_summary <- fold_change %>%
  group_by(functional_group_final) %>%
  summarise(max_fc = max(fc))

fold_summary

#AUC
auc_results <- data %>%
  group_by(functional_group_final, timepoint) %>%
  summarise(median_gamma = median(yH2AX), .groups = "drop") %>%
  arrange(functional_group_final, timepoint) %>%
  group_by(functional_group_final) %>%
  summarise(
    auc = sum(diff(as.numeric(timepoint)) *
                (head(median_gamma, -1) + tail(median_gamma, -1)) / 2)
  )

auc_results

#Consolidated results
library(dplyr)

celltype_summary <- data %>%
  group_by(functional_group_final, timepoint) %>%
  summarise(
    median_gamma = median(yH2AX, na.rm = TRUE),
    .groups = "drop"
  )
celltype_summary$timepoint <- factor(
  celltype_summary$timepoint,
  levels = c(0, 30, 60, 360)
)

library(ggplot2)

ggplot(celltype_summary,
       aes(x = timepoint,
           y = median_gamma,
           color = functional_group_final,
           group = functional_group_final)) +
  geom_line(linewidth = 1.0) +
  geom_point(size = 2) +
  theme_minimal() +
  labs(
    title = "Y-H2AX response over time across immune cell types",
    x = "Timepoint",
    y = "Median γH2AX",
    color = "Cell Type"
  )

