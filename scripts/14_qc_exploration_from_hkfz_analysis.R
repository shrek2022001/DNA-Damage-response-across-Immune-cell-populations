setwd("C:/Users/gaura/OneDrive/Documents/shreya/hkfz_data")

install.packages("BiocManager")
BiocManager::install("flowCore")

library(flowCore)

files <- c(
  "c01_DDR_Processed (1).fcs",
  "c02_DDR_Processed (1).fcs",
  "c03_DDR_Processed (1).fcs",
  "c04_DDR_Processed (1).fcs"
)

times <- c(0, 30, 60, 360)

# Load all FCS files
fcs_list <- lapply(files, function(f) {
  read.FCS(f, transformation = FALSE)
})

# Check event counts per file
event_counts <- sapply(fcs_list, function(f) {
  nrow(exprs(f))
})

event_counts

# Build QC-preserved combined dataset
all_data_qc <- list()

for (i in seq_along(files)) {
  f <- fcs_list[[i]]
  df <- as.data.frame(exprs(f))
  
  # Map marker names
  desc <- pData(parameters(f))
  new_names <- desc$desc
  idx <- is.na(new_names) | new_names == ""
  new_names[idx] <- desc$name[idx]
  colnames(df) <- new_names
  
  # Keep QC-relevant and biological columns, remove obvious junk
  keep_cols_qc <- colnames(df)[
    !grepl(
      "Time|Center|Width|Residual|Offset|Amplitude|Barcode|Bead|BCKG|IdU|208Pb",
      colnames(df),
      ignore.case = TRUE
    )
  ]
  
  df <- df[, keep_cols_qc]
  df$time <- times[i]
  
  all_data_qc[[i]] <- df
}

combined_qc <- do.call(rbind, all_data_qc)

dim(combined_qc)
colnames(combined_qc)
head(combined_qc)

library(ggplot2)

ggplot(combined_qc, aes(x = `191Ir_DNA1`)) +
  geom_histogram(bins = 100) +
  theme_bw()
# amuk amuk intensity che kiti cells ahet? kami intensity= dead cells
ggplot(combined_qc, aes(x = `191Ir_DNA1`, y = `193Ir_DNA2`)) +
  geom_point(alpha = 0.3, size = 0.5)

dna1 <- combined_qc$`191Ir_DNA1`
dna2 <- combined_qc$`193Ir_DNA2`

# Example thresholds (we'll refine)
keep <- dna1 > 500 & dna1 < 5000 & dna2 > 500 & dna2 < 5000

filtered <- combined_qc[keep, ]

dna1 <- combined_qc$`191Ir_DNA1`
dna2 <- combined_qc$`193Ir_DNA2`
ld   <- combined_qc$`103Rh_Live_Dead`
el   <- combined_qc$Event_length

dna_keep <- dna1 > quantile(dna1, 0.05) & dna1 < quantile(dna1, 0.99) &
  dna2 > quantile(dna2, 0.05) & dna2 < quantile(dna2, 0.99)

ld_keep <- ld < quantile(ld, 0.95)

el_keep <- el > quantile(el, 0.01) & el < quantile(el, 0.99)

keep_all <- dna_keep & ld_keep & el_keep

filtered <- combined_qc[keep_all, ]
dim(filtered)

n_before <- nrow(combined_qc)
n_after  <- nrow(filtered)
pct_removed <- 100 * (n_before - n_after) / n_before
pct_removed

library(ggplot2)

ggplot(filtered, aes(x = `191Ir_DNA1`, y = `193Ir_DNA2`)) +
  geom_point(alpha = 0.3, size = 0.5) +
  theme_bw()

ggplot(filtered, aes(x = `103Rh_Live_Dead`)) +
  geom_histogram(bins = 100) +
  theme_bw()

ggplot(filtered, aes(x = Event_length)) +
  geom_histogram(bins = 100) +
  theme_bw()

cor(combined_qc$`191Ir_DNA1`, combined_qc$`193Ir_DNA2`)
cor(filtered$`191Ir_DNA1`, filtered$`193Ir_DNA2`)

ggplot(analysis_data, aes(x = factor(time), y = `165Ho_yH2AX`)) +
  geom_violin(trim = FALSE) +
  theme_bw()

analysis_cols <- colnames(filtered)[
  !grepl(
    "Event_length|Live_Dead|DNA1|DNA2|194Pt",
    colnames(filtered),
    ignore.case = TRUE
  )
]

analysis_data <- filtered[, analysis_cols]

marker_cols <- setdiff(colnames(analysis_data), "time")
analysis_data[, marker_cols] <- asinh(analysis_data[, marker_cols] / 5)

ggplot(analysis_data, aes(x = factor(time), y = `165Ho_yH2AX`)) +
  geom_violin(trim = FALSE) +
  theme_bw()

library(dplyr)

analysis_data %>%
  group_by(time) %>%
  summarise(median_yH2AX = median(`165Ho_yH2AX`))

ggplot(analysis_data, aes(x = factor(time), y = `165Ho_yH2AX`)) +
  geom_violin(trim = FALSE, fill = "lightblue") +
  stat_summary(fun = median, geom = "point", color = "red", size = 2) +
  theme_bw()

#------------Preparing data for Clustering--------------

clustering_data <- analysis_data[, setdiff(colnames(analysis_data), "time")]
#Removing low information channels
clustering_data <- clustering_data[, 
                                   !grepl("Sn|Xe|Cs|Ba", colnames(clustering_data))
]
#-----------Step 2 Scale data ---------------------------

clustering_scaled <- scale(clustering_data)

#-----------Step 3 Run clustering K-means ---------------

set.seed(42)
k <- 8

kmeans_res <- kmeans(clustering_scaled, centers = k)

analysis_data$cluster <- as.factor(kmeans_res$cluster)

#-------------UMAP (visualization)-----------------------

install.packages("uwot")   # if not installed
library(uwot)

umap_res <- umap(clustering_scaled)
analysis_data$UMAP1 <- umap_res[,1]
analysis_data$UMAP2 <- umap_res[,2]

#-------------Plot clusters-----------------------------

ggplot(analysis_data, aes(x = UMAP1, y = UMAP2, color = cluster)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_bw()

#--------- naming the clusters--------------------------

library(dplyr)

cluster_summary <- analysis_data %>%
  group_by(cluster) %>%
  summarise(across(where(is.numeric), median))

cluster_summary

#------------------------------------------------------

cluster_summary %>%
  select(cluster,
         `170Er_CD3`,
         `144Nd_CD19`,
         `168Er_CD14`,
         `163Dy_CD56_NCAM`,
         `171Yb_CD20`,
         `173Yb_HLA-DR`,
         `165Ho_yH2AX`)
