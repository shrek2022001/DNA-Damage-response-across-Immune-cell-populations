#-------Build a marker annotation dataframe manually---------------------

t_markers <- c("CD3", "CD4", "CD8a", "TCRgd", "CD45RO", "CD45RA",
               "CCR7", "CD27", "CD28", "CD25", "CD127", "TIM3",
               "TIGIT", "CXCR3", "CCR4", "CCR6", "CXCR5", "CD161")

b_markers <- c("CD19", "CD20", "IgD", "CD38", "HLA-DR")

nk_markers <- c("CD56", "CD16", "NKG2A", "CD57", "CD294")

myeloid_markers <- c("CD14", "CD11c", "CD123", "CD66b", "HLA-DR")

exclude <- c(
  "Time","Event_length","Center","Width","Residual","Offset","Amplitude",
  "Barcode106","Barcode112","Barcode114","Barcode116","Barcode195","Barcode198",
  "Bead","Bead.1","190BCKG","194Pt","208Pb",
  "DNA1","DNA2","Puro","IdU","120Sn","131Xe","133Cs","138Ba",
  "Live_Dead","file_id","timepoint"
)

bio_markers <- setdiff(colnames(data), exclude)
bio_markers
t_markers <- c()
b_markers <- c()
nk_markers <- c()
myeloid_markers <- c()

grep("CD3|CD4|CD8|TCR|CCR|CD27|CD28|CD45|CD161", bio_markers, value = TRUE)

#Finding B related markers
grep("CD19|CD20|IgD", bio_markers, value = TRUE)

#Finding NK-related markers
grep("CD56|CD16|NKG2A|CD57|CD294", bio_markers, value = TRUE)

#Finding Myeloid related markers
grep("CD14|CD11c|CD123|CD66b|HLA", bio_markers, value = TRUE)

#Use these groups:
t_markers <- c("CD3", "CD4", "CD8a", "TCRgd", "CCR7", "CCR6", "CCR4",
               "CD27", "CD28", "CD45RA", "CD45RO")

b_markers <- c("CD19", "CD20", "IgD")

nk_markers <- c("CD56", "CD16", "NKG2A", "CD57", "CD294")

myeloid_markers <- c("CD14", "CD11c", "CD123", "CD66b", "HLA-DR")

#Safety check
t_markers <- intersect(t_markers, colnames(data))
b_markers <- intersect(b_markers, colnames(data))
nk_markers <- intersect(nk_markers, colnames(data))
myeloid_markers <- intersect(myeloid_markers, colnames(data))

t_markers
b_markers
nk_markers
myeloid_markers

data$T_score <- rowMeans(data[, t_markers], na.rm = TRUE)
data$B_score <- rowMeans(data[, b_markers], na.rm = TRUE)
data$NK_score <- rowMeans(data[, nk_markers], na.rm = TRUE)
data$Myeloid_score <- rowMeans(data[, myeloid_markers], na.rm = TRUE)

#Assigning each cell to the strongest group
score_mat <- data[, c("T_score", "B_score", "NK_score", "Myeloid_score")]

data$functional_group <- c("T_cell", "B_cell", "NK_cell", "Myeloid")[
  max.col(score_mat, ties.method = "first")
]

#See how many cells landed in each group
table(data$functional_group)
prop.table(table(data$functional_group))

#Is the data biologically sensible?
aggregate(data[, c("CD3", "CD19", "CD56", "CD14", "HLA-DR")],
          by = list(group = data$functional_group),
          median)
data$max_score <- apply(score_mat, 1, max)
hist(data$max_score, breaks = 50)
threshold <- 2.5
data$functional_group_refined <- ifelse(
  data$max_score < threshold,
  "Unassigned",
  data$functional_group
)
table(data$functional_group_refined)
prop.table(table(data$functional_group_refined))
aggregate(data[, c("CD3","CD19","CD56","CD14","HLA-DR")],
          by = list(group = data$functional_group_refined),
          median)

data$functional_group_final <- data$functional_group_refined

data$functional_group_final[
  data$functional_group_refined == "Unassigned" & data$CD3 > 2.5
] <- "T_cell"
table(data$functional_group_final)
prop.table(table(data$functional_group_final))

gamma_summary <- aggregate(
  data$yH2AX,
  by = list(functional_group = data$functional_group_final,
            timepoint = data$timepoint),
  median
)

colnames(gamma_summary)[3] <- "median_gamma"
gamma_summary

library(ggplot2)

ggplot(gamma_summary,
       aes(x = factor(timepoint), y = median_gamma,
           color = functional_group, group = functional_group)) +
  geom_line() +
  geom_point() +
  theme_minimal()

