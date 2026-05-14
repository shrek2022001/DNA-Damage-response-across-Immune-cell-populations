#------PHASE1 SETUP DATA AND MARKERS-------------------
setwd("C:/Users/gaura/OneDrive/Documents/shreya/hkfz_data")
files <- c("c01_DDR_Processed.fcs",
           "c02_DDR_Processed.fcs",
           "c03_DDR_Processed.fcs",
           "c04_DDR_Processed.fcs")

fs <- read.flowSet(files, transformation = FALSE)

#Adding metadata(timepoints)
timepoints <- c(0, 30, 60, 360)

class(fs)
length(fs)
sampleNames(fs)

expr_list <- lapply(seq_along(fs), function(i) {
  exprs(fs[[i]])
})
length(expr_list)
dim(expr_list[[1]])
dim(expr_list[[2]])

timepoints <- c(0, 30, 60, 360)
file_ids <- sampleNames(fs)

expr_list <- mapply(function(mat, tp, fid) {
  mat <- as.data.frame(mat)
  mat$timepoint <- tp
  mat$file_id <- fid
  mat
}, expr_list, timepoints, file_ids, SIMPLIFY = FALSE)
data <- do.call(rbind, expr_list)
dim(data)
head(data[, 1:10])
tail(colnames(data), 5)
table(data$timepoint)
table(data$file_id)

summary(data[, 1:5])

numeric_cols <- sapply(data, is.numeric)
numeric_cols["timepoint"] <- FALSE

data[, numeric_cols] <- asinh(as.matrix(data[, numeric_cols]) / 5)

colnames(data)

ff1 <- fs[[1]]

pdat <- pData(parameters(ff1))
pdat[, c("name", "desc")]

#Creating a clean marker table
marker_table <- pdat[, c("name", "desc")]
marker_table

#Removed technical and non-biological channels
exclude_channels <- c(
  "Time", "Event_length", "Center", "Width", "Residual", "Offset", "Amplitude",
  "Ce140Di", "Ce142Di",          # beads
  "Ir191Di", "Ir193Di",          # DNA
  "Pt194Di", "Pt195Di", "Pt198Di", # barcodes
  "Cd106Di", "Cd112Di", "Cd114Di", "Cd116Di", # barcodes
  "BCKG190Di", "Pb208Di"
)

marker_table_clean <- marker_table[!(marker_table$name %in% exclude_channels), ]
marker_table_clean
idx <- is.na(marker_table$desc)

marker_table$desc[idx] <- marker_table$name[idx]
marker_table$desc <- make.unique(marker_table$desc)
colnames(data) <- marker_table$desc[match(colnames(data), marker_table$name)]
colnames(data)
colnames(data) <- gsub("^[0-9]+[A-Za-z]+_", "", colnames(data))
original_cols <- colnames(data)

matched_desc <- marker_table$desc[match(original_cols, marker_table$name)]

new_cols <- ifelse(is.na(matched_desc), original_cols, matched_desc)
new_cols <- make.unique(new_cols)

colnames(data) <- new_cols
colnames(data)
colnames(data)[64] <- "timepoint"
colnames(data)[65] <- "file_id"

colnames(data)[colnames(data) == "CD196_CCR6"] <- "CCR6"
colnames(data)[colnames(data) == "CD123_IL-3R"] <- "CD123"
colnames(data)[colnames(data) == "CD194_CCR4"] <- "CCR4"
colnames(data)[colnames(data) == "CD25_IL-2Ra"] <- "CD25"
colnames(data)[colnames(data) == "CD183_CXCR3"] <- "CXCR3"
colnames(data)[colnames(data) == "CD185_CXCR5"] <- "CXCR5"
colnames(data)[colnames(data) == "CD56_NCAM"] <- "CD56"
colnames(data)[colnames(data) == "CD197_CCR7"] <- "CCR7"
colnames(data)[colnames(data) == "CD127_IL-7Ra"] <- "CD127"
colnames(data)[colnames(data) == "TGIT"] <- "TIGIT"

