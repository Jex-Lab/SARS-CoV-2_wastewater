# This script uses the mapping_table_MASTER.xlsx file to create:
# 1 - the mapping table (tab-separated format) used by the main WGS script
# 2 - the mapping table (Excel format) provided to DoH as a reference 

library(readxl)
library(writexl)

# import excel file
ay <- read_xlsx("mapping_table_MASTER.xlsx", sheet='AY', col_names=FALSE)
ay <- ay[,c(1,3)]
colnames(ay) <- c("Lineage","Mapping")

ba1 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA1', col_names=FALSE)
ba1 <- ba1[,c(1,3)]
colnames(ba1) <- c("Lineage","Mapping")

ba2 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA2', col_names=FALSE)
ba2 <- ba2[,c(1,3)]
colnames(ba2) <- c("Lineage","Mapping")

ba3 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA3', col_names=FALSE)
ba3 <- ba3[,c(1,3)]
colnames(ba3) <- c("Lineage","Mapping")

ba4 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA4', col_names=FALSE)
ba4 <- ba4[,c(1,3)]
colnames(ba4) <- c("Lineage","Mapping")

ba5 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA5', col_names=FALSE)
ba5 <- ba5[,c(1,3)]
colnames(ba5) <- c("Lineage","Mapping")

ba275x <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BA275x', col_names=FALSE)
ba275x <- ba275x[,c(1,3)]
colnames(ba275x) <- c("Lineage","Mapping")

bn <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BN', col_names=FALSE)
bn <- bn[,c(1,3)]
colnames(bn) <- c("Lineage","Mapping")

bq1x <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BQ1x', col_names=FALSE)
bq1x <- bq1x[,c(1,3)]
colnames(bq1x) <- c("Lineage","Mapping")

br2 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='BR2', col_names=FALSE)
br2 <- br2[,c(1,3)]
colnames(br2) <- c("Lineage","Mapping")

ch11 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='CH11', col_names=FALSE)
ch11 <- ch11[,c(1,3)]
colnames(ch11) <- c("Lineage","Mapping")

xbb <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBB', col_names=FALSE)
xbb <- xbb[,c(1,3)]
colnames(xbb) <- c("Lineage","Mapping")

xbb15 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBB15', col_names=FALSE)
xbb15 <- xbb15[,c(1,3)]
colnames(xbb15) <- c("Lineage","Mapping")

xbb191 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBB191', col_names=FALSE)
xbb191 <- xbb191[,c(1,3)]
colnames(xbb191) <- c("Lineage","Mapping")

xbb116 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBB116', col_names=FALSE)
xbb116 <- xbb116[,c(1,3)]
colnames(xbb116) <- c("Lineage","Mapping")

xbc <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBC', col_names=FALSE)
xbc <- xbc[,c(1,3)]
colnames(xbc) <- c("Lineage","Mapping")

xbf <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBF', col_names=FALSE)
xbf <- xbf[,c(1,3)]
colnames(xbf) <- c("Lineage","Mapping")

recomb <- read_xlsx("mapping_table_MASTER.xlsx", sheet='Recombinant', col_names=FALSE)
recomb <- recomb[,c(1,3)]
colnames(recomb) <- c("Lineage","Mapping")

xbb192 <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XBB192', col_names=FALSE)
xbb192 <- xbb192[,c(1,3)]
colnames(xbb192) <- c("Lineage","Mapping")

xcc <- read_xlsx("mapping_table_MASTER.xlsx", sheet='XCC', col_names=FALSE)
xcc <- xcc[,c(1,3)]
colnames(xcc) <- c("Lineage","Mapping")

all_maps <- rbind(ay, ba1, ba2, ba3, ba4, ba5, ba275x, bn, bq1x, br2, ch11, xbb, xbb15, xbb191, xbb116, xbc, xbf, recomb, xbb192, xcc)

write.table(all_maps, "mapping_table.txt", row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")

write_xlsx(all_maps, "mapping_table.xlsx", col_names=TRUE, format_headers=FALSE)

