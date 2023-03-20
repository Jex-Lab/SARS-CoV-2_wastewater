# This code transforms AmpSeq data exported from the FileMaker Pro
# database (via the 'VOC' button) into a format that is compatible
# with the DoH's Power-Bi system. The final Excel file can be
# directly uploaded to the DoH's SharePoint.

library(readxl)
library(reshape2)
library(writexl)

# import excel file
report_file <- read_xlsx("export.xlsx")

# reformat date
report_file$`Batch date` <- format(report_file$`Batch date`, format = "%d/%m/%Y")

# convert to long data format
report_file_transform <- melt(report_file, id=(colnames(report_file[grep("Mapping_",colnames(report_file),invert=T)])))

# remove "Mapping_" text
report_file_transform$variable <- gsub(x=report_file_transform$variable, pattern="Mapping_",replacement="")

# set column names
colnames(report_file_transform) <- c("RNA_ID",
                                     "AmpSeq Run",
                                     "DHHS SampleNo",
                                     "Batch Date",
                                     "RepeatNumber",
                                     "Filtered Reads",
                                     "Non-classified SNPs",
                                     "Non-classified SNPs comment",
                                     "Recommendations",
                                     "Technical notes",
                                     "Variant",
                                     "Variant reads")

# add columns
report_file_transform$Purpose <- c("VOC - Routine")
report_file_transform$Amplicon <- c("RBM")
report_file_transform$CodonRange <- c("440-503")
report_file_transform$VariantClassificationSNPs <- c("N440K, K444R/M/T, V445P/A, G446S/D/N, N450D, L452R/Q/M, N460K, T470N/A/I, G476C/S, S477N, T478K/R/E, E484A/R/V/T, G485D, F486P/S/I/V/L, F490V/S/L/P, Q493R, S494P, G496S, Q498R, N501Y")

# reorder columnns
report_file_transform <- report_file_transform[, c("Purpose",
                                                   "RNA_ID",
                                                   "AmpSeq Run",
                                                   "DHHS SampleNo",
                                                   "Batch Date",
                                                   "RepeatNumber",
                                                   "Amplicon",
                                                   "CodonRange",
                                                   "VariantClassificationSNPs",
                                                   "Filtered Reads",
                                                   "Variant",
                                                   "Variant reads",
                                                   "Non-classified SNPs",
                                                   "Non-classified SNPs comment",
                                                   "Recommendations",
                                                   "Technical notes")]

# export as excel file
write_xlsx(list(VOC_Results_v2 = report_file_transform), "2023####_###_VOC_Resultsv2.xlsx")
