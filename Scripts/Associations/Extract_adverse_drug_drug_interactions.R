set.seed(5081)


# Script to analyse the drug-drug interactions (DDIs) reported in DrugBank


# Load libraries
library(unixtools)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


# Read the DrugBank data
DrugBank_data <- readRDS("Databases/DrugBank/parsed_DrugBank_data.rds")


# Read the drug type information
DrugBank_drug_type <- DrugBank_data$drugs$general_information
DrugBank_drug_type <- DrugBank_drug_type[DrugBank_drug_type$type == "small molecule", ] # retain only small molecular drugs


#####


# Extract all possible drug interactions
DrugBank_ddi <- DrugBank_data$drugs$drug_interactions
colnames(DrugBank_ddi)[c(1,2,4)] <- c("Drug1_DrugBank_id", "Drug1_DrugBank_name", "Drug2_DrugBank_id")
DrugBank_ddi <- DrugBank_ddi[DrugBank_ddi$Drug1_DrugBank_id %in% DrugBank_drug_type$drugbank_id & DrugBank_ddi$Drug2_DrugBank_id %in% DrugBank_drug_type$drugbank_id, ] # retain only small molecular drugs
# DrugBank_ddi$Drug2_DrugBank_name <- DrugBank_drug_type$name[match(DrugBank_ddi$Drug2_DrugBank_id, DrugBank_drug_type$drugbank_id)]



# dt <- as.data.table(DrugBank_ddi)
# dt[, description := pbmapply(
#   
#   function(desc, d1_name, d2_name) {
#     
#     # Build pattern and replacement vectors
#     names_vec <- c(d1_name, d2_name)
#     repl_vec <- c("DRUG", "DRUG")
#     
#     # Order by name length descending (to prioritize longer match)
#     order_idx <- order(nchar(names_vec), decreasing = TRUE)
#     names_vec <- names_vec[order_idx]
#     repl_vec <- repl_vec[order_idx]
#     
#     # Use fixed() to avoid regex interpretation
#     for (i in seq_along(names_vec)) {
#       desc <- str_replace_all(desc, fixed(names_vec[i]), repl_vec[i])
#     }
#     return(desc)
#   },
#   description, Drug1_DrugBank_name, Drug2_DrugBank_name,
#   USE.NAMES = FALSE
# )]




#####


# Mark which drug combinations have serious ADR
DrugBank_ddi$ADR_status <- "unmarked"
DrugBank_ddi[grepl("risk or severity of .+toxicity can be increased|risk or severity of liver damage can be increased|risk or severity of adverse effects can be increased |increase the .+toxic activities", DrugBank_ddi$description, ignore.case = TRUE), "ADR_status"] <- "adr_positive"


#####


# Save
saveRDS(DrugBank_ddi, file = "IOFiles/Associations/DrugBank_DDI_processed.rds")


#####


print(warnings())