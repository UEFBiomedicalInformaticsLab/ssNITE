set.seed(5081)


# Script to split the TCGA-BRCA samples for ML analysis

# Load libraries
library(unixtools)
library(tidyverse)


# Set temporary directory
if(!dir.exists("tmp_dir/")){dir.create("tmp_dir/", recursive = TRUE)}
set.tempdir("tmp_dir/")


#####


print(paste0(Sys.time(), " | Reading the sample information"))

TCGA_data <- readRDS("IOFiles/GDC_data/vstCounts__DA_TCGA__DT_BreastCancer.rds")
sample_info <- TCGA_data$sample_info %>% filter(condition != "SolidTissueNormal_NA") %>% droplevels()
sample_info <- sample_info %>%
  select(c("barcode", "patient", "race", "paper_pathologic_stage", "paper_BRCA_Subtype_PAM50", "paper_age_at_initial_pathologic_diagnosis"))  %>%
  rename(c("sample_id" = "barcode",
           "Stage" = "paper_pathologic_stage",
           "Subtype" = "paper_BRCA_Subtype_PAM50",
           "Age at diagnosis" = "paper_age_at_initial_pathologic_diagnosis",
           "Race" = "race"))
row.names(sample_info) <- NULL

sample_info[which(is.na(sample_info$Race)), "Race"] <- "unknown"
sample_info[which(sample_info$Race == "not reported"), "Race"] <- "unknown"
remove_race <- names(which(table(sample_info$Race) < 10))
sample_info[which(sample_info$Race %in% remove_race), "Race"] <- "unknown"
rm(list = c("TCGA_data", "remove_race"))
sample_info[which(sample_info$Stage %in% "NA"), "Stage"] <- "unknown"


#####


# Add the survival information for the patients

print(paste0(Sys.time(), " | Adding the survival data for the samples"))

if(!dir.exists("Databases/Genomic_Data_Commons/")){dir.create("Databases/Genomic_Data_Commons/", recursive = TRUE)}
if(!file.exists("Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx")){
  download.file(url = "https://api.gdc.cancer.gov/data/1b5f413e-a8d1-4d10-92eb-7c4ae739ed81",
                destfile = "Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx",
                method = "wget")
}

TCGA_CDR <- readxl::read_xlsx("Databases/Genomic_Data_Commons/TCGA-CDR-SupplementalTableS1.xlsx", sheet = "TCGA-CDR", col_names = TRUE, guess_max = 10000)

TCGA_CDR <- TCGA_CDR %>%
  select(c("bcr_patient_barcode",
           "OS", "OS.time", "DSS", "DSS.time",
           "DFI", "DFI.time", "PFI", "PFI.time"))

sample_info <- sample_info %>%
  left_join(TCGA_CDR,
            by = c("patient" = "bcr_patient_barcode")) %>%
  mutate(across(c("OS", "DSS", "DFI", "PFI"),
                ~ ifelse(. == 1, "yes", "no"))) %>%
  mutate(across(c("OS", "DSS", "DFI", "PFI"),
                ~ ifelse(is.na(.), "unknown", .))) %>%
  mutate(across(c("OS.time", "DSS.time", "DFI.time", "PFI.time"),
                ~ ifelse(is.na(.), -1, .))) %>%
  rename(c("OS time" = "OS.time",
           "DSS time" = "DSS.time",
           "DFI time" = "DFI.time",
           "PFI time" = "PFI.time"))
rm(TCGA_CDR)


#####


# Add the information about histological marker status

print(paste0(Sys.time(), " | Adding the histological marker status for the samples"))

clinicalInfo_patient <- readRDS("IOFiles/GDC_data/TCGABRCA_data.rds")
clinicalInfo_patient <- clinicalInfo_patient$clinicalInfo_patient
clinicalInfo_patient <- clinicalInfo_patient %>%
  select(c("bcr_patient_barcode",
           "breast_carcinoma_estrogen_receptor_status",
           "breast_carcinoma_progesterone_receptor_status",
           "lab_proc_her2_neu_immunohistochemistry_receptor_status")) %>%
  rename(c("ER status" = "breast_carcinoma_estrogen_receptor_status",
           "PR status" = "breast_carcinoma_progesterone_receptor_status",
           "Her2 status" = "lab_proc_her2_neu_immunohistochemistry_receptor_status"))

sample_info <- sample_info %>%
  left_join(clinicalInfo_patient,
            by = c("patient" = "bcr_patient_barcode")) %>%
  mutate(across(c("ER status", "PR status", "Her2 status"),
                ~ ifelse(. %in% c("", "Indeterminate", "Equivocal"), "unknown", tolower(.)))) %>%
  mutate(across(c("ER status", "PR status", "Her2 status"),
                ~ ifelse(is.na(.), "unknown", .)))
rm(clinicalInfo_patient)


#####


# Add the tumour sample purity information

print(paste0(Sys.time(), " | Adding the tumour purity information for the samples"))

purity_data <- TCGAbiolinks::Tumor.purity %>%
  filter(`Cancer.type` %in% "BRCA") %>%
  mutate(across(c(ESTIMATE, ABSOLUTE, LUMP, IHC, CPE), function(x) {
    x <- as.character(x)
    x <- str_replace(x, ",", ".")
    as.numeric(x)
  }))

sample_info <- sample_info %>%
  mutate("tmp_id" = substr(sample_id, 1, 16)) %>%
  left_join(purity_data %>% select(!c("Cancer.type")),
            by = c("tmp_id" = "Sample.ID")) %>%
  select(!c("tmp_id")) %>%
  mutate(across(c("ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE"),
                ~ ifelse((is.nan(.) | is.na(.)), -1, .)))

rm(purity_data)

sample_info <- sample_info %>% select(!c("patient"))
# sample_info[sample_info == "NA"] <- NA
# sample_info[sample_info == "NaN"] <- NA


#####


classification_cols <- c("Race", "Stage", "Subtype", "ER status", "PR status", "Her2 status")
regression_cols <- c("Age at diagnosis", "OS time", "DSS time", "DFI time", "PFI time", "ESTIMATE", "ABSOLUTE", "LUMP", "IHC", "CPE")


#####


# Check for python virtual environment
env_name <- paste0("Environment/python_env_kudos")
source("Scripts/Functions/ssNITE/Function_ssNite__load_py.R")
load_python_env(env_name, create_missing_env = FALSE)


#####


# Generate the train and test data
print(paste0(Sys.time(), " | Stratifying samples data for all repeats+folds"))


n_repeats <- 10
n_folds <- 5

splitted_data <- list()

for(col_sel in classification_cols){
  
  train_data <- vector(mode = "list", length = n_repeats)
  test_data <-  vector(mode = "list", length = n_repeats)
  
  for(repeat_count in 1:n_repeats){
    
    # Split the data into folds
    folds_list <- sk$model_selection$StratifiedKFold(n_splits = as.integer(n_folds),
                                                     shuffle = TRUE,
                                                     random_state = as.integer(repeat_count))
    folds_list <- folds_list$split( X = r_to_py(sample_info$sample_id),
                                    y =  r_to_py(sample_info[[col_sel]]) )
    folds_list <- iterate(folds_list)
    folds_list <- lapply(folds_list, py_to_r)
    
    for(fold in 1:length(folds_list)){
      
      # Get the training and test data
      train_data[[repeat_count]][[fold]] <- sample_info$sample_id[folds_list[[fold]][[1]] + 1] # Using +1 as python index starts at 0
      test_data[[repeat_count]][[fold]]  <- sample_info$sample_id[folds_list[[fold]][[2]] + 1]
      
    }
  }
  splitted_data[[col_sel]] <- list("train_samples" = train_data, "test_samples" = test_data)
  
}



#####


print(paste0(Sys.time(), " | Saving sample list"))

tmp1 <- list("sample_info" = sample_info, 
             "classification_cols" = classification_cols, 
             "regression_cols" = regression_cols, 
             "stratifiedSamples" = splitted_data, 
             "n_repeats" = n_repeats, 
             "n_folds" = n_folds)

if(!dir.exists("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/")){dir.create("IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/", recursive = TRUE)}
saveRDS(tmp1, file = "IOFiles/ML_SSNtopol_vs_sampleInfo/TCGABRCA/ML_stratifiedSamples__DA_TCGA_DT_BreastCancer.rds")


#####


print(warnings())