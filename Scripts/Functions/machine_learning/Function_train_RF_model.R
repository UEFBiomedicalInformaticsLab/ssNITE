# Function to train ML model


func_train_model <- function( ml_data,  
                              ml_target_type = NULL, 
                              ml_sample_folds = NULL,
                              n_estimator = 250, 
                              max_features = "sqrt", 
                              n_cores = 1, 
                              n_jobs = 1 ){
  
  require(reticulate)
  require(parallel)
  
  if(!ml_target_type %in% c("regression", "classification")){ stop("'ml_target_type' should be either 'regression' or 'classification'.") }
  
  if(!all(c("sample_id", "ml_target") %in% colnames(ml_data))){ stop("'ml_data' must contain the columns 'sample_id' and 'ml_target'")}
  
  if( (ml_target_type %in% "classification") & (class(ml_data$ml_target) != "character")){
    stop("'ml_target' column should be a character column and not a factor.")
  }
  
  
  #####
  
  
  # Check for python virtual environment
  env_name <- paste0("Environment/python_env_kudos")
  source("Scripts/Functions/ssNITE/Function_ssNite__load_py.R")
  load_python_env(env_name, create_missing_env = FALSE)
  
  
  #####
  
  
  # Get the number of repeats and folds in the stratified sample list
  n_repeats <- unique(unlist(lapply(ml_sample_folds, length)))
  n_folds <- unique(unlist(lapply(ml_sample_folds, lengths)))
  
  
  # Check ml_data is not empty. If TRUE return NA for all repeat/fold
  check_empty_1 <- ml_data %>% select(!c("sample_id", "ml_target")) %>% as.matrix()
  check_empty_1 <- all(is.na(check_empty_1))
  if(check_empty_1){message("'ml_data' is empty")}
  
  check_empty_2 <- ml_data %>% select(!c("sample_id", "ml_target")) %>% as.matrix()
  check_empty_2 <- apply(check_empty_2, 2, var, na.rm = TRUE)
  check_empty_2 <- all( (is.na(check_empty_2)) | (check_empty_2 == 0) )
  if(check_empty_2){message("'ml_data' contains all zero variance cols or mix of zero variance and NA cols")}
  
  if(check_empty_1 | check_empty_2){
    
    repeat_fold_result <- expand.grid( "n_repeat" = 1:n_repeats, 
                                       "n_fold" = 1:n_folds, 
                                       "accuracy_score" = NA_real_, 
                                       "accuracy_metric" = ifelse(ml_target_type == "classification", "F1", "R2"), 
                                       stringsAsFactors = FALSE )
    
  }else{
    
    # Filter out NA cols
    na_cols <- ml_data %>% select(!c("sample_id", "ml_target")) %>% as.matrix()
    na_cols <- apply(na_cols, 2, function(x)sum(is.na(x)))
    na_cols <- names(which(na_cols == nrow(ml_data)))
    ml_data <- ml_data %>% select(!all_of(na_cols))
    if(length(na_cols)>0){message(paste0("Following columns will be removed from analysis due to zero variance: ", na_cols))}
    rm(na_cols)
    
    # Run model training/testing for all folds and repeats
    repeat_fold_grid <- expand.grid( repeat_count = 1:n_repeats, fold = 1:n_folds )
    repeat_fold_result <- mclapply(X = 1:nrow(repeat_fold_grid), 
                                   mc.preschedule = TRUE,
                                   mc.cores = n_cores, 
                                   FUN = function(j){
                                     
                                     repeat_count <- repeat_fold_grid$repeat_count[j]
                                     fold_count <- repeat_fold_grid$fold[j]
                                     
                                     # Get the list of samples
                                     train_samples <-  ml_sample_folds$train_samples[[repeat_count]][[fold_count]] 
                                     test_samples <-   ml_sample_folds$test_samples[[repeat_count]][[fold_count]] 
                                     
                                     # Get the target values for the training and test data
                                     y_train <- ml_data %>% filter(sample_id %in% train_samples) %>% select(c("sample_id", "ml_target"))
                                     y_test <- ml_data %>% filter(sample_id %in% test_samples) %>% select(c("sample_id", "ml_target"))
                                     
                                     # Get the input values of the training and test data
                                     x_train <- ml_data %>% filter(sample_id %in% train_samples) %>% select(!c("ml_target")) %>% column_to_rownames("sample_id")
                                     x_test <- ml_data %>% filter(sample_id %in% test_samples) %>% select(!c("ml_target")) %>% column_to_rownames("sample_id")
                                     
                                     if(any(is.na(x_train))){stop("Data contains NA")}
                                     if(any(is.na(x_test))){stop("Data contains NA")}
                                     
                                     
                                     # Train the random forest model
                                     
                                     if(ml_target_type == "classification"){
                                       rf <- sk$ensemble$RandomForestClassifier(n_estimators = as.integer(n_estimator), 
                                                                                max_features = max_features, 
                                                                                n_jobs = as.integer(n_jobs), 
                                                                                random_state = as.integer(5081))
                                     }else if(ml_target_type == "regression"){
                                       rf <- sk$ensemble$RandomForestRegressor(n_estimators = as.integer(n_estimator), 
                                                                               max_features = max_features, 
                                                                               n_jobs = as.integer(n_jobs), 
                                                                               random_state = as.integer(5081))
                                     }
                                     
                                     rf <- rf$fit( r_to_py(x_train), 
                                                   r_to_py(y_train %>% pull("ml_target"))  )
                                     
                                     # Generate predictions
                                     y_test$ml_predictions <- rf$predict(r_to_py( x_test )) %>% as.vector()
                                     
                                     # Calculate accuracy
                                     if(ml_target_type == "classification"){
                                       acc_score <- sk$metrics$f1_score( y_true = r_to_py(y_test %>% pull(ml_target)), 
                                                                         y_pred = r_to_py(y_test %>% pull(ml_predictions)), 
                                                                         average = "weighted" )
                                       acc_metric <- "F1"
                                     }else if(ml_target_type == "regression"){
                                       acc_score <- sk$metrics$r2_score( y_true = r_to_py(y_test %>% pull(ml_target)), 
                                                                         y_pred = r_to_py(y_test %>% pull(ml_predictions)) )
                                       acc_metric <- "R2"
                                     }
                                     
                                     
                                     # Calculate feature importance
                                     feature_imp <- as.data.frame(t(setNames( object = py_to_r(rf$feature_importances_), 
                                                                              nm = paste0("imp_", py_to_r(rf$feature_names_in_)) )))
                                     
                                     
                                     # prepare result for export
                                     tmp1 <- data.frame("n_repeat" = repeat_count,
                                                        "n_fold" = fold_count, 
                                                        "accuracy_score" = py_to_r(acc_score),
                                                        "accuracy_metric" = acc_metric,
                                                        feature_imp, 
                                                        check.names = FALSE)
                                     
                                     return(tmp1)
                                     
                                   })
    
    repeat_fold_result <- bind_rows(repeat_fold_result)
    
  }
  
  return(repeat_fold_result)
  
}
