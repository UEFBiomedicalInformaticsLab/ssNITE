# Function to load python environment via reticulate

load_python_env <- function( env_name, 
                             create_missing_env = FALSE ){
  
  require(reticulate)
  
  # Prevent pyc race conditions (job-local cache)
  Sys.setenv(PYTHONPYCACHEPREFIX = file.path(tempdir(), "pycache"))
  
  # Check if the virtual environment is installed.
  # If not, install if requested.
  if(!virtualenv_exists(env_name)){
    
    print(paste0(Sys.time(), " | Reticulate environment missing."))
    
    if(create_missing_env){
      
      print(paste0(Sys.time(), " | Preparing new reticulate environment!"))
      system("Rscript Environment/Reticulate_env_py.R")
      
    }else{
      
      stop("Reticulate environment needed to continue!")
      
    }
    
  }
  
  
  # Check if the environment is already loaded
  if( exists("py_fn", envir = .GlobalEnv) &&
      exists("pd", envir = .GlobalEnv) &&
      exists("np", envir = .GlobalEnv) &&
      exists("nz", envir = .GlobalEnv) ){
    
    print(paste0(Sys.time(), " | Python environment already loaded."))
    
  }else{
    
    print(paste0(Sys.time(), " | Using reticulate environment from: ", env_name))
    use_virtualenv(env_name, required = TRUE)
    pd <<- import("pandas", convert = FALSE)
    np <<- import("numpy", convert = FALSE)
    nz <<- import("netZooPy.bonobo", convert = FALSE)
    py_fn <<- import_builtins(convert = FALSE)
    
  }
  
}