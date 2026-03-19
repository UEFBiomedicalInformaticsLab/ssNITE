# Script to create python library


# Load library
library(reticulate)


# Create Python virtual environment

env_name <- paste0("Environment/python_env_", Sys.info()[["nodename"]])
py_version <- "3.11"

if(!dir.exists(env_name)){dir.create(env_name, recursive = TRUE)}

install_python(version = py_version, force = TRUE)

virtualenv_create(envname = env_name,
                  python = py_version,
                  force = TRUE)

print("---------------------------------------------------")

pkg_to_install <- c("ipython", "gnureadline", "zstandard")

virtualenv_install(envname = env_name,
                   packages = pkg_to_install,
                   pip_options = c("--no-cache-dir"))


print("---------------------------------------------------")

pkg_to_install <- c("igraph",
                    "matplotlib",
                    "networkx", "netZooPy", "numpy",
                    "pandas", "pyarrow", "pyreadr",
                    "rdata", "rpy2",
                    "scikit-learn",
                    "yacs")

virtualenv_install(envname = env_name,
                   packages = pkg_to_install,
                   pip_options = c("--no-cache-dir"))
