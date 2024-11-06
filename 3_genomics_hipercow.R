
# The pmcmc run by using Hipercow! #############################################
# https://mrc-ide.github.io/hipercow/articles/windows.html
library(hipercow)
library(parallel)
library(ids)
# sudo mount -a

# See filesystems and paths, CHANGE wd to those in /etc/fstab
# DO NOT CHANGE THE ARRANGEMENT OF THESE COMMANDS!!!

hipercow_init(driver = "windows")
hipercow_configure("windows", r_version = "4.4.0")
windows_authenticate() # authenticate by using DIDE account
windows_check()
hipercow_configuration() # for troubleshooting
# hipercow_hello() # test job

# Check the installed packages again by using hipercow_configuration()
# hipercow_configuration()

# If automatic install failed (don't know why), use pkgdepends.txt!
# install.packages("pkgdepends")
# hipercow_provision()

# https://mrc-ide.github.io/hipercow/reference/hipercow_resources.html
resources <- hipercow::hipercow_resources(cores = 20,
                                          max_runtime = "60d",
                                          memory_per_node = "64G",
)


# BactDating #######################################################################
hipercow_environment_create(sources = "genomics/1_treedater_bacdating_GPSC31.R")
hipercow_provision()
run_bactdating <- task_create_expr(run_bactdating(1e6), # try 1e8
                             resources = resources
)

# Something related to test the submitted job
task_status(run_bactdating)
task_result(run_bactdating)
task_log_show(run_bactdating)
task_info(run_bactdating)
task_info(run_bactdating)$times

# mleSky #######################################################################
hipercow_environment_create(sources = "genomics/2_mlesky.R")
hipercow_provision()
run_pop <- task_create_expr(run_mlesky(1e6), # try 1e8
                            resources = resources
)

# Something related to test the submitted job
task_status(run_pop)
task_result(run_pop)
task_log_show(run_pop)
task_info(run_pop)
task_info(run_pop)$times


# FINALLY I SURVIVED & SUCCEEEDED PARALLELLISATION!!!!!!!!!!!!!!!!!!!!!!!!!#####
# https://mrc-ide.github.io/hipercow/reference/hipercow_parallel.html

parallel <- hipercow::hipercow_parallel(method = "parallel",
                                        cores_per_process = 8)
cl <- parallel::makeCluster(getOption("cl.cores", 20))
id_parallel <- task_create_expr(parallel::clusterApply(NULL, 1:8, function(x) run_mlesky(1e6)),
                                # c(Sys.getpid(), hipercow_parallel_get_cores()),
                                # parallel = hipercow_parallel("future", cores_per_process = 20),
                                parallel = parallel,
                                resources = resources)

task_status(id_parallel)
task_result(id_parallel)
task_log_show(id_parallel)
task_info(id_parallel)
task_info(id_parallel)$times

# Keith J. Fraser:
# Sometimes adding more cores doesn't improve speed
# I find my parallelized code where operations on a large dataset 
# are split between multiple cores, with part of the dataset being sent to each core, 
# has a rough optimum number of cores to use. 
# I assume this is because with too many cores, 
# bringing the dataset back together becomes the rate determining step 
# (and/or all the threads get held back by whichever one is slowest).
