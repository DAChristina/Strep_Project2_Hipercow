library(magrittr)
library(ape)
library(mlesky)
library(tidyverse)
# 
# if (!require("BiocManager", quietly = TRUE)){
#   install.packages("BiocManager")
#   BiocManager::install("ggtree")
# }


# Functions ####################################################################
# I change the model to fit hipercow requirement by deleting ncpu usage

mod_parboot <-
  function (fit, nrep = 200, ncpu = 1, dd) 
  {
    if (missing(dd)) {
      if (Ntip(fit$tre) <= 500) 
        dd = F
      else dd = T
    }
    if (fit$adapt_time_axis) 
      stop("parboot not supported with adapt_time_axis==TRUE")
    af <- approxfun(fit$time, fit$ne, rule = 2)
    sts <- fit$sampleTimes
    if (is.null(fit$sampleTimes)) {
      sts <- ape::node.depth.edgelength(fit$tre)[1:ape::Ntip(fit$tre)]
      sts = sts - max(sts)
      names(sts) = fit$tre$tip.label
    }
    message("Simulating coalescent trees for parametric bootstrap: ")
    
    # Additional setup parallel cluster for annoying Windows
    cl <- parallel::makeCluster(ncpu)
    on.exit(parallel::stopCluster(cl))
    
    # parallel::clusterExport(cl, c("fit", "af", "sts", "dd", "mlskygrid", "ddSimCoal", "simCoal"))
    parallel::clusterExport(cl,
                            varlist = c("fit", "af", "sts", "dd", "mlskygrid", "ddSimCoal", "simCoal"),
                            envir = environment())
    
    res = parallel::parLapply(cl, 1:nrep, function(irep) {
      if (dd == T) 
        tr = ddSimCoal(sts, alphaFun = af, guessRootTime = min(c(min(sts), 
                                                                 min(fit$time))))
      else tr = simCoal(sts, alphaFun = af)
      f1 <- mlskygrid(tr, sampleTimes = sts, res = fit$res, 
                      tau = fit$tau, tau_tol = fit$tau_tol, ncross = fit$ncross, 
                      quiet = fit$quiet, NeStartTimeBeforePresent = fit$NeStartTimeBeforePresent, 
                      ne0 = median(fit$ne), adapt_time_axis = FALSE, formula = fit$formula, 
                      data = fit$data, ncpu = 1, model = fit$model)
      list(ne = f1$ne, beta = f1$beta, growthrate = f1$growthrate)
    }) #, mc.cores = ncpu)
    nemat <- do.call(cbind, lapply(res, "[[", "ne"))
    min_nemat <- apply(nemat, MARGIN = 1, min)
    max_nemat <- apply(nemat, MARGIN = 1, max)
    lognesd <- apply(log(nemat), MARGIN = 1, sd)
    fit$ne_ci <- cbind(nelb = exp(log(fit$ne) - 1.96 * lognesd), 
                       ne = fit$ne, neub = exp(log(fit$ne) + 1.96 * lognesd),
                       nelr = min_nemat,
                       neur = max_nemat)
    grmat <- do.call(cbind, lapply(res, "[[", "growthrate"))
    grsd <- apply(grmat, MARGIN = 1, sd)
    fit$growthrate_ci <- cbind(grlb = fit$growthrate - 1.96 * 
                                 grsd, gr = fit$growthrate, grub = fit$growthrate + 1.96 * 
                                 grsd)
    if (!is.null(fit$beta)) {
      betamat <- do.call(cbind, lapply(res, "[[", "beta"))
      fit$beta_ci <- cbind(betalb = apply(betamat, MARGIN = 1, 
                                          FUN = function(x) quantile(x, prob = 0.025)), beta = fit$beta, 
                           betaub = apply(betamat, MARGIN = 1, FUN = function(x) quantile(x, 
                                                                                          prob = 0.975)))
    }
    fit
  }

get_random_tree <-
  function(input_tree, proportion = 1.0) {
    random_taxa <-
      unique(sample(input_tree$tip.label,size = round(proportion*length(input_tree$tip.label)), replace = TRUE))
    
    # Get functional random tree
    maxdepth <- NA
    random_tree <- NA
    while (is.na(maxdepth)) {
      random_tree <-
        ape::keep.tip(input_tree,random_taxa)
      maxdepth <-
        max(ape::node.depth.edgelength(random_tree))
    }
    
    return(list("tree" = random_tree, "height" = maxdepth))
  }

fit_model <-
  function(input_tree, model_index = 1, tau = NULL, res = NULL, adapt = TRUE, threads = NULL) {
    
    # Fit model
    mlesky_obj <-
      mlskygrid(
        input_tree,
        tau = tau,
        res = res,
        quiet = FALSE,
        model = model_index,
        ncross = 5,
        tau_tol = 0.001,
        tau_lower = 0.001,
        tau_upper = 10000,
        adapt_time_axis = adapt,
        ncpu = threads
      )
    
    return(mlesky_obj)
  }

make_mlesky_df <- function(obj,input_tree,model_index,adapt,most_recent,strain) {
  
  # Set model names
  model_name <- c("skysigma","skygrid","skygrowth")
  
  mlesky_df <-
    tibble::tibble(
      "Strain" = strain,
      "Year" = most_recent+obj$time,
      "Ne" = obj$ne,
      "Model" = model_name[model_index],
      "Tau" = obj$tau,
      "Resolution" = obj$res,
      "Adapt" = ifelse(adapt,"Adaptive sampling by time","No adaptive sampling"),
      "Origin" = most_recent-max(ape::node.depth.edgelength(input_tree))
    )
  
  return(mlesky_df)
}

fit_mlesky_model <-
  function(input_tree, model_index = 1, adapt = TRUE, most_recent = 2022, strain = "", replicates = 0, proportion = 0.9, ncpu = 4) {
    
    # Fit model
    mlesky_obj <- fit_model(input_tree,
                            model_index = model_index,
                            adapt = adapt,
                            threads = ncpu)
    
    # Process output
    mlesky_df <- make_mlesky_df(mlesky_obj,
                                input_tree,
                                model_index,
                                adapt,
                                most_recent,
                                strain)
    
    # Calculate CIs
    mlesky_df <-
      calculate_CIs(mlesky_df,
                    mlesky_obj,
                    input_tree,
                    strain = strain,
                    model_index = 1,
                    adapt = adapt,
                    most_recent,
                    replicates = replicates,
                    proportion = proportion,
                    ncpu = ncpu,
                    threads = ncpu)
    
    # Return
    return(mlesky_df)
  }

calculate_CIs <-
  function(df,mlesky_obj,input_tree,strain = strain,model_index = 1,adapt = TRUE,most_recent,threads = ncpu,replicates = 0,proportion = 0.9,ncpu = 4) {
    
    if (!adapt) {
      boot_df <- mod_parboot(mlesky_obj,nrep = replicates,ncpu = ncpu, dd = FALSE)
      ci_df <-
        tibble::tibble(
          "Year" = boot_df$time + most_recent,
          "Lower_bound" = boot_df$ne_ci[,1],
          "Upper_bound" = boot_df$ne_ci[,3],
          "Lower_range" = boot_df$ne_ci[,4],
          "Upper_range" = boot_df$ne_ci[,5]
        )
    } else {
      input_tree_height = max(ape::node.depth.edgelength(input_tree))
      rtree_dfs <- list()
      for (i in 1:replicates) {
        rtree <- get_random_tree(input_tree, proportion = proportion)
        output_obj <- fit_model(rtree$tree,
                                tau = NULL,
                                res = round(rtree$height/(input_tree_height/df$Resolution[1])),
                                model_index = model_index,
                                adapt = adapt,
                                threads = ncpu)
        # output_obj <- fit_model(input_tree,
        #                         tau = NULL,
        #                         res = rpois(1,mlesky_obj$res),
        #                         model_index = model_index,
        #                         adapt = adapt,
        #                         threads = ncpu)
        output_df <- make_mlesky_df(output_obj,
                                    input_tree,
                                    model_index = model_index,
                                    adapt,
                                    most_recent,
                                    strain)
        rtree_dfs[[length(rtree_dfs)+1]] <- 
          tibble::tibble(
            "Year" = df$Year[df$Year >= min(output_df$Year) & df$Year <= max(output_df$Year)],
            "Ne" = spline(output_df$Year,output_df$Ne,xout = df$Year[df$Year >= min(output_df$Year) & df$Year <= max(output_df$Year)])$y
          )
      }
      ci_df <-
        dplyr::bind_rows(rtree_dfs) %>%
        dplyr::mutate(Ne = dplyr::if_else(Ne>0,Ne,0)) %>%
        dplyr::group_by(Year) %>%
        dplyr::summarise(
          "Lower_bound" = quantile(Ne,probs = 0.025),
          "Upper_bound" = quantile(Ne,probs = 0.975),
          "Lower_range" = min(Ne),
          "Upper_range" = max(Ne)
        )
    }
    
    # Join to main DF
    df %<>%
      dplyr::left_join(ci_df,
                       by = c("Year"))
    
    return(df)
    
  }


# Load #########################################################################
# args <- commandArgs(trailingOnly = TRUE)
# if (length(args) != 4) {
#   quit(1)
# }

# Slightly edit the code to fit Hipercow upload!
run_mlesky <- function(rep){
  bactdating_obj <- readRDS("outputs/genomics/choosen_n703/method_strictgamma_1e6/mcmc_bacdating.rds")
  # bactdating_obj <- readRDS(args[1])
  
  # model_index <- as.integer(bactdating_obj[3]) # because mcmc_bacdating[1] is from gubbins' & mcmc_bacdating[2] = tree
  model_index <-
    ifelse(bactdating_obj[[3]] == "strictgamma", 1, NA_integer_)
  
  adapt_val <- ifelse(bactdating_obj[4] == "TRUE",TRUE,FALSE)
  # threads <- as.integer(bactdating_obj[5])
  threads <- 1
  
  # Get strain names
  strain_name <- c("GPSC31")
  # ifelse(grepl("GPSC31", bactdating_obj[2]),
  #        "GPSC31",
  #        "GPSC2"
  # )
  
  most_recent_date <- list()
  most_recent_date[["GPSC31"]] <- 2021.989041
  # most_recent_date[["GPSC2"]] <- 2022.221918
  
  message(paste("Strain:",strain_name," Model index: ",model_index," Adaptation: ",adapt_val))
  
  model_output <- fit_mlesky_model(bactdating_obj$tree,
                                   strain = strain_name,
                                   model_index = model_index,
                                   most_recent = most_recent_date[[strain_name]],
                                   adapt = adapt_val,
                                   ncpu = threads,
                                   replicates = rep)
  
  # Write output
  output_fn = paste0("outputs/genomics/", strain_name,"_model_",model_index,"_adaptation_",ifelse(adapt_val,"true","false"),"_mlesky.csv")
  write.csv(model_output,
            file = output_fn,
            quote = FALSE,
            row.names = FALSE
  )
  
}

