#' Abundance Removal Model
#' A wrapper function run JAGS-based removal models for stream data.
#' Takes habitat, transect, and count data, processes it, runs the models for each species, and saves the output.
#'
#' @param hab Path to the habitat data file (CSV format).
#' @param tran Path to the transect habitat data file (CSV format)
#' @param count Path to the count data file (CSV format).
#' @param output_dir Directory to save the outputs.
#' @return A list containing count-habitat model summaries and plots.
#' @export
abundance_removal_model <- function(count, output_dir){
  library(dplyr)
  library(tidyr)
  library(jagsUI)
  library(ggplot2)
  library(coda)

  # Step 1: Load data
  count <- count

  # Format count data
  unique_passes <- unique(count$DepletionPass)
  unique_years <- unique(count$Year)
  unique_subunits <- unique(count$Subunit)

  yGIPU <- yPOSO <- yAGCH <- yCYFO <- array(NA, dim = c(length(unique_subunits), length(unique_passes), length(unique_years)))

  for(i in seq_along(unique_subunits)){
    for(j in seq_along(unique_passes)){
      for(t in seq_along(unique_years)){
        subunit <- unique_subunits[i]
        pass <- unique_passes[j]
        year <- unique_years[t]

        data_row <- count %>%
          filter(Subunit == subunit, DepletionPass == pass, Year == year)

        if(nrow(data_row)>0){
          yGIPU[i,j,t] <- data_row$GIPU
          yPOSO[i,j,t] <- data_row$POSO
          yAGCH[i,j,t] <- data_row$AGCH
          yCYFO[i,j,t] <- data_row$CYFO
        }
      }
    }
  }

  # Jags Model
  # POSO Models
  ################################################################################
  modelFilename = "POSO.Removal.Model.txt"
  cat("
  model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions

  # Observation model: Each variable is conditionally distributed around its mean

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega + rhat[i] * (t - 1)
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  for(i in 1:nsites){
    rhat[i] ~ dunif(-1, 1)
  }

  for(i in 1:nsites){
    for(t in 1:nyears){
      latent.g[i, t] ~ dunif(1,npasses)
      g[i, t] <- round(latent.g[i, t])
      p[i, t] ~ dbeta(1, 1)
    }
  }

  for(i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:npasses){
        q[i, j, t] <- 1 - pow((1-p[i, t]), 1/g[i, t])
        #q[i, j, t] ~ dbeta(1, 1)
      }
    }
  }

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      p1_Y[i, t] ~ dbin(q[i, 1, t], N[i, t])
      p2_Y[i, t] ~ dbin(q[i, 2, t], (N[i, t] - p1_Y[i, t]))
      p3_Y[i, t] ~ dbin(q[i, 3, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t]))
      p4_Y[i, t] ~ dbin(q[i, 4, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t]))
      p5_Y[i, t] ~ dbin(q[i, 5, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t] - p4_Y[i, t]))
    }
  }
}
", fill = TRUE,file=modelFilename)


  # Parameters to save
  jags_params <- c(
    "omega",
    "rhat",
    "N"
  )

  # Initialize
  nsites <- dim(yPOSO)[1]
  npasses <- dim(yPOSO)[2]
  nyears <- dim(yPOSO)[3]

  N.inits <- array(NA, dim = c(nsites, nyears))
  p1_Y <- yPOSO[1:nsites,1,1:nyears]
  p2_Y <- yPOSO[1:nsites,2,1:nyears]
  p3_Y <- yPOSO[1:nsites,3,1:nyears]
  p4_Y <- yPOSO[1:nsites,4,1:nyears]
  p5_Y <- yPOSO[1:nsites,5,1:nyears]

  for (i in 1:nsites) {
    for (t in 1:nyears) {
      # Try summing the non-NA values across datasets
      summed_value <- sum(p1_Y[i, t], p2_Y[i, t], p3_Y[i, t],
                          p4_Y[i, t], p5_Y[i, t], na.rm = TRUE)

      # If all values were NA, the sum would be zero (or NA), so check this
      if (is.na(summed_value) || summed_value == 0) {
        # If the summed value is NA or zero, assign a Poisson value
        N.inits[i, t] <- round(rpois(1, lambda = 1))
      } else {
        # Otherwise, use the summed value multiplied by 10
        N.inits[i, t] <- round(summed_value)
      }
    }
  }

  # Bundle data
  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y
  )

  # Initial values
  init_vals <- function() list(omega = rnorm(1, 0, 0.001),
                               N = N.inits)

  # MCMC settings
  ni <- 300000
  nc <- 3
  nb <- 20000
  nt <- 1

  set.seed(1234)

  # Fit Model
  M1_POSO <- jagsUI::jags(data = jags_data,
                          parameters.to.save = jags_params,
                          model.file = modelFilename,
                          n.chains = nc,
                          n.iter = ni,
                          n.burnin = nb,
                          n.thin = nt,
                          inits = init_vals)

  mcmc_output <- M1_POSO$samples

  summary_df <- as.data.frame(M1_POSO$summary)

  rhat_df <- summary_df[grepl("rhat\\[", rownames(summary_df)),]

  rhat_df$site <- unique_subunits

  output_file <- file.path(outputdir, "POSO_Trend_Results.csv")

  write.csv(rhat_df, file = output_file, row.names = FALSE)

  poso_trend <- ggplot(rhat_df, aes(x = site, y = mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "green", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "red", alpha = 0.1) +
    geom_point(stat = "identity", color = "black", size = 5) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "POSO Trend by Site (Green = Increasing, Red = Declining)",
      x = "Site",
      y = "Trend"
    ) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_file <- file.path(outputdir, "POSO_Trend_plot.png")

  ggsave(output_file, poso_trend, width = 10, height = 6, dpi = 300)




  ################################################################################
  # Jags Model
  # GIPU Models
  modelFilename = "GIPU.Removal.Model.txt"
  cat("
  model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions

  # Observation model: Each variable is conditionally distributed around its mean

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega + rhat[i] * (t - 1)
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  for(i in 1:nsites){
    rhat[i] ~ dunif(-1, 1)
  }

  for(i in 1:nsites){
    for(t in 1:nyears){
      latent.g[i, t] ~ dunif(1,npasses)
      g[i, t] <- round(latent.g[i, t])
      }
  }

  for(i in 1:nsites){
    p[i] ~ dbeta(1, 1)
  }

  for(i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:npasses){
        q[i, j, t] <- 1 - pow((1-p[i]), 1/g[i, t])
        #q[i, j, t] ~ dbeta(1, 1)
      }
    }
  }

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      p1_Y[i, t] ~ dbin(q[i, 1, t], N[i, t])
      p2_Y[i, t] ~ dbin(q[i, 2, t], (N[i, t] - p1_Y[i, t]))
      p3_Y[i, t] ~ dbin(q[i, 3, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t]))
      p4_Y[i, t] ~ dbin(q[i, 4, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t]))
      p5_Y[i, t] ~ dbin(q[i, 5, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t] - p4_Y[i, t]))
    }
  }
}
", fill = TRUE,file=modelFilename)




  # Parameters to save
  jags_params <- c(
    "omega",
    "rhat",
    "N"
  )

  # Initialize
  nsites <- dim(yGIPU)[1]
  npasses <- dim(yGIPU)[2]
  nyears <- dim(yGIPU)[3]

  N.inits <- array(NA, dim = c(nsites, nyears))
  p1_Y <- yGIPU[1:nsites,1,1:nyears]
  p2_Y <- yGIPU[1:nsites,2,1:nyears]
  p3_Y <- yGIPU[1:nsites,3,1:nyears]
  p4_Y <- yGIPU[1:nsites,4,1:nyears]
  p5_Y <- yGIPU[1:nsites,5,1:nyears]

  for (i in 1:nsites) {
    for (t in 1:nyears) {
      # Try summing the non-NA values across datasets
      summed_value <- sum(p1_Y[i, t], p2_Y[i, t], p3_Y[i, t],
                          p4_Y[i, t], p5_Y[i, t], na.rm = TRUE)

      # If all values were NA, the sum would be zero (or NA), so check this
      if (is.na(summed_value) || summed_value == 0) {
        # If the summed value is NA or zero, assign a Poisson value
        N.inits[i, t] <- round(rpois(1, lambda = 1))
      } else {
        # Otherwise, use the summed value multiplied by 10
        N.inits[i, t] <- round(summed_value)
      }
    }
  }

  # Bundle data
  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y
  )


  init_vals <- function() list(omega = rnorm(1, 0, 0.001),
                               N = N.inits)
  # MCMC settings
  ni <- 300000
  nc <- 3
  nb <- 20000
  nt <- 1

  set.seed(1234)

  # Fit Model
  M1_GIPU <- jagsUI::jags(data = jags_data,
                          parameters.to.save = jags_params,
                          model.file = modelFilename,
                          n.chains = nc,
                          n.iter = ni,
                          n.burnin = nb,
                          n.thin = nt,
                          inits = init_vals)

  mcmc_output <- M1_GIPU$samples

  summary_df <- as.data.frame(M1_GIPU$summary)

  rhat_df <- summary_df[grepl("rhat\\[", rownames(summary_df)),]

  rhat_df$site <- unique_subunits

  output_file <- file.path(outputdir, "GIPU_Trend_Results.csv")

  write.csv(rhat_df, file = output_file, row.names = FALSE)

  gipu_trend <- ggplot(rhat_df, aes(x = site, y = mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "green", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "red", alpha = 0.1) +
    geom_point(stat = "identity", color = "black", size = 5) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "GIPU Trend by Site (Green = Increasing, Red = Declining)",
      x = "Site",
      y = "Trend"
    ) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_file <- file.path(outputdir, "GIPU_Trend_plot.png")

  ggsave(output_file, gipu_trend, width = 10, height = 6, dpi = 300)



  # Jags Model
  # AGCH Models
  modelFilename = "AGCH.Removal.Model.txt"
  cat("
  model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions

  # Observation model: Each variable is conditionally distributed around its mean

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega + rhat[i] * (t - 1)
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  for(i in 1:nsites){
    rhat[i] ~ dunif(-1, 1)
  }

  for(i in 1:nsites){
      latent.g[i] ~ dunif(1,npasses)
      g[i] <- round(latent.g[i])
    }

  for(i in 1:nsites){
    p[i] ~ dbeta(1, 1)
    }


  for(i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:npasses){
        q[i, j, t] <- 1 - pow((1-p[i]), 1/g[i])
        #q[i, j, t] ~ dbeta(1, 1)
      }
    }
  }

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      p1_Y[i, t] ~ dbin(q[i, 1, t], N[i, t])
      p2_Y[i, t] ~ dbin(q[i, 2, t], (N[i, t] - p1_Y[i, t]))
      p3_Y[i, t] ~ dbin(q[i, 3, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t]))
      p4_Y[i, t] ~ dbin(q[i, 4, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t]))
      p5_Y[i, t] ~ dbin(q[i, 5, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t] - p4_Y[i, t]))
    }
  }
}
", fill = TRUE,file=modelFilename)



  # Parameters to save
  jags_params <- c(
    "omega",
    "rhat",
    "N"
  )


  # Initialize
  nsites <- dim(yAGCH)[1]
  npasses <- dim(yAGCH)[2]
  nyears <- dim(yAGCH)[3]

  N.inits <- array(NA, dim = c(nsites, nyears))
  p1_Y <- yAGCH[1:nsites,1,1:nyears]
  p2_Y <- yAGCH[1:nsites,2,1:nyears]
  p3_Y <- yAGCH[1:nsites,3,1:nyears]
  p4_Y <- yAGCH[1:nsites,4,1:nyears]
  p5_Y <- yAGCH[1:nsites,5,1:nyears]


  for (i in 1:nsites) {
    for (t in 1:nyears) {
      # Try summing the non-NA values across datasets
      summed_value <- sum(p1_Y[i, t], p2_Y[i, t], p3_Y[i, t],
                          p4_Y[i, t], p5_Y[i, t], na.rm = TRUE)

      # If all values were NA, the sum would be zero (or NA), so check this
      if (is.na(summed_value) || summed_value == 0) {
        # If the summed value is NA or zero, assign a Poisson value
        N.inits[i, t] <- round(rpois(1, lambda = 1))
      } else {
        # Otherwise, use the summed value multiplied by 10
        N.inits[i, t] <- round(summed_value)
      }
    }
  }

  # Bundle data
  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y
  )


  init_vals <- function() list(omega = rnorm(1, 0, 0.001),
                               N = N.inits)

  # MCMC settings
  ni <- 300000
  nc <- 3
  nb <- 20000
  nt <- 1

  # Fit Model
  M1_AGCH <- jagsUI::jags(data = jags_data,
                          parameters.to.save = jags_params,
                          model.file = modelFilename,
                          n.chains = nc,
                          n.iter = ni,
                          n.burnin = nb,
                          n.thin = nt,
                          inits = init_vals)

  mcmc_output <- M1_AGCH$samples

  summary_df <- as.data.frame(M1_AGCH$summary)

  rhat_df <- summary_df[grepl("rhat\\[", rownames(summary_df)),]

  rhat_df$site <- unique_subunits

  output_file <- file.path(outputdir, "AGCH_Trend_Results.csv")

  write.csv(rhat_df, file = output_file, row.names = FALSE)

  agch_trend <- ggplot(rhat_df, aes(x = site, y = mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "green", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "red", alpha = 0.1) +
    geom_point(stat = "identity", color = "black", size = 5) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "AGCH Trend by Site (Green = Increasing, Red = Declining)",
      x = "Site",
      y = "Trend"
    ) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_file <- file.path(outputdir, "AGCH_Trend_plot.png")

  ggsave(output_file, agch_trend, width = 10, height = 6, dpi = 300)




  # Jags Model
  # CYFO Models
  modelFilename = "CYFO.Removal.Model.txt"
  cat("
  model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions

  # Observation model: Each variable is conditionally distributed around its mean

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega + rhat[i] * (t - 1)
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  for(i in 1:nsites){
    rhat[i] ~ dunif(-1, 1)
  }

  for(i in 1:nsites){
      latent.g[i] ~ dunif(1,npasses)
      g[i] <- round(latent.g[i])
    }

  for(i in 1:nsites){
    p[i] ~ dbeta(1, 1)
    }


  for(i in 1:nsites){
    for(t in 1:nyears){
      for(j in 1:npasses){
        q[i, j, t] <- 1 - pow((1-p[i]), 1/g[i])
        #q[i, j, t] ~ dbeta(1, 1)
      }
    }
  }

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      p1_Y[i, t] ~ dbin(q[i, 1, t], N[i, t])
      p2_Y[i, t] ~ dbin(q[i, 2, t], (N[i, t] - p1_Y[i, t]))
      p3_Y[i, t] ~ dbin(q[i, 3, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t]))
      p4_Y[i, t] ~ dbin(q[i, 4, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t]))
      p5_Y[i, t] ~ dbin(q[i, 5, t], (N[i, t] - p1_Y[i, t] - p2_Y[i, t] - p3_Y[i, t] - p4_Y[i, t]))
    }
  }
}
", fill = TRUE,file=modelFilename)



  # Parameters to save
  jags_params <- c(
    "omega",
    "rhat",
    "N"
  )


  # Initialize
  nsites <- dim(yCYFO)[1]
  npasses <- dim(yCYFO)[2]
  nyears <- dim(yCYFO)[3]

  N.inits <- array(NA, dim = c(nsites, nyears))
  p1_Y <- yCYFO[1:nsites,1,1:nyears]
  p2_Y <- yCYFO[1:nsites,2,1:nyears]
  p3_Y <- yCYFO[1:nsites,3,1:nyears]
  p4_Y <- yCYFO[1:nsites,4,1:nyears]
  p5_Y <- yCYFO[1:nsites,5,1:nyears]


  for (i in 1:nsites) {
    for (t in 1:nyears) {
      # Try summing the non-NA values across datasets
      summed_value <- sum(p1_Y[i, t], p2_Y[i, t], p3_Y[i, t],
                          p4_Y[i, t], p5_Y[i, t], na.rm = TRUE)

      # If all values were NA, the sum would be zero (or NA), so check this
      if (is.na(summed_value) || summed_value == 0) {
        # If the summed value is NA or zero, assign a Poisson value
        N.inits[i, t] <- round(rpois(1, lambda = 1))
      } else {
        # Otherwise, use the summed value multiplied by 10
        N.inits[i, t] <- round(summed_value)
      }
    }
  }

  # Bundle data
  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y
  )


  init_vals <- function() list(omega = rnorm(1, 0, 0.001),
                               N = N.inits)

  # MCMC settings
  ni <- 300000
  nc <- 3
  nb <- 20000
  nt <- 1

  # Fit Model
  M1_CYFO <- jagsUI::jags(data = jags_data,
                          parameters.to.save = jags_params,
                          model.file = modelFilename,
                          n.chains = nc,
                          n.iter = ni,
                          n.burnin = nb,
                          n.thin = nt,
                          inits = init_vals)

  mcmc_output <- M1_CYFO$samples

  summary_df <- as.data.frame(M1_CYFO$summary)

  rhat_df <- summary_df[grepl("rhat\\[", rownames(summary_df)),]

  rhat_df$site <- unique_subunits

  output_file <- file.path(outputdir, "CYFO_Trend_Results.csv")

  write.csv(rhat_df, file = output_file, row.names = FALSE)

  cyfo_trend <- ggplot(rhat_df, aes(x = site, y = mean)) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = 0, ymax = Inf),
              fill = "green", alpha = 0.1) +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = 0),
              fill = "red", alpha = 0.1) +
    geom_point(stat = "identity", color = "black", size = 5) +
    geom_errorbar(aes(ymin = `2.5%`, ymax = `97.5%`), width = 0.2, color = "black") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "red", linewidth = 1) +
    labs(
      title = "CYFO Trend by Site (Green = Increasing, Red = Declining)",
      x = "Site",
      y = "Trend"
    ) +
    theme_bw(base_size = 18) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  output_file <- file.path(outputdir, "CYFO_Trend_plot.png")

  ggsave(output_file, cyfo_trend, width = 10, height = 6, dpi = 300)

}
