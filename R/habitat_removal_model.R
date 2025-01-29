#' Habitat Removal Model
#' A wrapper function run JAGS-based removal models for stream data.
#' Takes habitat, transect, and count data, processes it, runs the models for each species, and saves the output.
#'
#' @param hab Path to the habitat data file (CSV format).
#' @param tran Path to the transect habitat data file (CSV format)
#' @param count Path to the count data file (CSV format).
#' @param output_dir Directory to save the outputs.
#' @return A list containing count-habitat model summaries and plots.
#' @export
habitat_removal_model <- function(hab, tran, count, output_dir){
  library(dplyr)
  library(tidyr)
  library(jagsUI)
  library(ggplot2)
  library(coda)

  # Step 1: Load data
  hab <- hab
  tran <- tran
  count <- count

  check_data_consistency <- function(hab, tran, count) {
    # Unique sites and years in each dataset
    hab_sites <- unique(hab$Subunit)
    tran_sites <- unique(tran$Subunit)
    count_sites <- unique(count$Subunit)

    hab_years <- unique(hab$Year)
    tran_years <- unique(tran$Year)
    count_years <- unique(count$Year)

    # Check if all datasets have the same sites
    if (!setequal(hab_sites, tran_sites) || !setequal(hab_sites, count_sites)) {
      stop("Error: Site names (Subunit) are inconsistent across datasets.\n",
           "Habitat Sites: ", paste(hab_sites, collapse = ", "), "\n",
           "Transect Sites: ", paste(tran_sites, collapse = ", "), "\n",
           "Count Sites: ", paste(count_sites, collapse = ", "))
    }

    # Check if all datasets have the same years
    if (!setequal(hab_years, tran_years) || !setequal(hab_years, count_years)) {
      stop("Error: Years are inconsistent across datasets.\n",
           "Habitat Years: ", paste(hab_years, collapse = ", "), "\n",
           "Transect Years: ", paste(tran_years, collapse = ", "), "\n",
           "Count Years: ", paste(count_years, collapse = ", "))
    }

    # Ensure complete `Subunit-Year` combinations
    hab_combinations <- expand.grid(Subunit = hab_sites, Year = hab_years)
    tran_combinations <- expand.grid(Subunit = tran_sites, Year = tran_years)
    count_combinations <- expand.grid(Subunit = count_sites, Year = count_years)

    if (!setequal(hab_combinations, tran_combinations) ||
        !setequal(hab_combinations, count_combinations)) {
      missing_hab <- setdiff(tran_combinations, hab_combinations)
      missing_tran <- setdiff(hab_combinations, tran_combinations)
      missing_count <- setdiff(hab_combinations, count_combinations)

      stop("Error: Inconsistent Subunit-Year combinations across datasets.\n",
           if (nrow(missing_hab) > 0) paste("Missing in Habitat Data: ", paste(missing_hab, collapse = ", "), "\n") else "",
           if (nrow(missing_tran) > 0) paste("Missing in Transect Data: ", paste(missing_tran, collapse = ", "), "\n") else "",
           if (nrow(missing_count) > 0) paste("Missing in Count Data: ", paste(missing_count, collapse = ", "), "\n") else "")
    }

    message("Data consistency checks passed!")
  }

  check_data_consistency(hab,tran,count)






  # Step 2: Preprocess habitat data
  hab_clean <- hab %>%
    group_by(Year, Subunit) %>%
    filter(!all(is.na(select(., starts_with("Velocity")))) &
             !all(is.na(select(., starts_with("Depth")))))

  hab_clean_filtered <- hab_clean %>%
    group_by(Year, Subunit) %>%
    filter(!all(is.na(Velocity1) & is.na(Velocity2) & is.na(Velocity3) & is.na(Velocity4))) %>%
    ungroup()  # Ungroup to avoid grouping-related issues

  # Calculate each metric with NAs removed by transect
  hab_summary <- hab_clean_filtered %>%
    mutate(
      MeanFlow = rowMeans(select(., starts_with("Velocity")), na.rm = TRUE),
      MeanDepth = rowMeans(select(., starts_with("Depth")), na.rm = TRUE),
      CrossSectionalArea = ChannelWidth * MeanDepth,
      Volume = CrossSectionalArea * 25, # 25 m channel unit length (m)
      ChannelUnitArea = ChannelWidth * 25 # 25 m channel unit length (m)
    )

  # Calculate each metric and remaining variables by Subunit and Year
  hab_summary <- hab_summary %>%
    group_by(Subunit, Year) %>%
    summarise(
      Mean_ChannelWidth = mean(ChannelWidth, na.rm = TRUE),
      Mean_Temp = mean(TransectTemp, na.rm = TRUE),
      Mean_Conductivity = mean(TransectConductivity, na.rm = TRUE),
      Mean_DO = mean(TransectDissolvedOxygen, na.rm = TRUE),
      Mean_pH = mean(TransectpH, na.rm = TRUE),
      Mean_Fine = mean(Fine, na.rm = TRUE),
      Mean_Gravel = mean(Gravel, na.rm = TRUE),
      Mean_Cobble = mean(Cobble, na.rm = TRUE),
      Mean_Boulder = mean(Boulder, na.rm = TRUE),
      Mean_Bedrock = mean(Bedrock , na.rm = TRUE),
      Mean_ChannelUnit = mean(ChannelUnit, na.rm = TRUE),
      Mean_MeanFlow = mean(MeanFlow, na.rm = TRUE),
      Mean_MeanDepth = mean(MeanDepth, na.rm = TRUE),
      Mean_CrossSectionalArea = mean(CrossSectionalArea, na.rm = TRUE),
      Mean_ChannelUnitArea = mean(ChannelUnitArea, na.rm = TRUE),
      Total_Volume = sum(Volume, na.rm = TRUE)
    )

  year_range <- range(hab_summary$Year, na.rm = TRUE)

  all_combinations <- expand.grid(
    Subunit = unique(hab_summary$Subunit),
    Year = seq(year_range[1], year_range[2])
  )

  hab_summary <- all_combinations %>%
    left_join(hab_summary, by = c("Subunit", "Year")) %>%
    arrange(Subunit, Year)


  # Create a matrix for each variable
  unique_subunits <- unique(hab_summary$Subunit)
  unique_years <- unique(hab_summary$Year)


  channelwidth <- wtemp <- wcond <- doxygen <- pH <- fine <- gravel <- cobble <- boulder <- bedrock <- channelunit <- flow <- wdepth <- crossSectionalArea <- channelunitArea <- volume <- matrix(NA, nrow = length(unique_subunits), ncol = length(unique_years))

  for (i in seq_along(unique_subunits)) {
    for (j in seq_along(unique_years)) {
      subunit <- unique_subunits[i]
      year <- unique_years[j]

      # Filter the data for the specific Subunit and Year
      data_row <- hab_summary %>%
        filter(Subunit == subunit, Year == year)

      # Fill in the values if the row exists
      if (nrow(data_row) > 0) {
        channelwidth[i, j] <- data_row$Mean_ChannelWidth
        wtemp[i, j] <- data_row$Mean_Temp
        wcond[i, j] <- data_row$Mean_Conductivity
        doxygen[i, j] <- data_row$Mean_DO
        pH[i, j] <- data_row$Mean_pH
        fine[i, j] <- data_row$Mean_Fine
        gravel[i, j] <- data_row$Mean_Gravel
        cobble[i, j] <- data_row$Mean_Cobble
        boulder[i, j] <- data_row$Mean_Boulder
        bedrock[i, j] <- data_row$Mean_Bedrock
        channelunit[i, j] <- data_row$Mean_ChannelUnit
        flow[i, j] <- data_row$Mean_MeanFlow
        wdepth[i, j] <- data_row$Mean_MeanDepth
        crossSectionalArea[i, j] <- data_row$Mean_CrossSectionalArea
        channelunitArea[i, j] <- data_row$Mean_ChannelUnitArea
        volume[i, j] <- data_row$Total_Volume
      }
    }
  }


  # Format transect/stream reach data
  unique_years <- unique(tran$Year)
  unique_subunits <- unique(tran$Subunit)

  netters <- woodpieces <- wturb <- efreq <- evolt <- matrix(NA, nrow = length(unique_subunits), ncol = length(unique_years))

  for(i in seq_along(unique_subunits)){
    for(t in seq_along(unique_years)){
      subunit <- unique_subunits[i]
      year <- unique_years[t]

      data_row <- tran %>%
        filter(Subunit == subunit, Year == year)

      if(nrow(data_row)>0){
        netters[i,t] <- data_row$NumberofNetters
        woodpieces[i,t] <- data_row$NumberofWoodPieces / 25
        wturb[i,t] <- data_row$Turbidity
        efreq[i,t] <- data_row$Frequency
        evolt[i,t] <- data_row$Voltage
      }
    }
  }

  # Scale variables
  channelwidth <- apply(channelwidth, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  wtemp <- apply(wtemp, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  wcond <- apply(wcond, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  doxygen <- apply(doxygen, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  pH <- apply(pH, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  fine <- apply(fine, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  gravel <- apply(gravel, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  cobble <- apply(cobble, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  boulder <- apply(boulder, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  bedrock <- apply(bedrock, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  channelunit <- apply(channelunit, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  flow <- apply(flow, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  wdepth <- apply(wdepth, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  crossSectionalArea <- apply(crossSectionalArea, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  channelunitArea <- apply(channelunitArea, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  volume <- apply(volume, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

  netters <- apply(netters, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  woodpieces <- apply(woodpieces, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  wturb <- apply(wturb, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  efreq <- apply(efreq, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))
  evolt <- apply(evolt, 2, function(x)(x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE))

  # Aggregate coordinates for each Subunit
  subunit_coords <- tran %>%
    group_by(Subunit) %>%
    summarise(
      mean_latitude = mean(Latitude, na.rm = TRUE),
      mean_longitude = mean(Longitude, na.rm = TRUE)
    )

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
  modelFilename = "POSO.Habitat.Removal.Model.txt"
  cat("
model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions
sigma_channelwidth ~ dunif(0, 100)
tau_channelwidth <- 1 / pow(sigma_channelwidth, 2)

sigma_flow ~ dunif(0, 100)
tau_flow <- 1 / pow(sigma_flow, 2)

sigma_wtemp ~ dunif(0, 100)
tau_wtemp <- 1 / pow(sigma_wtemp, 2)

sigma_wcond ~ dunif(0, 100)
tau_wcond <- 1 / pow(sigma_wcond, 2)

sigma_pH ~ dunif(0, 100)
tau_pH <- 1 / pow(sigma_pH, 2)

sigma_cobble ~ dunif(0, 100)
tau_cobble <- 1 / pow(sigma_cobble, 2)

sigma_wdepth ~ dunif(0, 100)
tau_wdepth <- 1 / pow(sigma_wdepth, 2)

sigma_volume ~ dunif(0, 100)
tau_volume <- 1 / pow(sigma_volume, 2)

sigma_channelunit ~ dunif(0, 100)
tau_channelunit <- 1 / pow(sigma_channelunit, 2)

# Observation model: Each variable is conditionally distributed around its mean
for (i in 1:nsites) {
  for (t in 1:nyears) {
    # Variables and their observation-level likelihoods
    channelwidth[i, t] ~ dnorm(mu_channelwidth[i, t], tau_channelwidth)
    flow[i, t] ~ dnorm(mu_flow[i, t], tau_flow)
    wtemp[i, t] ~ dnorm(mu_wtemp[i, t], tau_wtemp)
    wcond[i, t] ~ dnorm(mu_wcond[i, t], tau_wcond)
    pH[i, t] ~ dnorm(mu_pH[i, t], tau_pH)
    cobble[i, t] ~ dnorm(mu_cobble[i, t], tau_cobble)

    wdepth[i, t] ~ dnorm(mu_wdepth[i, t], tau_wdepth)
    volume[i, t] ~ dnorm(mu_volume[i, t], tau_volume)
    channelunit[i, t] ~ dnorm(mu_channelunit[i, t], tau_channelunit)

  }
}


  # Priors for the mean parameters
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      mu_channelwidth[i, t] ~ dnorm(0, 0.001)
      mu_flow[i, t] ~ dnorm(0, 0.001)
      mu_wtemp[i, t] ~ dnorm(0, 0.001)
      mu_wcond[i, t] ~ dnorm(0, 0.001)
      mu_pH[i, t] ~ dnorm(0, 0.001)
      mu_cobble[i, t] ~ dnorm(0, 0.001)

      mu_wdepth[i, t] ~ dnorm(0, 0.001)
      mu_volume[i, t] ~ dnorm(0, 0.001)
      mu_channelunit[i, t] ~ dnorm(0, 0.001)
    }
  }

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega +
                          beta.channelwidth * channelwidth[i, t] +
                          beta.flow * flow[i, t] +
                          beta.wtemp * wtemp[i, t] +
                          beta.wcond * wcond[i, t] +
                          beta.cobble * cobble[i, t]
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  beta.channelwidth ~ dnorm(0, 0.01)
  beta.flow ~ dnorm(0, 0.01)
  beta.wtemp ~ dnorm(0, 0.01)
  beta.wcond ~ dnorm(0, 0.01)
  beta.cobble ~ dnorm(0, 0.01)

  theta0 ~ dnorm(0, 0.01)
  theta.depth ~ dnorm(0, 0.01)
  theta.volume ~ dnorm(0, 0.01)
  theta.chunit ~ dnorm(0, 0.01)
  theta.pH ~ dunif(0, 0.01)

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (j in 1:npasses) {
        logit(q[i, j, t]) <- theta0 + theta.depth * wdepth[i, t] + theta.volume * volume[i, t] + theta.chunit * channelunit[i, t] + theta.pH * pH[i, t]
      }
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
    "beta.channelwidth",
    "beta.flow",
    "beta.wtemp",
    "beta.wcond",
    "beta.cobble",
    "theta0",
    "theta.depth",
    "theta.volume",
    "theta.chunit",
    "theta.pH"
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
    p5_Y = p5_Y,
    channelwidth = channelwidth,
    flow = flow,
    wtemp = wtemp,
    wcond = wcond,
    pH = pH,
    cobble = cobble,
    wdepth = wdepth,
    volume = volume,
    channelunit = channelunit
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

  out2 <- capture.output(print(M1_POSO, digits=3))

  output_file <- file.path(outputdir, "POSO_Habitat_Results.csv")

  write.csv(out2, file=output_file, row.names = FALSE, quote = FALSE)

  mcmc_output <- M1_POSO$samples

  beta_params <- mcmc_output[, grep("^beta", varnames(mcmc_output))]
  theta_params <- mcmc_output[, grep("^theta", varnames(mcmc_output))]


  # Function to calculate means and credible intervals
  summary_stats <- function(mcmc_list) {
    # Combine all chains into a single matrix
    combined <- do.call(rbind, lapply(mcmc_list, as.data.frame))

    # Calculate summary statistics for each parameter
    summary <- apply(combined, 2, function(x) {
      mean_x <- mean(x)
      lower <- quantile(x, 0.025)
      upper <- quantile(x, 0.975)
      c(mean = mean_x, lower = lower, upper = upper)
    })

    # Convert to a data frame
    summary_df <- as.data.frame(t(summary))
    summary_df$Parameter <- rownames(summary_df)
    summary_df <- summary_df %>% arrange(desc(mean))
    rownames(summary_df) <- NULL
    return(summary_df)
  }

  # Get summary statistics for beta parameters
  beta_summary <- summary_stats(beta_params)

  colnames(beta_summary) <- c("mean", "lower", "upper", "Parameter")

  # Coefficient plot
  p_beta_poso <- ggplot(beta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Beta Parameters for POSO",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(outputdir, "POSO_beta_coefficient_plot.png")

  ggsave(output_file, p_beta_poso, width = 8, height = 6, dpi = 300)

  theta_summary <- summary_stats(theta_params)
  colnames(theta_summary) <- c("mean", "lower", "upper", "Parameter")

  theta_summary <- theta_summary[theta_summary$Parameter != "theta0", ]

  # Coefficient plot
  p_theta_poso <- ggplot(theta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Theta Parameters for POSO",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(outputdir, "POSO_theta_coefficient_plot.png")

  ggsave(output_file, p_theta_poso, width = 8, height = 6, dpi = 300)


  ################################################################################
  # Jags Model
  # GIPU Models
  modelFilename = "GIPU.Habitat.Removal.Model.txt"
  cat("
model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions
sigma_gravel ~ dunif(0, 100)
tau_gravel <- 1 / pow(sigma_gravel, 2)

sigma_flow ~ dunif(0, 100)
tau_flow <- 1 / pow(sigma_flow, 2)

sigma_wtemp ~ dunif(0, 100)
tau_wtemp <- 1 / pow(sigma_wtemp, 2)

sigma_pH ~ dunif(0, 100)
tau_pH <- 1 / pow(sigma_pH, 2)

sigma_cobble ~ dunif(0, 100)
tau_cobble <- 1 / pow(sigma_cobble, 2)

sigma_wdepth ~ dunif(0, 100)
tau_wdepth <- 1 / pow(sigma_wdepth, 2)

sigma_volume ~ dunif(0, 100)
tau_volume <- 1 / pow(sigma_volume, 2)

sigma_channelunit ~ dunif(0, 100)
tau_channelunit <- 1 / pow(sigma_channelunit, 2)

sigma_boulder ~ dunif(0, 100)
tau_boulder <- 1 / pow(sigma_boulder, 2)

# Observation model: Each variable is conditionally distributed around its mean
for (i in 1:nsites) {
  for (t in 1:nyears) {
    # Variables and their observation-level likelihoods
    boulder[i, t] ~ dnorm(mu_boulder[i, t], tau_boulder)
    flow[i, t] ~ dnorm(mu_flow[i, t], tau_flow)
    wtemp[i, t] ~ dnorm(mu_wtemp[i, t], tau_wtemp)
    gravel[i, t] ~ dnorm(mu_gravel[i, t], tau_gravel)
    pH[i, t] ~ dnorm(mu_pH[i, t], tau_pH)
    cobble[i, t] ~ dnorm(mu_cobble[i, t], tau_cobble)
    wdepth[i, t] ~ dnorm(mu_wdepth[i, t], tau_wdepth)
    volume[i, t] ~ dnorm(mu_volume[i, t], tau_volume)
    channelunit[i, t] ~ dnorm(mu_channelunit[i, t], tau_channelunit)

  }
}


  # Priors for the mean parameters
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      mu_boulder[i, t] ~ dnorm(0, 0.001)
      mu_flow[i, t] ~ dnorm(0, 0.001)
      mu_wtemp[i, t] ~ dnorm(0, 0.001)
      mu_gravel[i, t] ~ dnorm(0, 0.001)
      mu_pH[i, t] ~ dnorm(0, 0.001)
      mu_cobble[i, t] ~ dnorm(0, 0.001)
      mu_wdepth[i, t] ~ dnorm(0, 0.001)
      mu_volume[i, t] ~ dnorm(0, 0.001)
      mu_channelunit[i, t] ~ dnorm(0, 0.001)
    }
  }

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega +
                          beta.flow * flow[i, t] +
                          beta.boulder * boulder[i, t] +
                          beta.gravel * gravel[i, t] +
                          beta.wtemp * wtemp[i, t] +
                          beta.cobble * cobble[i, t]
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  beta.flow ~ dnorm(0, 0.01)
  beta.boulder ~ dnorm(0, 0.01)
  beta.gravel ~ dnorm(0, 0.01)
  beta.wtemp ~ dnorm(0, 0.01)
  beta.cobble ~ dnorm(0, 0.01)

  theta0 ~ dnorm(0, 0.01)
  theta.depth ~ dnorm(0, 0.01)
  theta.volume ~ dnorm(0, 0.01)
  theta.chunit ~ dnorm(0, 0.01)
  theta.pH ~ dunif(0, 0.01)

  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (j in 1:npasses) {
        logit(q[i, j, t]) <- theta0 + theta.depth * wdepth[i, t] + theta.volume * volume[i, t] + theta.chunit * channelunit[i, t]
      }
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
    "beta.boulder",
    "beta.flow",
    "beta.wtemp",
    "beta.cobble",
    "beta.gravel",
    "theta0",               # Save Bayesian indicators
    "theta.depth",
    "theta.volume",
    "theta.chunit"
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

  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y,
    flow = flow,
    gravel = gravel,
    boulder = boulder,
    wtemp = wtemp,
    pH = pH,
    cobble = cobble,
    wdepth = wdepth,
    volume = volume,
    channelunit = channelunit
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

  out2 <- capture.output(print(M1_GIPU, digits = 3))

  output_file <- file.path(output_dir, "GIPU_Habitat_Results.csv")

  write.csv(out2, file=output_file, row.names = FALSE, quote = FALSE)

  mcmc_output <- M1_GIPU$samples

  beta_params <- mcmc_output[, grep("^beta", varnames(mcmc_output))]
  theta_params <- mcmc_output[, grep("^theta", varnames(mcmc_output))]

  # Function to calculate means and credible intervals
  summary_stats <- function(mcmc_list) {
    # Combine all chains into a single matrix
    combined <- do.call(rbind, lapply(mcmc_list, as.data.frame))

    # Calculate summary statistics for each parameter
    summary <- apply(combined, 2, function(x) {
      mean_x <- mean(x)
      lower <- quantile(x, 0.025)
      upper <- quantile(x, 0.975)
      c(mean = mean_x, lower = lower, upper = upper)
    })

    # Convert to a data frame
    summary_df <- as.data.frame(t(summary))
    summary_df$Parameter <- rownames(summary_df)
    summary_df <- summary_df %>% arrange(desc(mean))
    rownames(summary_df) <- NULL
    return(summary_df)
  }

  # Get summary statistics for beta parameters
  beta_summary <- summary_stats(beta_params)

  colnames(beta_summary) <- c("mean", "lower", "upper", "Parameter")

  # Coefficient plot
  p_beta_gipu <- ggplot(beta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Beta Parameters for GIPU",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(output_dir, "GIPU_beta_coefficient_plot.png")

  ggsave(output_file, p_beta_gipu, width = 8, height = 6, dpi = 300)

  theta_summary <- summary_stats(theta_params)
  colnames(theta_summary) <- c("mean", "lower", "upper", "Parameter")

  theta_summary <- theta_summary[theta_summary$Parameter != "theta0", ]

  # Coefficient plot
  p_theta_gipu <- ggplot(theta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Theta Parameters for GIPU",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(output_dir, "GIPU_theta_coefficient_plot.png")

  ggsave(output_file, p_theta_gipu, width = 8, height = 6, dpi = 300)


  # Jags Model
  # AGCH Models
  modelFilename = "AGCH.Habitat.Removal.Model.txt"
  cat("
  model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions
sigma_doxygen ~ dunif(0, 100)
tau_doxygen <- 1 / pow(sigma_doxygen, 2)

sigma_fine ~ dunif(0, 100)
tau_fine <- 1 / pow(sigma_fine, 2)

sigma_wdepth ~ dunif(0, 100)
tau_wdepth <- 1 / pow(sigma_wdepth, 2)

sigma_volume ~ dunif(0, 100)
tau_volume <- 1 / pow(sigma_volume, 2)

sigma_channelunit ~ dunif(0, 100)
tau_channelunit <- 1 / pow(sigma_channelunit, 2)



# Observation model: Each variable is conditionally distributed around its mean
for (i in 1:nsites) {
  for (t in 1:nyears) {
    # Variables and their observation-level likelihoods
    doxygen[i, t] ~ dnorm(mu_doxygen[i, t], tau_doxygen)
    fine[i, t] ~ dnorm(mu_fine[i, t], tau_fine)
    wdepth[i, t] ~ dnorm(mu_wdepth[i, t], tau_wdepth)
    volume[i, t] ~ dnorm(mu_volume[i, t], tau_volume)
    channelunit[i, t] ~ dnorm(mu_channelunit[i, t], tau_channelunit)
  }
}


  # Priors for the mean parameters
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      mu_doxygen[i, t] ~ dnorm(0, 0.001)
      mu_fine[i, t] ~ dnorm(0, 0.001)
      mu_wdepth[i, t] ~ dnorm(0, 0.001)
      mu_volume[i, t] ~ dnorm(0, 0.001)
      mu_channelunit[i, t] ~ dnorm(0, 0.001)
    }
  }

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega +
                          beta.doxygen * doxygen[i, t] +
                          beta.fine * fine[i, t] +
                          beta.volume * volume[i, t]
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  beta.doxygen ~ dnorm(0, 0.01)
  beta.fine ~ dnorm(0, 0.01)
  beta.volume ~ dnorm(0, 0.01)

  theta0 ~ dnorm(0, 0.01)
  theta.depth ~ dnorm(0, 0.01)
  theta.volume ~ dnorm(0, 0.01)
  theta.chunit ~ dnorm(0, 0.01)


  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (j in 1:npasses) {
        logit(q[i, j, t]) <- theta0 + theta.depth * wdepth[i, t] + theta.volume * volume[i, t] + theta.chunit * channelunit[i, t]
      }
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
    "beta.doxygen",
    "beta.fine",
    "beta.volume",
    "theta0",               # Save Bayesian indicators
    "theta.depth",
    "theta.volume",
    "theta.chunit"
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

  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y,
    doxygen = doxygen,
    fine = fine,
    wdepth = wdepth,
    volume = volume,
    channelunit = channelunit
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

  out2 <- capture.output(print(M1_AGCH, digits = 3))

  output_file <- file.path(output_dir, "AGCH_Habitat_Results.csv")

  write.csv(out2, file=output_file, row.names = FALSE, quote = FALSE)

  mcmc_output <- M1_AGCH$samples

  beta_params <- mcmc_output[, grep("^beta", varnames(mcmc_output))]
  theta_params <- mcmc_output[, grep("^theta", varnames(mcmc_output))]

  # Function to calculate means and credible intervals
  summary_stats <- function(mcmc_list) {
    # Combine all chains into a single matrix
    combined <- do.call(rbind, lapply(mcmc_list, as.data.frame))

    # Calculate summary statistics for each parameter
    summary <- apply(combined, 2, function(x) {
      mean_x <- mean(x)
      lower <- quantile(x, 0.025)
      upper <- quantile(x, 0.975)
      c(mean = mean_x, lower = lower, upper = upper)
    })

    # Convert to a data frame
    summary_df <- as.data.frame(t(summary))
    summary_df$Parameter <- rownames(summary_df)
    summary_df <- summary_df %>% arrange(desc(mean))
    rownames(summary_df) <- NULL
    return(summary_df)
  }

  # Get summary statistics for beta parameters
  beta_summary <- summary_stats(beta_params)

  colnames(beta_summary) <- c("mean", "lower", "upper", "Parameter")

  # Coefficient plot
  p_beta_agch <- ggplot(beta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Beta Parameters for AGCH",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars


  output_file <- file.path(output_dir, "AGCH_beta_coefficient_plot.png")

  ggsave(output_file, p_beta_agch, width = 8, height = 6, dpi = 300)

  theta_summary <- summary_stats(theta_params)
  colnames(theta_summary) <- c("mean", "lower", "upper", "Parameter")

  theta_summary <- theta_summary[theta_summary$Parameter != "theta0", ]

  # Coefficient plot
  p_theta_agch <- ggplot(theta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Theta Parameters for AGCH",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars


  output_file <- file.path(output_dir, "AGCH_theta_coefficient_plot.png")

  ggsave(output_file, p_theta_agch, width = 8, height = 6, dpi = 300)



  # Jags Model
  # CYFO Models
  modelFilename = "CYFO.Habitat.Removal.Model.txt"
  cat("
model {

  # Priors for hyperparameters (precision parameters)
  # Priors for variance components and their precisions
sigma_flow ~ dunif(0, 100)
tau_flow <- 1 / pow(sigma_flow, 2)

sigma_wcond ~ dunif(0, 100)
tau_wcond <- 1 / pow(sigma_wcond, 2)

sigma_cobble ~ dunif(0, 100)
tau_cobble <- 1 / pow(sigma_cobble, 2)

sigma_wdepth ~ dunif(0, 100)
tau_wdepth <- 1 / pow(sigma_wdepth, 2)

sigma_volume ~ dunif(0, 100)
tau_volume <- 1 / pow(sigma_volume, 2)

sigma_channelunit ~ dunif(0, 100)
tau_channelunit <- 1 / pow(sigma_channelunit, 2)



# Observation model: Each variable is conditionally distributed around its mean
for (i in 1:nsites) {
  for (t in 1:nyears) {
    # Variables and their observation-level likelihoods
    flow[i, t] ~ dnorm(mu_flow[i, t], tau_flow)
    wcond[i, t] ~ dnorm(mu_wcond[i, t], tau_wcond)
    cobble[i, t] ~ dnorm(mu_cobble[i, t], tau_cobble)
    wdepth[i, t] ~ dnorm(mu_wdepth[i, t], tau_wdepth)
    volume[i, t] ~ dnorm(mu_volume[i, t], tau_volume)
    channelunit[i, t] ~ dnorm(mu_channelunit[i, t], tau_channelunit)
  }
}


  # Priors for the mean parameters
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      mu_flow[i, t] ~ dnorm(0, 0.001)
      mu_wcond[i, t] ~ dnorm(0, 0.01)
      mu_cobble[i, t] ~ dnorm(0, 0.01)
      mu_wdepth[i, t] ~ dnorm(0, 0.001)
      mu_volume[i, t] ~ dnorm(0, 0.001)
      mu_channelunit[i, t] ~ dnorm(0, 0.001)
    }
  }

  # Process model for abundance
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      N[i, t] ~ dpois(lambda[i, t])
      log(lambda[i, t]) <- omega +
                          beta.flow * flow[i, t] +
                          beta.wcond * wcond[i, t] +
                          beta.cobble * cobble[i, t]
    }
  }

  # Priors for process model parameters
  omega ~ dnorm(0, 0.01)
  beta.flow ~ dnorm(0, 0.01)
  beta.wcond ~ dnorm(0, 0.01)
  beta.cobble ~ dnorm(0, 0.01)
  theta0 ~ dnorm(0, 0.01)
  theta.depth ~ dnorm(0, 0.01)
  theta.volume ~ dnorm(0, 0.01)
  theta.chunit ~ dnorm(0, 0.01)


  # Observation model for removal data
  for (i in 1:nsites) {
    for (t in 1:nyears) {
      for (j in 1:npasses) {
        logit(q[i, j, t]) <- theta0 + theta.depth * wdepth[i, t] + theta.volume * volume[i, t] + theta.chunit * channelunit[i, t]
      }
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
    "beta.flow",
    "beta.wcond",
    "beta.cobble",
    "theta0",
    "theta.depth",
    "theta.volume",
    "theta.chunit"
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

  jags_data <- list(
    nsites = nsites,
    nyears = nyears,
    npasses = npasses,
    p1_Y = p1_Y,
    p2_Y = p2_Y,
    p3_Y = p3_Y,
    p4_Y = p4_Y,
    p5_Y = p5_Y,
    flow = flow,
    wcond = wcond,
    cobble = cobble,
    wdepth = wdepth,
    volume = volume,
    channelunit = channelunit
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

  out2 <- capture.output(print(M1_CYFO, digits = 3))

  output_file <- file.path(output_dir, "CYFO_Habitat_Results.csv")

  write.csv(out2, file=output_file, row.names = FALSE, quote = FALSE)

  mcmc_output <- M1_CYFO$samples

  beta_params <- mcmc_output[, grep("^beta", varnames(mcmc_output))]
  theta_params <- mcmc_output[, grep("^theta", varnames(mcmc_output))]

  # Function to calculate means and credible intervals
  summary_stats <- function(mcmc_list) {
    # Combine all chains into a single matrix
    combined <- do.call(rbind, lapply(mcmc_list, as.data.frame))

    # Calculate summary statistics for each parameter
    summary <- apply(combined, 2, function(x) {
      mean_x <- mean(x)
      lower <- quantile(x, 0.025)
      upper <- quantile(x, 0.975)
      c(mean = mean_x, lower = lower, upper = upper)
    })

    # Convert to a data frame
    summary_df <- as.data.frame(t(summary))
    summary_df$Parameter <- rownames(summary_df)
    summary_df <- summary_df %>% arrange(desc(mean))
    rownames(summary_df) <- NULL
    return(summary_df)
  }

  # Get summary statistics for beta parameters
  beta_summary <- summary_stats(beta_params)

  colnames(beta_summary) <- c("mean", "lower", "upper", "Parameter")

  # Coefficient plot
  p_beta_cyfo <- ggplot(beta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Beta Parameters for CYFO",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(output_dir, "CYFO_beta_coefficient_plot.png")

  ggsave(output_file, p_beta_cyfo, width = 8, height = 6, dpi = 300)

  theta_summary <- summary_stats(theta_params)
  colnames(theta_summary) <- c("mean", "lower", "upper", "Parameter")

  theta_summary <- theta_summary[theta_summary$Parameter != "theta0", ]

  # Coefficient plot
  p_theta_cyfo <- ggplot(theta_summary, aes(x = reorder(Parameter, mean), y = mean)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "black", linewidth = 1) +
    geom_point(size = 4, color = "blue") +  # Larger points with color
    geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.3, color = "black") +  # Black error bars
    theme_classic(base_size = 14) +  # Increase base font size for readability
    theme(
      axis.text = element_text(size = 12),  # Increase axis text size
      axis.title = element_text(size = 14, face = "bold"),  # Bold axis titles
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5)  # Center the title
    ) +
    labs(
      title = "Posterior Estimates of Theta Parameters for CYFO",
      x = "Parameter",
      y = "Posterior Mean (95% Credible Interval)"
    ) +
    coord_flip()  # Flip coordinates for horizontal error bars

  output_file <- file.path(output_dir, "CYFO_theta_coefficient_plot.png")

  ggsave(output_file, p_theta_cyfo, width = 8, height = 6, dpi = 300)


}
