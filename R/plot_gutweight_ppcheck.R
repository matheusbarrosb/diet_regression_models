plot_gutweight_ppcheck <- function(
    post,           
    stan_data,      
    spps = c("Pinfish", "Croaker", "Silver perch"),
    palette = "Bay",
    title,
    upperX = 0.4,
    strip_text_size
) {
  library(dplyr)
  library(ggplot2)
  
  n_iter <- dim(post$alpha)[1]
  n_sites <- dim(post$beta_site)[3]
  n_species <- length(spps)
  
  # Use representative TL and FW (mean from stan_data)
  rep_TL <- mean(stan_data$TL, na.rm = TRUE)
  rep_FW <- mean(stan_data$FW, na.rm = TRUE)
  
  # Calculate mean site effect for each iteration and species
  mean_site_effect <- array(NA, dim = c(n_iter, n_species))
  for (s in seq_along(spps)) {
    mean_site_effect[, s] <- apply(post$beta_site[, s, ], 1, mean)
  }
  
  # Generate posterior predictions for each species (average over status)
  pred_df <- data.frame()
  for (s in seq_along(spps)) {
    # Linear predictor for each posterior draw: average over statuses
    lp <- post$alpha[, s] +
      mean_site_effect[, s] +
      apply(post$beta_RS[, s, , drop=FALSE], 1, mean) + # mean over statuses
      post$beta_TL[, s] * rep_TL +
      post$beta_FW[, s] * rep_FW
    pred_df <- rbind(pred_df, data.frame(
      species = spps[s],
      draw = 1:n_iter,
      log_GW = lp,
      GW = exp(lp)
    ))
  }
  
  # Observed data: GW and species
  obs_df <- data.frame(
    GW = stan_data$GW,
    species = spps[stan_data$species]
  )
  
  # Plot densities
  ggplot() +
    geom_density(data = pred_df, aes(x = GW, fill = "Predicted"), alpha = 0.5, fill = pnw_palette(palette, 2)[1]) +
    geom_density(data = obs_df, aes(x = GW, fill = "Observed"), alpha = 0.5, fill = pnw_palette(palette, 2)[2]) +
    facet_wrap(~ species, scales = "free_y") +
#    scale_fill_manual(values = pnw_palette(palette, 2)) +
    labs(
      title = title,
      x = "Gut Weight",
      y = "Density",
      fill = ""
    ) +
    custom_theme() +
    theme(legend.position = "top",
          strip.text = element_blank()) +
    xlim(0, upperX)
}
