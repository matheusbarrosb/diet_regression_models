plot_betas <- function(
    post,
    spps = c("Pinfish", "Croaker", "Silver perch"),
    site_labels = c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"),
    status_labels = c("CT", "LS"),
    param_panels = c("beta_RS", "beta_site", "beta_TL", "beta_FW", "beta_site_RS"),
    title = "",
    palette = "Bay",
    ncol
) {
  
  # Setup species colors
  n_species <- length(spps)
  species_colors <- pnw_palette(palette, n = n_species)
  
  # Site presence rules (site names!)
  present_sites <- list(
    "Pinfish" = site_labels,
    "Croaker" = c("CI", "NEPaP", "SA"),
    "Silver perch" = c("CI", "NEPaP", "SA")
  )
  
  beta_long <- data.frame()
  
  for (param in param_panels) {
    if (param == "beta_RS" && !is.null(post$beta_RS)) {
      n_status <- dim(post$beta_RS)[3]
      n_iter   <- dim(post$beta_RS)[1]
      for (s in seq_along(spps)) {
        for (tt in seq_len(n_status)) {
          df <- data.frame(
            draw = seq_len(n_iter),
            value = post$beta_RS[, s, tt],
            species = spps[s],
            index = "",
            parameter = "beta_RS"
          )
          beta_long <- bind_rows(beta_long, df)
        }
      }
    } else if (param == "beta_site" && !is.null(post$beta_site)) {
      n_sites <- dim(post$beta_site)[3]
      n_iter  <- dim(post$beta_site)[1]
      for (s in seq_along(spps)) {
        # Loop through present site names and get their index in site_labels
        for (site_name in present_sites[[spps[s]]]) {
          i <- which(site_labels == site_name)
          if (length(i) == 1 && i <= n_sites) {
            df <- data.frame(
              draw = seq_len(n_iter),
              value = post$beta_site[, s, i],
              species = spps[s],
              index = site_labels[i],
              parameter = "beta_site"
            )
            beta_long <- bind_rows(beta_long, df)
          }
        }
      }
    } else if (param == "beta_TL" && !is.null(post$beta_TL)) {
      n_iter <- dim(post$beta_TL)[1]
      for (s in seq_along(spps)) {
        df <- data.frame(
          draw = seq_len(n_iter),
          value = post$beta_TL[, s],
          species = spps[s],
          index = "",
          parameter = "beta_TL"
        )
        beta_long <- bind_rows(beta_long, df)
      }
    } else if (param == "beta_FW" && !is.null(post$beta_FW)) {
      n_iter <- dim(post$beta_FW)[1]
      for (s in seq_along(spps)) {
        df <- data.frame(
          draw = seq_len(n_iter),
          value = post$beta_FW[, s],
          species = spps[s],
          index = "",
          parameter = "beta_FW"
        )
        beta_long <- bind_rows(beta_long, df)
      }
    } else if (param == "beta_site_RS" && !is.null(post$beta_site_RS)) {
      n_sites  <- dim(post$beta_site_RS)[2]
      n_status <- dim(post$beta_site_RS)[3]
      n_iter   <- dim(post$beta_site_RS)[1]
      for (s in seq_along(spps)) {
        for (site_name in present_sites[[spps[s]]]) {
          i <- which(site_labels == site_name)
          if (length(i) == 1 && i <= n_sites) {
            for (tt in seq_len(n_status)) {
              df <- data.frame(
                draw = seq_len(n_iter),
                value = post$beta_site_RS[, s, i, tt],
                species = spps[s],
                index = paste(site_labels[i], status_labels[tt], sep = " × "),
                parameter = "beta_site_RS"
              )
              beta_long <- bind_rows(beta_long, df)
            }
          }
        }
      }
    }
  }
  
  if (nrow(beta_long) == 0) stop("No valid parameters found in 'post' object for selected param_panels.")
  
  # Summarize posterior means and intervals
  beta_summary <- beta_long %>%
    group_by(parameter, species, index) %>%
    summarise(
      mean = mean(value),
      q025 = quantile(value, 0.2),
      q975 = quantile(value, 0.8),
      .groups = "drop"
    )
  
  pretty_labels <- c(
    "beta_RS" = "beta[RS]",
    "beta_site" = "beta[site]",
    "beta_TL" = "beta[TL]",
    "beta_FW" = "beta[FW]",
    "beta_site_RS" = "beta[site*','*RS]"
  )
  beta_summary$parameter <- factor(beta_summary$parameter, levels = names(pretty_labels), labels = pretty_labels)
  
  p <- ggplot(beta_summary, aes(x = mean, y = index, color = species)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_errorbarh(aes(xmin = q025, xmax = q975), height = 0, position = position_dodge(width = 0.5)) +
    facet_wrap(~ parameter, scales = "free", ncol = ncol, labeller = label_parsed) +
    scale_color_manual(values = species_colors, name = "Species") +
    labs(
      y = NULL,
      x = "Posterior Mean (95% Interval)",
      title = title
    ) +
    custom_theme() +
    theme(
      legend.title = element_blank(),
      legend.position = "top",
      strip.text = element_text(face = "bold", size = 14),
      axis.title.y = element_text(angle = 0, vjust = 0.5)
    )
  return(p)
}
