library(dplyr)
library(tidyr)
library(stringr)
library(ggmcmc)

## Pinfish ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_pinfish, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "SAV", "Tanaidacea"))
status_names = sort(c("LS", "CT"))  # Adjust if statuses are ordered differently

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = status_names[status_idx]
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.30),
    upper = quantile(value, 0.70),
    mean  = mean(value),
    .groups = "drop"
  )

# ---- Test interval overlap for each prey ----
# For each prey, compare intervals for all pairs of statuses.
overlap_results = ci_summary %>%
  group_by(Prey) %>%
  summarise(
    status1 = Status[1],
    status2 = Status[2],
    lower1 = lower[1], upper1 = upper[1],
    lower2 = lower[2], upper2 = upper[2],
    mean1 = mean[1], mean2 = mean[2],
    overlap = !(upper1 < lower2 | upper2 < lower1), # TRUE if intervals overlap
    overlap_pct = (min(upper1, upper2) - max(lower1, lower2)) /
      (max(upper1, upper2) - min(lower1, lower2)),
    .groups = "drop"
  ) %>%
  #round numbers
  mutate(across(where(is.numeric), ~ round(., 2)))

overlap_results = overlap_results %>%
  mutate(
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  select(Prey, lower1, upper1, mean1, lower2, upper2, mean2, overlap, overlap_pct)
  
colnames(overlap_results) <- c("Prey","LowerCT", "UpperCT", "MeanCT",
                               "LowerLS", "UpperLS", "MeanLS", "Overlap", "Overlap_%")

res_dir = file.path(here::here(), "res", "figures", "supplementary", "tables")
write.csv(overlap_results, file.path(res_dir, "pinfish_PO_overlap.csv"), row.names = FALSE)

## Croaker ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_micund, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))
status_names = sort(c("CT", "LS"))  # Adjust if statuses are ordered differently

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = status_names[status_idx]
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.10),
    upper = quantile(value, 0.90),
    mean  = mean(value),
    .groups = "drop"
  )

# ---- Test interval overlap for each prey ----
# For each prey, compare intervals for all pairs of statuses.
overlap_results = ci_summary %>%
  group_by(Prey) %>%
  summarise(
    status1 = Status[1],
    status2 = Status[2],
    lower1 = lower[1], upper1 = upper[1],
    lower2 = lower[2], upper2 = upper[2],
    mean1 = mean[1], mean2 = mean[2],
    overlap = !(upper1 < lower2 | upper2 < lower1), # TRUE if intervals overlap
    overlap_pct = (min(upper1, upper2) - max(lower1, lower2)) /
      (max(upper1, upper2) - min(lower1, lower2)),
    .groups = "drop"
  ) %>%
  #round numbers
  mutate(across(where(is.numeric), ~ round(., 2)))

overlap_results = overlap_results %>%
  mutate(
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  select(Prey, lower1, upper1, mean1, lower2, upper2, mean2, overlap, overlap_pct)

colnames(overlap_results) <- c("Prey", "LowerCT", "UpperCT", "MeanCT",
                               "LowerLS", "UpperLS", "MeanLS", "Overlap", "Overlap_%")

res_dir = file.path(here::here(), "res", "figures", "supplementary", "tables")
write.csv(overlap_results, file.path(res_dir, "croaker_PO_overlap.csv"), row.names = FALSE)

## Silver perch ##

# ---- Extract posterior draws for p_status ----
df_pstatus = ggmcmc::ggs(fit_silverperch, family = "p_status")

prey_names = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))
status_names = sort(c("LS", "CT"))  # Adjust if statuses are ordered differently

# Parse indices
df_parsed = df_pstatus %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = status_names[status_idx]
  ) %>%
  select(-indices, -prey_idx, -status_idx)

# ---- Compute 80% credible intervals ----
ci_summary = df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    lower = quantile(value, 0.20),
    upper = quantile(value, 0.80),
    mean  = mean(value),
    .groups = "drop"
  )

# ---- Test interval overlap for each prey ----
# For each prey, compare intervals for all pairs of statuses.
overlap_results = ci_summary %>%
  group_by(Prey) %>%
  summarise(
    status1 = Status[1],
    status2 = Status[2],
    lower1 = lower[1], upper1 = upper[1],
    lower2 = lower[2], upper2 = upper[2],
    mean1 = mean[1], mean2 = mean[2],
    overlap = !(upper1 < lower2 | upper2 < lower1), # TRUE if intervals overlap
    overlap_pct = (min(upper1, upper2) - max(lower1, lower2)) /
      (max(upper1, upper2) - min(lower1, lower2)),
    .groups = "drop"
  ) %>%
  #round numbers
  mutate(across(where(is.numeric), ~ round(., 2)))

overlap_results = overlap_results %>%
  mutate(
    overlap_pct = ifelse(overlap, overlap_pct, 0)
  ) %>%
  select(Prey, lower1, upper1, mean1, lower2, upper2, mean2, overlap, overlap_pct)

colnames(overlap_results) <- c("Prey", "LowerCT", "UpperCT", "MeanCT", "Overlap",
                               "LowerLS", "UpperLS", "MeanLS", "Overlap_%")

res_dir = file.path(here::here(), "res", "figures", "supplementary", "tables")
write.csv(overlap_results, file.path(res_dir, "silverperch_PO_overlap.csv"), row.names = FALSE)

