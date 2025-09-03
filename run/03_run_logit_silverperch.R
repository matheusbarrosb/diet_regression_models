# Packages ---------------------------------------------------------------------
pack_list = c("dplyr", "rstan", "tidyr", "ggplot2", "purrr", "readr",
              "stringr", "here", "reshape2", "PNWColors", "ggpubr", "ggmcmc")
if (!all(pack_list %in% rownames(installed.packages()))) {
  install.packages(pack_list[!pack_list %in% rownames(installed.packages())])
}
lapply(pack_list, library, character.only = TRUE)

# Source functions -------------------------------------------------------------
fun_dir = here::here("R")
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

# Load data --------------------------------------------------------------------
data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))

# Data wrangling ---------------------------------------------------------------
spps  = sort(c("BAICHR"))
sites = sort(c("CI", "NEPaP", "SA"))
prey  = sort(c("Amphipod", "Crustacean", "Fish", "Isopod", "Polychaete", "Tanaidacea"))
RS_levels = sort(c("CT", "LS"))

wg_data =
  raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  filter(Group %in% prey) %>%
  mutate(present = 1) %>%
  pivot_wider(names_from = Group, values_from = present, values_fill = list(present = 0))

# Gather data for stan ---------------------------------------------------------
prey_mat = 
  wg_data %>%
  select(Species.code, Site, Treatment, Amphipod, Crustacean, Fish, Isopod, Polychaete, Tanaidacea) %>%
  mutate(Species.code = as.numeric(factor(Species.code, levels = spps))) %>%
  mutate(Site = as.numeric(factor(Site, levels = sites))) %>%
  mutate(Treatment = as.numeric(factor(Treatment, levels = RS_levels))) %>%
  select(-Species.code, -Site, -Treatment) %>%
  as.matrix()

site   = as.numeric(as.factor(sort(wg_data$Site)))
status = as.numeric(as.factor(sort(wg_data$Treatment)))
spp    = as.numeric(as.factor(sort(wg_data$Species.code)))
TL     = wg_data$Length; TL[is.na(TL)] = mean(TL, na.rm = TRUE)

stan_data = list(
  N = nrow(prey_mat),
  S = length(unique(sites)),
  K = length(unique(RS_levels)),
  G = length(unique(prey)),
  site     = site,
  status   = status,
  TL       = TL,
  prey_mat = prey_mat,
  meanTLs_site = tapply(TL, site, mean),
  meanTLs_status = tapply(TL, status, mean)
)

# fit the model ----------------------------------------------------------------
model_dir = here::here("stan/")
stanc(paste0(model_dir, "logit_int.stan"))
rstan_options(auto_write = TRUE)

fit_silverperch = stan(
  file   = paste0(model_dir, "logit_int.stan"),  
  data   = stan_data,
  chains = 3,        
  iter   = 10000,      
  warmup = 1000,       
  cores  = 3,        
  seed   = 444       
)

### 1. Probabilities of occurrence ###
post = rstan::extract(fit_silverperch, pars = "p_status")$p_status

df =
  ggmcmc::ggs(fit_silverperch, family = "p_status")

prey_names = prey
status_names = c("CT", "LS")

df_parsed = df %>%
  mutate(
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 2)[,1]),
    status_idx = as.integer(str_split_fixed(indices, ",", 2)[,2]),
    Prey   = prey_names[prey_idx],
    Status = status_names[status_idx]
  ) %>%
  select(-indices, -prey_idx, -status_idx)

head(df_parsed)

baichr_probs = 
  df_parsed %>%
  group_by(Prey, Status) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.2),
    upper = quantile(value, 0.8)
  ) %>%
  
  ggplot(aes(x = Prey, y = mean, color = Status)) +
  geom_point(position = position_dodge(width = 0.5)) +
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0,
                position = position_dodge(width = 0.5)) +
  custom_theme() +
  coord_flip() +
  ylab("Probability of encounter") +
  xlab("") +
  theme(legend.title = element_blank(),
        legend.position = "none",
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank(),
        axis.title.y = element_blank()) +
  scale_color_manual(values = pnw_palette("Bay", 2)) 

## 2. beta_tl plot ##
## beta_TL is indexed by prey only
post = rstan::extract(fit_silverperch, pars = "beta_TL")$beta_TL
df =
  ggmcmc::ggs(fit_silverperch, family = "beta_TL") %>%
  mutate(
    # Extract indices from Parameter: beta_TL[prey]
    indices = str_extract(Parameter, "\\[(.*?)\\]"),
    indices = str_remove_all(indices, "\\[|\\]"),
    prey_idx   = as.integer(str_split_fixed(indices, ",", 1)[,1]),
    Prey   = prey_names[prey_idx]
  ) %>%
  select(-indices, -prey_idx) %>%
  
  group_by(Prey) %>%
  summarise(
    mean = mean(value),
    lower = quantile(value, 0.2),
    upper = quantile(value, 0.8)
  )

baichr_beta_TL =
  df %>%
  ggplot(aes(x = Prey, y = mean)) +
  geom_point() +
  geom_errorbar(aes(ymin = lower, ymax = upper), width = 0) +
  custom_theme() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
  coord_flip() +
  xlab("") +
  ylab(expression(beta[TL]))

baichr_logit_plot = 
ggarrange(
  baichr_beta_TL,
  baichr_probs,
  ncol = 2,
  widths = c(1, 0.8)
)

#### Arrange for final plot ####
final_plot = 
  ggarrange(
    pinfish_logit_plot,
    micund_logit_plot,
    baichr_logit_plot,
    nrow = 3, common.legend = TRUE,
    labels = c("A", "B", "C")
  )

ggsave(
  filename = here::here("res", "figures", "logit_res.pdf"),
  plot = final_plot,
  width = 5, height = 6, dpi = 300
)

# Posterior predictive check ---------------------------------------------------
obs_df <- as.data.frame(prey_mat)
obs_df$site <- site

obs_freq <- obs_df %>%
  group_by(site) %>%
  summarise(across(everything(), mean)) %>%
  pivot_longer(-site, names_to = "Prey", values_to = "Observed") %>%
  mutate(Site = sites[site]) %>%
  select(Site, Prey, Observed)

post_pred <- rstan::extract(fit_micund, pars = "p_site")$p_site 

post_pred_df = reshape2::melt(post_pred)
colnames(post_pred_df) <- c("Iteration", "Prey_idx", "Site_idx", "Modelled")
post_pred_df$Prey <- prey[post_pred_df$Prey_idx]
post_pred_df$Site <- sites[post_pred_df$Site_idx]

ppc_summary <- post_pred_df %>%
  group_by(Prey, Site) %>%
  summarise(
    Model_Mean  = mean(Modelled),
    Model_Lower = quantile(Modelled, 0.025),
    Model_Upper = quantile(Modelled, 0.975),
    .groups = "drop"
  )

silverperch_ppcheck =
  obs_freq %>%
  left_join(ppc_summary, by = c("Prey", "Site")) %>%
  
  ggplot(aes(x = Observed, y = Model_Mean, color = Prey)) +
  geom_point(size = 1) +
  geom_errorbar(aes(ymin = Model_Lower, ymax = Model_Upper), width = 0.0) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey60") +
  labs(
    x = "",
    y = ""
  ) +
  theme(legend.title = element_blank()) +
  custom_theme() +
  xlim(0,1) + ylim(0,1); print(silverperch_ppcheck)

#### Arrange ppcheck plots #####

ggarrange(
  pinfish_ppcheck,
  croaker_ppcheck,
  silverperch_ppcheck,
  nrow = 1, common.legend = TRUE
)

ggsave(
  filename = here::here("res", "figures", "ppcheck.pdf"),
  plot = last_plot(),
  width = 8, height = 3, dpi = 300
)




