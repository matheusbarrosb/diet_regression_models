# Packages ---------------------------------------------------------------------
pack_list = c("dplyr", "rstan", "tidyr", "ggplot2", "purrr", "readr",
              "stringr", "here", "reshape2", "PNWColors", "ggpubr")
if (!all(pack_list %in% rownames(installed.packages()))) {
  install.packages(pack_list[!pack_list %in% rownames(installed.packages())])
}
lapply(pack_list, library, character.only = TRUE)

# Source functions -------------------------------------------------------------
fun_dir = here::here("R")
# source multiple functions in the R/ folder
fun_files = list.files(fun_dir, pattern = "\\.R$", full.names = TRUE)
lapply(fun_files, source)

# Load data --------------------------------------------------------------------
data_dir = here::here("data")
raw_data = read.csv(file.path(data_dir, "rawData.csv"))

# Data wrangling ---------------------------------------------------------------
# Species and site definitions as before
spps  = sort(c("LAGRHO", "MICUND", "BAICHR"))
sites = sort(c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"))
RS_levels = sort(c("CT", "LS")) # or use unique(wg_data$Treatment)

wg_data = raw_data %>%
  filter(Species.code %in% spps) %>%
  filter(Site %in% sites) %>%
  filter(Treatment %in% RS_levels) %>%
  group_by(Fish_ID_year)


spp  = as.numeric(factor(wg_data$Species.code, levels = spps))
site = as.numeric(factor(wg_data$Site, levels = sites))
RS   = as.numeric(factor(wg_data$Treatment, levels = RS_levels))
GW   = as.numeric(wg_data$Gut.weight)
FW   = as.numeric(wg_data$Wet.Weight)
TL   = as.numeric(wg_data$Length)

# put everything into a dataframe and exclude NA values
df =
data.frame(
  spp  = spp,
  site = site,
  RS   = RS,
  GW   = GW,
  FW   = FW,
  TL   = TL
) %>%
  filter(!is.na(spp) & !is.na(site) & !is.na(RS) & !is.na(GW) & !is.na(FW) & !is.na(TL))

# put into list for stan
stan_data = list(
  N       = nrow(df),
  S       = length(spps),         
  I       = length(sites),        
  TT      = length(RS_levels),    
  species = df$spp,               
  site    = df$site,             
  status  = df$RS,                
  GW      = df$GW,
  FW      = df$FW,
  TL      = df$TL
)

# fit --------------------------------------------------------------------------
model_dir = here::here("stan/")
stanc(paste0(model_dir, "GF_centered.stan"))
rstan_options(auto_write = TRUE)

fit <- stan(
  file   = paste0(model_dir, "GF_centered.stan"),  
  data   = stan_data,
  chains = 2,        
  iter   = 5000,      
  warmup = 1000,       
  cores  = 2,        
  seed   = 444       
)

# plotting ---------------------------------------------------------------------

### PROBABILITIES OF LS > CT ###
species_names <- c("Pinfish", "Croaker", "Silver perch")
status_names <- c("CT", "LS") # CT is reference (status==1), LS is non-reference

post_diffs_plot =
plot_post_diffs(
  beta_RS_array = post$beta_RS,
  species_names = species_names,
  status_names = status_names,
  title = "A",
  prob_text_size = 3,
  strip_text_size = 10
) +
  xlim(-1,2)

### PREDICTIVE CHECK ###
post <- rstan::extract(fit, pars = c("alpha", "beta_site", "beta_RS", "beta_TL", "beta_FW"))
ppcheck_plot =
plot_gutweight_ppcheck(post, stan_data, title = "B", palette = "Bay", strip_text_size = 10)

### PLOT COEFFICIENTS ###
coef_plot =
plot_betas(
  post = post,
  spps = c("Pinfish", "Croaker", "Silver perch"),
  site_labels = c("AM", "DR", "HWP", "LB", "NEPaP", "SA", "CI"),
  status_labels = c("CT", "LS"),
  param_panels = c("beta_site", "beta_RS", "beta_TL", "beta_FW", "beta_site_RS"),
  title = "C",
  palette = "Bay",
  ncol = 2
)

# arrange 
plot1 =
ggarrange(
  post_diffs_plot,
  ppcheck_plot,
  nrow = 2,
  label.y = "Density"
); print(plot1)

plot2 = 
ggarrange(
  plot1,
  coef_plot,
  ncol = 2
  )

fig_dir = here::here("res", "figures")

ggsave(
  filename = file.path(fig_dir, "gut_fullness.pdf"),
  plot = plot2,
  width = 10,
  height = 6,
  units = "in",
  device = cairo_pdf
)
