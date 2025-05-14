library(lubridate)
library(tidyverse)
library(neonUtilities)
library(patchwork)
library(sp)
library(raster)
library(tidybayes)
library(modelr)
library(loo)

# Analysis R Script 
library(brms)
library(ggridges)
library(ggdist)
library(rlang)
library(stringr)

library(posterior)   # as_draws_df()
library(tidyr)
library(stringr)
library(ggridges)
library(forcats)

library(kableExtra)
library(gridExtra)  # for grid.table()

#setwd("~/Desktop/Thesis_25")

df_flr_final_summary <- read_rds("Data/df_flr_final_summary.rds")
df_flr_final_filtered <- read_rds("Data/df_flr_final_filtered.rds")
dim(df_flr_final_filtered)


## Prepping data for model ---- 

data <- df_flr_final_filtered %>% dplyr::select(latitude, longitude, species, spring_temp,
                                               spring_precip, elevation, doy, life_history)
data <- na.omit(data)
unique(data$species)

#scaling 
doy_num <- as.numeric(data$doy)
doy_center <- mean(doy_num)
doy_scale <- sd(doy_num)
data$doy_sc <- (doy_num - doy_center) / doy_scale

latitude_num <- as.numeric(data$latitude)
latitude_center <- mean(latitude_num)
latitude_scale <- sd(latitude_num)
data$latitude_sc <- (latitude_num - latitude_center) / latitude_scale

stemp_num <- as.numeric(data$spring_temp)
stemp_center <- mean(stemp_num)
stemp_scale <- sd(stemp_num)
data$stemp_sc <- (stemp_num - stemp_center) / stemp_scale


sprecip_num <- as.numeric(data$spring_precip)
sprecip_center <- mean(sprecip_num)
sprecip_scale <- sd(sprecip_num)
data$sprecip_sc <- (sprecip_num - sprecip_center) / sprecip_scale     

elevation_num <- as.numeric(data$elevation)
elevation_center <- mean(elevation_num)
elevation_scale <- sd(elevation_num)
data$elevation_sc <- (elevation_num - elevation_center) / elevation_scale

# Define a few alternative models that essentialluy pose different hypotheses

# "No species-specific responses model"
# fit0:  no species-specific responses; just varying-intercept
formula0 <- doy_sc ~ 1 + stemp_sc * latitude_sc + stemp_sc * elevation_sc +
  pprecip_sc +   (1 | species)

# "No interactions model"
formula1 <- doy_sc ~ 1 + ptemp_sc + latitude_sc + elevation_sc + pprecip_sc + life_history + 
  (1 + latitude_sc + ptemp_sc + elevation_sc + pprecip_sc | species)

# FULL MODEL
formula_full <- doy_sc ~ 1 + latitude_sc * stemp_sc * elevation_sc * sprecip_sc * life_history +
  (1 + latitude_sc * stemp_sc * elevation_sc * sprecip_sc | species)

# only 2-way interactions
formula_two_way <- doy_sc ~ 1 + 
  # Main effects
  latitude_sc + stemp_sc + elevation_sc + sprecip_sc + life_history +
  # Two-way interactions
  latitude_sc:stemp_sc + 
  latitude_sc:elevation_sc + 
  latitude_sc:sprecip_sc +
  latitude_sc:life_history +
  stemp_sc:elevation_sc + 
  stemp_sc:sprecip_sc +
  stemp_sc:life_history +
  elevation_sc:sprecip_sc +
  elevation_sc:life_history +
  sprecip_sc:life_history +
  # Random effects - only main effects and two-way interactions
  (1 + 
     # Main effects
     latitude_sc + stemp_sc + elevation_sc + sprecip_sc + 
     # Two-way interactions
     latitude_sc:stemp_sc + 
     latitude_sc:elevation_sc + 
     latitude_sc:sprecip_sc +
     stemp_sc:elevation_sc + 
     stemp_sc:sprecip_sc +
     elevation_sc:sprecip_sc | 
     species)

# Run FULL model ----
fit <- brm(
  formula = formula_full,
  data = data,
  family = gaussian(),  
  #control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  #init = "0"
)
saveRDS(fit, file = "Data/fit_full.RDS")

# Run 2-way model ----
fit_two_way <- brm(
  formula = formula_two_way,
  data = data,
  family = gaussian(),  
  #control = list(adapt_delta = 0.99, max_treedepth = 15),
  save_pars = save_pars(all = TRUE),
  chains = 4,
  iter = 4000,
  warmup = 2000,
  thin = 3,
  cores = 4,
  #init = "0"
)

saveRDS(fit_two_way, file = "Data/fit_two_way.RDS")





# Read in fitted models ----
fit_two_way <- readRDS("Data/fit_two_way.RDS") 
fit_full <- readRDS("Data/fit_full.RDS") 

### model selection ----
loo1 <- loo(fit_full)
loo1
loo2 <- loo(fit_two_way)
loo_compare(loo1, loo2)


# A) Plot FULL MODEL -----
fit <- fit_full
fit<- fit_fullselect1
summary(fit)

original_data <- fit$data 
colnames(original_data)
spp <- unique(original_data$species)
spp
life_history <- unique(original_data$life_history)
life_history



## ridge plots of posterior distributions ----
get_variables(fit)

# 1) Define which terms get species‐specific random slopes vs overall only
species_main   <- c("latitude_sc","stemp_sc","elevation_sc","sprecip_sc")
species_int    <- c(
  "latitude_sc:stemp_sc", 
  "latitude_sc:elevation_sc", 
  "latitude_sc:sprecip_sc",
  "stemp_sc:elevation_sc", 
  "stemp_sc:sprecip_sc", 
  "elevation_sc:sprecip_sc"
)
species_terms  <- c(species_main, species_int)

# overall also includes life_historyperennial
overall_terms  <- c(species_main, "life_historyperennial", species_int)

# 2) grab the full posterior as a tibble
draws_df <- as_draws_df(fit)

# Sample 1% of the dataset for testing
draws_df_sample <- draws_df[sample(nrow(draws_df), size = 0.01 * nrow(draws_df)), ]

# 3) species‐specific totals
fixed_long <- draws_df_sample %>%
  tidyr::pivot_longer(
    cols      = starts_with("b_"),
    names_to  = "param",
    values_to = "fixed"
  ) %>%
  mutate(param = str_remove(param, "^b_")) %>%
  filter(param %in% species_terms)


random_long <- draws_df_sample %>%
  tidyr::pivot_longer(
    cols      = starts_with("r_species"),
    names_to  = "raw",
    values_to = "random"
  ) %>%
  mutate(
    species = str_extract(raw, "(?<=\\[)[^,]+"),
    param   = str_extract(raw, "(?<=,)[^\\]]+")
  ) %>%
  dplyr::select(-raw) %>%
  filter(param %in% species_terms)

species_draws <- dplyr::inner_join(
  fixed_long, random_long,
  by = c(".chain", ".iteration", ".draw", "param")
) %>%
  transmute(
    .draw, species, param,
    total = fixed + random
  )

# 4) overall fixed‐effect draws (including life_historyperennial)
overall_draws <- draws_df %>%
  tidyr::pivot_longer(
    cols      = starts_with("b_"),
    names_to  = "param",
    values_to = "estimate"
  ) %>%
  mutate(param = str_remove(param, "^b_")) %>%
  filter(param %in% overall_terms)

# 5) binning helper
make_bins <- function(df, est_col, group_cols) {
  df %>%
    group_by(across(all_of(group_cols))) %>%
    summarize(
      mean_eff = mean(.data[[est_col]]),
      p_gt0    = mean(.data[[est_col]] > 0),
      p_lt0    = mean(.data[[est_col]] < 0),
      .groups  = "drop"
    ) %>%
    mutate(
      direction = if_else(mean_eff > 0, "Positive", "Negative"),
      prob      = pmax(p_gt0, p_lt0),
      cert_bin  = case_when(
        prob > 0.95 ~ ">95%",
        prob > 0.90 ~ "90-95%",
        TRUE        ~ "<90%"
      ),
      effect_cat = factor(
        paste0(direction, " (", cert_bin, ")"),
        levels = c(
          "Positive (>95%)",
          "Positive (90-95%)",
          "Negative (>95%)",
          "Negative (90-95%)",
          "Positive (<90%)",
          "Negative (<90%)"
        )
      )
    )
}

species_bins  <- make_bins(species_draws,  "total", c("species","param"))
overall_bins  <- make_bins(overall_draws, "estimate", "param")

# 6) join bins back on
species_draws2 <- species_draws %>%
  dplyr::left_join(
    dplyr::select(species_bins, species, param, mean_eff, effect_cat),
    by = c("species","param")
  )

overall_draws2 <- overall_draws %>%
  dplyr::left_join(
    dplyr::select(overall_bins, param, mean_eff, effect_cat),
    by = "param"
  )

# 7) relabel and define panel ordering
labmap <- c(
  stemp_sc              = "Temperature",
  sprecip_sc            = "Precipitation",
  latitude_sc           = "Latitude",
  elevation_sc          = "Elevation",
  life_historyperennial = "Life history"
)
relabel <- function(t){
  if (str_detect(t,":")){
    parts <- str_split(t,":",simplify=TRUE)
    return(str_c(labmap[parts], collapse=" × "))
  }
  labmap[[t]] %||% t
}

# factor levels: mains first (with temp & precip at very top), then interactions starting with temp × precip
main_labs <- c("Temperature","Precipitation","Latitude","Elevation","Life history")
int_labs  <- c(
  "Temperature × Precipitation",
  "Temperature × Elevation",
  "Latitude × Temperature",
  "Latitude × Elevation",
  "Latitude × Precipitation",
  "Elevation × Precipitation"
)
term_levels <- c(main_labs, int_labs)
term_levels

species_draws2 <- species_draws2 %>%
  mutate(
    term_lab = map_chr(param, relabel),
    term_lab = factor(term_lab, levels = term_levels)
  ) %>%
  dplyr::select(.draw, species, term_lab, total, mean_eff, effect_cat)

overall_draws2 <- overall_draws2 %>%
  mutate(
    term_lab = map_chr(param, relabel),
    term_lab = factor(term_lab, levels = term_levels)
  ) %>%
  dplyr::select(.draw, term_lab, estimate = estimate, mean_eff, effect_cat)


saveRDS(species_draws2, file="Data/species_draws2_woelevation.rds") #for fit_fullselect1
saveRDS(species_draws2, file="Data/species_draws2.rds") #for fit_full 



# 1 - Overall parameters plot ----
# Define shared palette
shared_pal <- c(
  "Positive (>95%)"   = "#00008B",
  "Positive (90–95%)" = "#4169E1",
  "Negative (>95%)"   = "#8B0000",
  "Negative (90–95%)" = "#CD5C5C",
  "Positive (<90%)"   = "#A9A9A9",
  "Negative (<90%)"   = "#A9A9A9"
)

# Climate panel (no legend)
p_climate <- overall_draws2 %>% 
  filter(term_lab != "Life history") %>%
  ggplot(aes(x = estimate, y = term_lab, fill = effect_cat)) +
  geom_density_ridges(
    scale           = 1.2,
    rel_min_height  = 0.01,
    alpha           = 0.8,
    quantile_lines  = TRUE,
    quantiles       = c(0.025, 0.5, 0.975)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  labs(title = "Climate terms", x = "Estimate", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y     = element_text(face = "bold"),
    legend.position = "none"
  )
p_climate

# Life‑history panel (no legend, squashed vertically)
p_life <- overall_draws2 %>% 
  filter(term_lab == "Life history") %>%
  ggplot(aes(x = estimate, y = 1, fill = effect_cat)) +
  geom_density_ridges(
    scale           = 1.2,
    rel_min_height  = 0.01,
    alpha           = 0.8,
    quantile_lines  = TRUE,
    quantiles       = c(0.025, 0.5, 0.975)
  ) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_y_continuous(breaks = NULL) +
  labs(title = "Life history", x = "Estimate", y = NULL) +
  theme_minimal() +
  theme(
    axis.text.y     = element_blank(),
    legend.position = "none",
    plot.margin     = margin(t = 5, r = 5, b = 5, l = 5),
    aspect.ratio    = 1.8     # squashes it to half the height
  )
p_life


# Stitch together and collect one shared legend, on the right
combined <- (p_climate | p_life) +
  plot_layout(
    widths = c(3, 1),
    guides = "collect"
  ) &
  scale_fill_manual(
    values = shared_pal,
    name   = "Direction & certainty"
  ) &
  theme(
    legend.position  = "right",
    legend.key.size  = unit(0.6, "lines"),
    legend.text      = element_text(size = 8),
    legend.title     = element_text(size = 10)
  )
combined
ggsave(plot=combined,"Analysis_Images/full_model/overall_params.pdf", width=7, height=3.5)

# 2 - Species‐specific ridge plot----
p_species <- ggplot(species_draws2,
                    aes(x=total, y=fct_reorder(species, mean_eff), fill=effect_cat)) +
  geom_density_ridges(
    scale          = 0.9,
    rel_min_height = 0.01,
    alpha          = 0.7,
    quantile_lines = TRUE,
    quantiles      = c(0.025, 0.5, 0.975)
  ) +
  geom_vline(xintercept = 0, linetype="dashed") +
  facet_wrap(~ term_lab, scales="free_x", ncol=4) +
  scale_fill_manual(values = c(
    "Positive (>95%)"  = "#00008B",
    "Positive (90–95%)"= "#4169E1",
    "Negative (>95%)"  = "#8B0000",
    "Negative (90–95%)"= "#CD5C5C",
    "Positive (<90%)"  = "#A9A9A9",
    "Negative (<90%)"  = "#A9A9A9"
  )) +
  labs(
    title    = "Species‑specific slopes",
    subtitle = "Main effects & climate×climate interactions",
    x        = "Estimate (fixed + random)",
    y        = "Species",
    fill     = "Direction & certainty"
  ) +
  theme_minimal() +
  theme(
    legend.position = "bottom",
    strip.text      = element_text(face="bold"),
    axis.text.y     = element_text(size=6)
  )
p_species
ggsave(plot=p_species,"Analysis_Images/full_model/species_params.pdf", width=8, height=10)

# Parameter Tables ----

head(species_draws2)
head(overall_draws2)

# --- helper to compute summary stats & bins for any draw column ---
summarize_terms <- function(df, value, group_vars) {
  df %>%
    group_by(across(all_of(group_vars))) %>%
    summarise(
      mean_est   = mean({{value}}),
      median_est = median({{value}}),
      lower_95   = quantile({{value}}, 0.025),
      upper_95   = quantile({{value}}, 0.975),
      p_gt0      = mean({{value}} > 0),
      p_lt0      = mean({{value}} < 0),
      .groups    = "drop"
    ) %>%
    mutate(
      direction = if_else(mean_est > 0, "Positive", "Negative"),
      prob      = pmax(p_gt0, p_lt0),
      cert_bin  = case_when(
        prob > 0.95             ~ ">95%",
        prob > 0.90 & prob <= 0.95 ~ "90–95%",
        TRUE                    ~ "<90%"
      ),
      effect_cat = paste0(direction, " (", cert_bin, ")")
    ) %>%
    dplyr::select(-p_gt0, -p_lt0, -prob, -direction, -cert_bin)
}

# --- 1. Overall parameter summaries ---
overall_summary <- summarize_terms(
  overall_draws2,
  value = estimate,
  group_vars = "term_lab"
) 

# --- 2. Species‐specific summaries ---
species_summary <- summarize_terms(
  species_draws2,
  value = total,
  group_vars = c("species","term_lab")
)

# --- 3. Display with kable ---
kable(
  overall_summary,
  caption = "Overall posterior summaries by parameter",
  digits  = 3
) %>%
  kable_styling(full_width = FALSE)

kable(
  species_summary,
  caption = "Species‐specific posterior summaries",
  digits  = 3
) %>%
  kable_styling(full_width = FALSE)


# 1) Save to CSV
write.csv(overall_summary, "overall_summary.csv", row.names = FALSE)
write.csv(species_summary, "species_summary.csv", row.names = FALSE)

# 2) Save overall_summary to PDF
library(gridExtra)

# 1a) Make a grob with tiny text
tbl_grob <- tableGrob(
  overall_summary,
  rows     = NULL,
  theme    = ttheme_minimal(
    core   = list(fg_params = list(cex = 0.6)),
    colhead = list(fg_params = list(cex = 0.7, fontface = "bold"))
  )
)
pdf("overall_summary.pdf", width = 8.5, height = 11)  # letter size
grid.draw(tbl_grob)
dev.off()


# 3) Save species_summary to PDF
#    (you may need to increase height if you have many rows)
n_per_page <- 25
chunks     <- split(
  species_summary,
  (seq_len(nrow(species_summary))-1) %/% n_per_page
)

pdf("species_summary_paged.pdf", width = 8.5, height = 11)
for(i in seq_along(chunks)) {
  grid.newpage()
  grid.table(
    chunks[[i]],
    rows = NULL,
    theme = ttheme_minimal(
      core   = list(fg_params = list(cex = 0.7)),
      colhead = list(fg_params = list(cex = 0.8, fontface = "bold"))
    )
  )
}
dev.off()


