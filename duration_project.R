library(haven)
library(dplyr)
library(tidyr)
library(survival)
library(survminer)
library(survival)
library(tidyverse)

udat <- read_dta("~/Documents/S2/duartions_models/project_exam/fi.dta")
view(udat)

sum(is.na(udat))

udat <- udat[,1:134] ## get rid of some superfluous variables - there was 166 and we select 134
udat <- as_tibble(udat)

dim(udat)  ## N = 225,821

glimpse(udat[,1:36])

udat$type

# -11 = PDB + RR / -10 = PDB / -1 = RR / 0 = control 

## Computation of average spells when durations are truncated at 104 weeks
udat %>%
  mutate(dur104 = dur,
         dur104 = ifelse(dur104 > 104, 104, dur104)) ->
  udat


udat <- udat %>%
  mutate(
    type_name = factor(type, 
                       levels = c(0, -1, -10, -11), 
                       labels = c("Control", "RR", "PBD", "PBD+RR"))
  )



# Method 1: Algebraic Aggregation (Exact Table 4 Replication)
didd_in_diff <- udat %>%
  group_by(type, after) %>%
  summarize(Y_bar = mean(dur104, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = after, values_from = Y_bar, names_prefix = "Period_") %>%
  rename(Before = Period_0, After = Period_1) %>%
  mutate(Difference = After - Before)

control_diff <- didd_in_diff %>% 
  filter(type == 0) %>% 
  pull(Difference)

didd_in_diff <- didd_in_diff %>%
  mutate(DiD_Estimator = Difference - control_diff)

print(didd_in_diff)



library(survival)
library(survminer)
library(dplyr)

# Define a production-grade plotting function to ensure standardized axes and scales
generate_km_plot <- function(dataset, target_type, title_text) {
  
  # Isolate the specific treatment or control group
  strata_data <- dataset %>% filter(type == target_type)
  
  # Compute the Kaplan-Meier estimator
  # Surv(dur, uncc) rigorously processes right-censored data
  km_fit <- survfit(Surv(dur, uncc) ~ after, data = strata_data)
  
  # Render the survival curve
  plot <- ggsurvplot(
    km_fit,
    data = strata_data,
    conf.int = FALSE,           # Display 95% confidence intervals to assess statistical overlap
    pval = TRUE,               # Execute Log-Rank test to confirm global statistical difference
    risk.table = FALSE,
    title = title_text,
    legend.title = "Reform Status",
    legend.labs = c("Pre-Reform (after=0)", "Post-Reform (after=1)"),
    xlab = "Unemployment Duration (Weeks)",
    ylab = "Survival Probability S(t)",
    xlim = c(0, 104),          # Force x-axis standardization across all plots for comparability
    palette = c("grey", "black")
  )
  
  return(plot)
}

# 1. Generate the plot for the Potential Benefit Duration extension group (PBD)
plot_pbd <- generate_km_plot(udat, -10, "Kaplan-Meier: ePBD Group (Duration Extension)")

# 2. Generate the plot for the Replacement Rate increase group (RR)
plot_rr <- generate_km_plot(udat, -1, "Kaplan-Meier: eRR Group (Financial Increase)")

# 3. Generate the plot for the Control group
plot_control <- generate_km_plot(udat, 0, "Kaplan-Meier: Control Group")

# 4 PBD + RR group 
plot_pbd_rr <- generate_km_plot(udat, -11, "Kaplan-Meier: ePBD + eRR Group")


# Output the plots
print(plot_pbd)
print(plot_rr)
print(plot_control)
print(plot_pbd_rr)





library(survival)
library(KernSmooth)
library(dplyr)
library(ggplot2)

# 1. Define the Hazard Extraction and Smoothing Engine
generate_smooth_hazard <- function(dataset, target_type, bandwidth_val = 2.5) {
  
  # Initialize storage for pre- and post-reform data
  smoothed_data <- data.frame()
  
  for (period in c(0, 1)) {
    # Isolate strata
    strata_data <- dataset %>% filter(type == target_type, after == period)
    
    # Fit base Kaplan-Meier
    km_fit <- survfit(Surv(dur, uncc) ~ 1, data = strata_data)
    
    # Extract step-function survival probabilities strictly at integer weeks
    grid_weeks <- 1:104
    km_summary <- summary(km_fit, times = grid_weeks, extend = TRUE)
    
    S_t <- km_summary$surv
    t_grid <- km_summary$time
    
    # Compute discrete density and empirical hazard
    # f(t) = S(t-1) - S(t); lambda(t) = f(t) / S(t-1)
    S_t_minus_1 <- c(1, S_t[-length(S_t)]) 
    f_t <- S_t_minus_1 - S_t
    raw_hazard <- f_t / S_t_minus_1
    
    # Apply Local Polynomial Regression Smoothing (locpoly)
    # degree = 1 (Local Linear) is preferred over Nadaraya-Watson (degree = 0) to reduce boundary bias
    smooth_fit <- locpoly(x = t_grid, 
                          y = raw_hazard, 
                          bandwidth = bandwidth_val, 
                          degree = 1, 
                          gridsize = 104)
    
    # Store bounded results (hazard cannot be negative)
    temp_df <- data.frame(
      time = smooth_fit$x,
      hazard = pmax(smooth_fit$y, 0),
      Reform = factor(period, levels = c(0, 1), labels = c("Pre-Reform (after=0)", "Post-Reform (after=1)"))
    )
    
    smoothed_data <- bind_rows(smoothed_data, temp_df)
  }
  
  return(smoothed_data)
}

# 2. Execute Smoothing for Policy Groups
hazard_pbd <- generate_smooth_hazard(udat, -10, bandwidth_val = 3)
hazard_rr  <- generate_smooth_hazard(udat, -1, bandwidth_val = 3)

# 3. Construct the Production-Grade Visualizations
plot_smoothed_hazard <- function(data_input, title_text) {
  ggplot(data_input, aes(x = time, y = hazard, color = Reform, linetype = Reform)) +
    geom_line(linewidth = 1.2) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 0.12)) + # Focus on the critical 60-week window
    scale_x_continuous(breaks = seq(0, 190, by = 20)) +
    scale_color_manual(values = c("#1f77b4", "#d62728")) +
    labs(title = title_text,
         x = "Unemployment Duration (Weeks)",
         y = "Smoothed Hazard Rate \u03bb(t)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(face = "bold"))
}

plot_haz_pbd <- plot_smoothed_hazard(hazard_pbd, "Smoothed Exit Hazard: ePBD Group")
plot_haz_rr  <- plot_smoothed_hazard(hazard_rr, "Smoothed Exit Hazard: eRR Group")

print(plot_haz_pbd)
print(plot_haz_rr)

udat %>%
  mutate(all = tr * (t39 + t52) ) ->
  udat

breaks <- seq(from=3,to=59, by=4)
labels <- paste("(", c(0,breaks), ",", c(breaks,104), "]",sep="")

gux <- survSplit(Surv(dur104,uncc) ~., data=udat, cut = breaks,
                 end = "time", event="death", start="start", episode="interval")

gux %>%
  mutate(exposure = time - start,
         interval=factor(interval+1, labels = labels) ) ->
  gux


# 1. Estimation du modèle complet avec les paramètres causaux
# On inclut les effets fixes de groupe, l'effet macro (after), 
# et les 5 interventions via les interactions temporelles (DiD)
pwe <- glm(death ~ interval + 
             t39 + t52 + tr + all +                  # Effets fixes de groupe
             after +                                 # Paramètre 1 : Effet macro-économique
             after:t39 + after:t52 + after:tr +      # Paramètres 2, 3, 4 : Effets isolés
             after:t39:tr + after:t52:tr +           # Paramètres 5, 6 : Effets synergiques
             offset(log(exposure)), 
           family = poisson(link = "log"), 
           data = gux)

# 2. Affichage strict des effets de traitement (Les paramètres Delta)
library(stargazer)

stargazer(pwe, 
          dep.var.caption="", dep.var.labels="",
          omit = "interval", # On masque les 15 intervalles pour voir les paramètres causaux
          omit.table.layout = "n", star.cutoffs = NA,
          keep.stat=c("n", "ll"), no.space=TRUE,
          header=FALSE,
          title="The Causal Treatment Effects (Delta Parameters)", type="text")
