# CSIA Code V 0.1 August 28, 2023 Robin B. Trayler
# EVERYTHING is Case Sensitive!

# load required libraries -----------------------------------------------------
library(tidyverse)
library(broom)
library(viridis)
library(cowplot)
theme_set(theme_minimal())

# source required functions ---------------------------------------------------
source('./R/assign_peaks.R')
source('./R/parse_GC_csv.R')

# file path -------------------------------------------------------------------
file_name = '~/Box Sync/Data Repository/GCC/GCC20221012.csv'

# set up required tables ------------------------------------------------------
# define retention times 
retention_time <- tribble(
  ~retention_time, ~peak_name,
        236,        "CO2_ref",
        274,        "CO2_ref",
        314,        "CO2_ref",
        623,        "CO2_ref",
        750,          "Ala",
        809,          "Val",
        862,          "Gly",
        887,          "Ile",
        898,          "Leu",
        945,          "Nor",
        1050,         "Pro",
        1160,       "CO2_ref",
        1238,         "Asp",
        1260,         "Thr",
        1376,         "Glu",
        1421,         "Phe",
        1718,     "CO2_ref",
        1757,     "CO2_ref",
        1797,     "CO2_ref")

# define standard isotope compositions
standards <- tribble(~sample,    ~peak_name,  ~d13C_true,
                     "L-AA Mix",      "Ala",     -18.43,
                     "L-AA Mix",      "Val",     -11.88,
                     "L-AA Mix",      "Gly",     -42.79,
                     "L-AA Mix",      "Ile",     -12.17,
                     "L-AA Mix",      "Leu",     -28.33,
                     "L-AA Mix",      "Nor",     -27.76,
                     "L-AA Mix",      "Pro",     -10.48,
                     "L-AA Mix",      "Asp",     -22.76,
                     "L-AA Mix",      "Thr",     -10.86,
                     "L-AA Mix",      "Phe",     -13.11,
                     "L-AA Mix",      "Glu",      -11.3)

# get amino acid names ----------------
AA_names <- standards |> 
  pull(peak_name)

# define standards for corrections --------------------------------------------
standard_used      = c('L-AA Mix')
drift_standard     = c('L-AA Mix')
linearity_standard = c('L-AA Mix')

# What Kind of Corrections to do ----------------------------------------------
drift_correction     = TRUE
linearity_correction = FALSE

# read in the data ------------------------------------------------------------
data = read_csv(file = file_name, col_names = TRUE) |> 
  parse_GC_csv() |> 
  assign_peaks(retention_time = retention_time) |> 
  filter(!(peak_name %in% 'CO2_ref')) |> 
  # filter(!is.na(peak_name)) |>  # drop unassigned peaks
  full_join(standards, by = c('sample', 'peak_name')) |> 
  group_by(sample, peak_name) |> 
  mutate(d13C_dev = d13C_measured - mean(d13C_measured)) |> 
  ungroup()

# Check for problems ----------------------------------------------------------
# count the number of peaks in each group
peak_no <- data |> 
  group_by(id1) |> 
  summarize(first(sample), 
            n = n())

# check if any samples have too many or two few peaks
if(!any(peak_no$n == length(AA_names))) {
  warning('Warning! Some samples have a different number peaks than specified retention times.')
  # stop()
}

# check if any samples are missing Amino Acid assignments 
ids <- unique(data$id1)
for(i in seq_along(ids)) {
  AA <- data |> 
    filter(id1 == ids[i]) |> 
    pull(peak_name)
  if(sum(!(AA_names %in% AA)) > 0) {
    paste('WARNING! Sample with Identifer 1 of',
          ids[i], 'is missing the',
          AA_names[!(AA_names %in% AA)], 
          'peak') |> 
      warning()
    
  }
  # stop()
}

# Visualize Peak Assignments --------------------------------------------------
# visualize the peak assignments
data |> 
  # filter(peak_name %in% AA_names) |> 
  ggplot(mapping = aes(x = retention_time,
                       y = sample,
                       color = peak_name)) + 
  geom_point(position = position_dodge2(width = 0.5)) +
  geom_linerange(mapping = aes(xmin = start,
                               xmax = end),
                 position  = position_dodge2(width = 0.5)) +
  # facet_wrap(~peak_name,
  #            scales = 'free_x') +
  geom_vline(data = retention_time |>  
               filter(peak_name != 'CO2_ref'),
             mapping = aes(xintercept = retention_time),
             linewidth = 0.75,
             lineend = 'round',
             color = 'darkgrey',
             linetype = 'solid') + 
  scale_color_viridis(discrete = TRUE,
                      option = 'plasma',
                      na.value = 'black',
                      end = 0.8) +
  theme(legend.position = 'top',
        legend.title = element_blank()) + 
  xlab('retention time (s)') + 
  ylab('sample name') + 
  guides(color=guide_legend(ncol=11)) + 
  ggtitle('Labeled Peaks')

# Zoom in to check retention time alignments
data |> 
  filter(peak_name %in% AA_names) |>
  ggplot(mapping = aes(x = retention_time,
                       y = sample,
                       color = peak_name)) + 
  geom_point(position = position_dodge2(width = 0.5)) +
  geom_linerange(mapping = aes(xmin = start,
                               xmax = end),
                 position  = position_dodge2(width = 0.5)) +
  facet_wrap(~peak_name,
             scales = 'free_x') +
  geom_vline(data = retention_time |>  
               filter(peak_name != 'CO2_ref'),
             mapping = aes(xintercept = retention_time),
             linewidth = 0.75,
             lineend = 'round',
             color = 'darkgrey',
             linetype = 'solid') + 
  scale_color_viridis(discrete = TRUE,
                      option = 'plasma',
                      na.value = 'black',
                      end = 0.8) +
  theme(legend.position = 'top',
        legend.title = element_blank()) + 
  xlab('retention time (s)') + 
  ylab('sample name') + 
  guides(color=guide_legend(ncol=11)) + 
  ggtitle('Labeled Peaks')

# Remove unlabeled peaks ------------------------------------------------------
data <- data |> 
  filter(peak_name %in% AA_names)

# Drift Correction ------------------------------------------------------------
if(drift_correction) {
  # check for drift Amino acid by Amino acid
  pre_drift <- data |> 
    ggplot(mapping = aes(x = seq_nr, 
                         y = d13C_dev,
                         color = sample)) + 
    geom_point() + 
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~peak_name) + 
    theme(legend.position = 'top') + 
    ylab(expression(delta^13*C[deviation])) + 
    xlab('sequence number') + 
    scale_color_viridis(discrete = TRUE, 
                        option = 'plasma', 
                        end = 0.8) + 
    ggtitle('not drift corrected')
  
  # fit a linear regression to the standards for each AA
  drift <- data |>
    filter(sample %in% standard_used) |>
    nest(data = -c(peak_name)) |> 
    mutate(fit = map(data, ~lm(d13C_dev ~ seq_nr, data = .x)),
           tidied = map(fit, broom::tidy)) |> 
    unnest(tidied) |> 
    select(
      peak_name, 
      term, 
      estimate) |> 
    mutate(term = case_when(term == '(Intercept)' ~ 'intercept',
                            term == 'seq_nr' ~ 'slope')) 
  
  # rotate the data and isolate the slopes and intercept  
  drift <- drift |> 
    group_by(peak_name) |> 
    pivot_wider(names_from = term, 
                values_from = estimate) |> 
    select( 
      drift_intercept = intercept, 
      drift_slope     = slope) |> 
    ungroup()
  
  # add the drift slope and intercept to the data and apply it
  data <- drift |>
    full_join(data, by = 'peak_name') |> 
    mutate(d13C_measured = d13C_measured - 
             (seq_nr * drift_slope - drift_intercept)) |> 
    group_by(peak_name, sample) |> 
    mutate(d13C_dev = d13C_measured - mean(d13C_measured)) |> 
    ungroup()
  
  # post drift plot
  post_drift <- data |> 
    ggplot(mapping = aes(x = seq_nr, 
                         y = d13C_dev,
                         color = sample)) + 
    geom_point() + 
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~peak_name) + 
    scale_color_viridis(discrete = TRUE, 
                        option = 'plasma', 
                        end = 0.8) + 
    ylab(expression(delta^13*C[deviation])) + 
    xlab('sequence number') + 
    theme(legend.position = 'top') + 
    ggtitle('drift corrected')
  
  # plot it 
  plot_grid(pre_drift, post_drift)
} 

# Linearity correction --------------------------------------------------------
if(linearity_correction) {
  pre_lin <- data |>
    # filter(sample %in% standard_used) |>
    ggplot(mapping = aes(x = amp_44,
                         y = d13C_dev,
                         color = sample)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~peak_name, scales = 'free_x') + 
    scale_color_viridis(discrete = TRUE, 
                        option = 'plasma',
                        end = 0.8) + 
    theme(legend.position = 'top') + 
    ylab(expression(delta^13*C[deviation])) + 
    xlab('Mass 44 Amplitude') + 
    ggtitle('not linearity corrected')
  
  # calculate slopes and intercepts 
  linearity <- data |>
    filter(sample %in% standard_used) |>
    nest(data = -peak_name) |>
    mutate(fit = map(data, ~lm(d13C_dev ~ amp_44, 
                               data = .x)),
           tidied = map(fit, broom::tidy)) |>
    unnest(tidied) |>
    select(peak_name,
           term,
           estimate) |>
    mutate(term = case_when(term == '(Intercept)' ~ 'intercept',
                            term == 'amp_44' ~ 'slope'))
  
  # rotate the data and isolate the slopes and intercept  
  linearity <- linearity |>
    group_by(peak_name) |>
    pivot_wider(names_from = term,
                values_from = estimate) |>
    select(linearity_intercept = intercept,
           linearity_slope     = slope) |>
    ungroup()
  
  # add the drift slope and intercept to the data and apply it
  data <- linearity |> 
    full_join(data, by = 'peak_name') |>
    mutate(d13C_measured = d13C_measured - 
             (amp_44 * linearity_slope - linearity_intercept)) |>
    group_by(peak_name, 
             sample) |>
    mutate(d13C_dev = d13C_measured - mean(d13C_measured)) |>
    ungroup()
  
  # post linearity plot
  post_lin <- data |>
    ggplot(mapping = aes(x = amp_44,
                         y = d13C_dev,
                         color = sample)) +
    geom_point() +
    geom_smooth(method = 'lm', se = FALSE) + 
    facet_wrap(~peak_name, scales = 'free_x') + 
    scale_color_viridis(discrete = TRUE,
                        option = 'plasma', 
                        end = 0.8) + 
    theme(legend.position = 'top') + 
    ylab(expression(delta^13*C[deviation])) + 
    xlab('Mass 44 Amplitude') + 
    ggtitle('linearity corrected')
  
  plot_grid(pre_lin, post_lin)
}

# Standard Correction ---------------------------------------------------------
# define fractions. 
# These are from a spreadsheet provided by Mario Hernandez.
true_fraction <- tribble(~peak_name, ~fraction,
                         'Ala',        0.50,
                         'Val',        0.63,
                         'Gly',        0.40,
                         'Ile',        0.67,
                         'Leu',        0.67,
                         'Nor',        0.67,
                         'Pro',        0.63,
                         'Asp',        0.50,
                         'Thr',        0.50,
                         'Phe',        0.75,
                         'Glu',        0.55)

# apply the correction --------------------------------------------------------
data <- data |> 
  # add in the fractions
  full_join(true_fraction, by = 'peak_name') |> 
  # use the standards to determine the correction factor
  filter(sample %in% standard_used) |> 
  mutate(d13C_corr_fact = (d13C_measured - fraction * d13C_true) / (1 - fraction)) |> 
  # calculate the average factor by AA
  group_by(peak_name) |> 
  summarize(d13C_corr_fact = mean(d13C_corr_fact),
            fraction = unique(fraction)) |> 
  # add the full data back in
  full_join(data, by = c('peak_name')) |> 
  # apply the correction
  mutate(d13C_calculated = (d13C_measured - ((1 - fraction) * d13C_corr_fact)) / fraction)

# summarize the data ----------------------------------------------------------
final_data <- data |> 
  select(identifier_1 = id1,
         sample = sample,
         amino_acid = peak_name,
         amount = amount,
         amp_44 = amp_44,
         d13C_corrected = d13C_calculated)

# summary statistics ----------------------------------------------------------
summary_stats <- final_data |> 
  group_by(sample, amino_acid) |> 
  summarize(d13C_mean = mean(d13C_corrected), 
            d13C_sd   = sd(d13C_corrected), 
            n = n()) |> 
  mutate_if(is.numeric, round, 2)

# plot the data to look for anomalies -----------------------------------------
final_data |> 
  ggplot(mapping = aes(x = amino_acid,
                       y = d13C_corrected,
                       color = sample)) + 
  geom_jitter(size = 3,
              alpha = 0.75,
              width = 0.1,
              height = 0) +
  scale_color_viridis(discrete = TRUE, option = 'plasma', 
                      end = 0.8) +
  xlab('amino acid') + 
  ylab(expression(delta^13*C[corrected])) + 
  theme(legend.position = 'top')
