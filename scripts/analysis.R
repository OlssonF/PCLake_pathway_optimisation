#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Example problem scenario comparisons - read in, plot, and summarise the 
## Date: 2026-01-23
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

library(tidyverse)
library(ggpubr)
library(GGally)
library(ggh4x)


# Read in data ------------------------------------------------------------
setwd(here::here())
possible_measures <- read_csv('possible_measures.csv')


## Rdata output of the optimisation
rds_files <-  list.files('output', '.RData', full.names = T)

for (input_file in rds_files) {
  assign(gsub('.RData', '', basename(input_file)),
         read_rds(input_file))
}

## Best member of each iteration
bestmemit <- lapply(list.files('output', 'bestmemit', full.names = TRUE),
                    read_csv, show_col_types = F)
names(bestmemit) <- gsub('.csv', '', 
                         basename(list.files('output', 'bestmemit', full.names = TRUE)))

## Last iteration population, obj_function output
lastpop <- lapply(list.files('output', 'lastpop_', full.names = TRUE),
                  read_csv, show_col_types = F)
names(lastpop) <- gsub('.csv', '', 
                       basename(list.files('output', 'lastpop_', full.names = TRUE)))


## Last iteration population, obj_function output
lastpopstate <- lapply(list.files('output', 'lastpopstate_', full.names = TRUE),
                       read_delim, show_col_types = F)
names(lastpopstate) <- gsub('.csv', '', 
                            basename(list.files('output', 'lastpopstate_', full.names = TRUE)))


## Last iteration population, obj_function output
lastpoppathways <- lapply(list.files('output', 'lastpoppathways_', full.names = TRUE),
                          read_csv, show_col_types = F)

names(lastpoppathways) <- gsub('.csv', '', 
                               basename(list.files('output', 'lastpoppathways_', full.names = TRUE)))


# Aesthetics --------------------------------------------------------------
labels_measures <- c(expression(atop('P load', (gP~m^-2~d^-1))),
                     expression(atop('P load', (gP~m^-2~d^-1))),
                     # expression(paste("P load\n(mg ", L^-1, ")")),
                     # expression(paste("P load\n(mg ", L^-1, ")")), 
                     expression(atop("Marsh area", "(fraction)")),
                     expression(atop("Vegetation removed", "(fraction)")),
                     expression(atop("Day of vegetation", "removal (day of year)")),
                     "cDredInterval", "cDredStart",
                     expression(atop("Start P load", "(year)")),
                     expression(atop("Start P load", "(year)")),
                     expression(atop("Start marsh", " area (year)")),
                     expression(atop("Start vegetation", "removal (year)")))



names(labels_measures) <- c(possible_measures$parameter)
cols_measures <- c("#3E4A89FF",
                   "#3E4A89FF",
                   "#6DCD59FF",
                   "#F89441FF",
                   "grey",
                   "grey", "grey",
                   "grey",
                   "grey",
                   "grey",
                   "grey")


labels_states <- c(oChlaEpi =expression(atop('Chlorophyll-a concentration', (mu*g~L^-1))),
                   aDSubVeg =expression(atop('Submerged vegetation biomass', (gDW~m^-2))),
                   aDFish   =expression(atop("Benthivorous fish biomass", (gDW~m^-2))))

labels_states_str <- c(oChlaEpi = "atop('Chlorophyll-a concentration', (mu*g~L^-1))",
                       aDSubVeg = "atop('Submerged vegetation biomass', (gDW~m^-2))",
                       aDFish   = "atop('Benthivorous fish biomass', (gDW~m^-2))")
#
# PS1 - Single ----------------------------------------------------------------
summary_single$iter
nrow(lastpop$lastpop_single)

lastpop$lastpop_single |> 
  filter(fn_out <= 0) |> 
  pivot_longer(mPLoadEpi:fMarsh_lag, names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

# Plot all measures as scatter 
single_scatter <- lastpop$lastpop_single |> 
  filter(fn_out <= 0) |> 
  ggplot(aes(x=mPLoadEpi_lag, y = mPLoadEpi, 
             colour = fMarsh, size = fMarsh_lag)) + 
  geom_point(alpha = 0.8) +
  scale_colour_viridis_c(option  ='viridis', begin = 0.2, end = 1) +
  theme_bw() +
  labs(dictionary = labels_measures) 

# plot all measures as a parallel coordinate plot
single_parallel <- lastpop$lastpop_single |> 
  filter(fn_out <= 0) |> 
  ggparcoord(scale = 'std', # rescaling
             columns = 1:4, # measures
             showPoints = T, alphaLines = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = labels_measures, name = 'Measure') +
  scale_y_continuous(name = 'Normalised value') + 
  theme(axis.text.x = element_text(vjust = -0.5, hjust = 0.5))

ggarrange(single_scatter, single_parallel, widths = c(1,0.7), align = 'h',
          labels = c("A)", "B)"), label.x = -0.01) |> 
  ggsave(filename = './output/plots/figure_2.png', width = 24, height = 12, unit = 'cm')

# PS2 - Two horizon --------------------------------------------
summary_twohorizon2$desired_states
summary_twohorizon2$iter
nrow(lastpop$lastpop_twohorizon2)
lastpop$lastpop_twohorizon2 |> 
  filter(fn_out <= 0) 

lastpop$lastpop_twohorizon2 |> 
  filter(fn_out <= 0) |> 
  mutate(mPLoadEpi_lag_first = ifelse(mPLoadEpi_lag < mPLoadEpi_lag2, mPLoadEpi_lag, mPLoadEpi_lag2),
         mPLoadEpi_lag_second = ifelse(mPLoadEpi_lag > mPLoadEpi_lag2, mPLoadEpi_lag, mPLoadEpi_lag2),
         mPLoadEpi_first = ifelse(mPLoadEpi_lag < mPLoadEpi_lag2, mPLoadEpi, mPLoadEpi2),
         mPLoadEpi_second = ifelse(mPLoadEpi_lag > mPLoadEpi_lag2, mPLoadEpi, mPLoadEpi2)) #|> 
# pivot_longer(cols = c("fMarsh", "fMarsh_lag", "mPLoadEpi_lag_first", "mPLoadEpi_lag_second", "mPLoadEpi_first", "mPLoadEpi_second"),
#              names_to = 'measure') |> 
# reframe(.by = 'measure', 
#         n = n(),
#         mean = mean(value),
#         sd = sd(value),
#         min = min(value),
#         max = max(value)) |> 
# mutate(cv = sd/mean)

# Plot successful pathways
# do.call(rbind.data.frame, summary_twohorizon2$desired_states[1]) #  pivot_longer(cols = everything(), names_pattern = "^([^.]+).(.*)$", names_to = c('var', 'thing'))
twohorizon2_ds <- data.frame(opt_var = unique(lastpopstate$lastpopstate_twohorizon2$opt_var),
                             lower_range = 0,
                             upper_range = c(100,20))


# Which pathways were successful?
successful_pathways <- lastpopstate$lastpopstate_twohorizon2 |> 
  #pivot_wider(id_cols = ID, names_from = variable, values_from = output) |>
  # pivot_longer(cols = paste(names(desired_states1), c(5,30),sep = '_'),
  #              names_to = 'opt_var', values_to = 'out') |> 
  full_join(twohorizon2_ds,
            by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> # check within range
  filter(n() == nrow(twohorizon2_ds),  
         .by = ID) |> 
  distinct(ID) |> 
  pull(ID)

pathways_2h <- lastpoppathways$lastpoppathways_twohorizon2 |> 
  filter(ID %in% successful_pathways) |>
  mutate(.by = year, ID = row_number()) |> # renumber the pathways
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0)) |>
  # this whole convoluted, nested ifelse() generates the timeseries based on the lags of mPLoadEpi
  mutate(mPLoadEpi_use = ifelse(year < mPLoadEpi_lag &
                                  year < mPLoadEpi_lag2, 0.01,
                                ifelse(year < mPLoadEpi_lag &
                                         year > mPLoadEpi_lag2, mPLoadEpi2,
                                       ifelse(year > mPLoadEpi_lag &
                                                year < mPLoadEpi_lag2, mPLoadEpi,
                                              ifelse(year > mPLoadEpi_lag &
                                                       year > mPLoadEpi_lag2 &
                                                       mPLoadEpi_lag > mPLoadEpi_lag2, mPLoadEpi, mPLoadEpi2))))) |>
  mutate(fMarsh = fMarsh_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  
  select(all_of(c('ID', 'year', 'oChlaEpi', 'mPLoadEpi', 'fMarsh'))) |>
  pivot_longer(cols = !any_of(c('ID','year'))) |>
  # filter(ID %in% 7:9) |>
  ggplot(aes(x=year, y = value,
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y',  nrow = 9,
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'),
                           strip = strip_nested(text_y = list(element_text(),
                                                              element_text(colour = 'white',
                                                                           size = 1)),
                                                background_y = list(element_rect(),
                                                                    element_blank()),
                                                by_layer_y = TRUE)) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.title = element_text(hjust = 0.5),
        panel.spacing.y = unit(c( rep( c( rep(0.3,2), 0.7), 1), rep(0.3,2)), "lines")) +
  facetted_pos_scales(y = list(name == "fMarsh" ~ scale_y_continuous(limits = c(0,1),
                                                                     n.breaks = 2),
                               name == "mPLoadEpi" ~ scale_y_continuous(limits = c(0,0.01),
                                                                        n.breaks = 2),
                               name == "oChlaEpi" ~ scale_y_continuous(limits = c(0,200),
                                                                       n.breaks = 2,
                                                                       minor_breaks = c(20,100)))) +
  scale_colour_manual(values = c(cols_measures, 'grey20'), 
                      name = 'Lake state', 
                      breaks = c(names(labels_measures),'oChlaEpi'), 
                      labels = c(labels_measures, labels_states)) +
  guides(colour = guide_legend(direction = 'vertical', ncol = 2))

ggsave(pathways_2h, filename = './output/plots/figure_3.png', 
       width = 10, height = 15, unit = 'cm')

# PS3 - Multi --------------------------------------------
summary_multi$desired_states
summary_multi$iter
nrow(lastpop$lastpop_multi)

lastpop$lastpop_multi |> 
  filter(fn_out <= 0) |> 
  distinct(runID) |>  nrow()

lastpop$lastpop_multi |> 
  filter(fn_out <= 0) |> 
  pivot_longer(mPLoadEpi:fMarsh_lag, names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

# Plot all measures as scatter 
multi_scatter <- lastpop$lastpop_multi |> 
  filter(fn_out <= 0) |> 
  ggplot(aes(x=mPLoadEpi_lag, y = mPLoadEpi, 
             colour = fMarsh, size = fMarsh_lag)) + 
  geom_point(alpha = 0.8) +
  scale_colour_viridis_c(option  ='viridis', begin = 0.2, end = 1) +
  theme_bw()+
  labs(dictionary = labels_measures)  +
  scale_y_continuous(breaks = seq(0.0001, 0.0006, 0.0001), labels = scales::label_number(accuracy = 0.0001)) 

# plot all measures as a parallel coordinate plot
multi_parallel <- lastpop$lastpop_multi |> 
  filter(fn_out <= 0) |> 
  ggparcoord(scale = 'std', # rescaling
             columns = 1:4, # measures
             showPoints = T, alphaLines = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = labels_measures, name = 'Measure') +
  scale_y_continuous(name = 'Normalised value') + 
  theme(axis.text.x = element_text(vjust = -0.5, hjust = 0.5))

ggarrange(multi_scatter, multi_parallel, widths = c(1,0.8), align = 'h',
          labels = c("A)", "B)"), label.x = -0.01) |> 
  ggsave(filename = './output/plots/figure_S3.png', width = 24, height = 12, unit = 'cm')

# PS4 - MultiES -----------------------------------------
summary_multiES$iter
nrow(lastpop$lastpop_multiES)

lastpop$lastpop_multiES |> 
  #   # filter(fn_out <= 0) |> 
  pivot_longer(mPLoadEpi:fManVeg_lag, names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

# objectively "best pathway"
lastpop$lastpop_multiES |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_multiES) # what were the associated states

best_multiES <- lastpop$lastpop_multiES |> 
  slice_min(fn_out) |> # lowest objective function
  select(runID) |> 
  left_join(lastpoppathways$lastpoppathways_multiES, by = join_by(runID == ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.01))  |> 
  
  mutate(fManVeg = fManVeg_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  
  select(all_of(c('ID', 'year', 'aDSubVeg', 'aDFish', 'oChlaEpi', 'mPLoadEpi', 'fManVeg'))) |> #'aDFish'
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  ggplot(aes(x=year, y = value, 
             # size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(name), scales = 'free_y',  nrow = 16,
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black')) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top', 
        strip.background = element_blank(),
        strip.text = element_blank()
        # panel.spacing.y = unit(c( rep( c( rep(0.2,3), 0.6), 3), rep(0.2,3)),"lines")
  ) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, title = '')) +
  facetted_pos_scales(y = list(name == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                      n.breaks = 2),
                               name == "mPLoadEpi" ~ scale_y_continuous(limits = c(0,0.01),
                                                                        n.breaks = 2),
                               name == "oChlaEpi" ~ scale_y_continuous(limits = c(0,200),
                                                                       n.breaks = 2),
                               name == "aDSubVeg" ~ scale_y_continuous(limits = c(0,100),
                                                                       n.breaks = 2),
                               name == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                     n.breaks = 2))) +
  scale_colour_manual(values = c('grey20', 'seagreen','gold', cols_measures), 
                      name = 'Lake state',
                      breaks = c(names(labels_states), names(labels_measures)), 
                      labels = c(labels_states, labels_measures))
ggsave(best_multiES, filename = 'output/plots/figure_S4.png',
       height = 18, width = 18, units = 'cm')

# out of all the pathways how many achieve each of the indicator targets
summary_multiES$desired_states

multiES_ds <- data.frame(opt_var = names(summary_multiES$desired_states),
                         lower_range = sapply(summary_multiES$desired_states,
                                              function(x) min(x$target)),
                         upper_range = sapply(summary_multiES$desired_states,
                                              function(x) max(x$target))) |> 
  mutate(opt_var_val = row_number())

lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |>
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = sum(opt_var_val)) |> # 1 = 1, 3 = 1+2, 4 = 1+3, 5 = 2+3, 6 = 1+2+3
  reframe(.by = c(total_achieve, ind_ID),
          n = n())|> 
  mutate(prop = n/nrow(lastpop$lastpop_multiES))

ggarrange(lastpopstate$lastpopstate_multiES |>
            full_join(multiES_ds, by = join_by(opt_var)) |>
            ggplot(aes(y = out, x = fManVeg, colour = fManVeg_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'plasma', begin = 0.9, end  = 0.5) +
            facet_wrap(~opt_var, 
                       scales   = "free_y",
                       labeller = labeller(opt_var = as_labeller(labels_states_str, 
                                                                 label_parsed))) +
            theme_bw() +
            scale_x_continuous(name = "Vegetation removed (fraction)") +
            labs(dictionary = labels_measures),
          
          lastpopstate$lastpopstate_multiES |> 
            full_join(multiES_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = mPLoadEpi, colour = mPLoadEpi_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_x_continuous(n.breaks = 4, labels = scales::label_comma()) +
            scale_colour_viridis_c(option = 'magma', begin = 0.1, end = 0.6) +
            facet_wrap(~opt_var, 
                       scales   = "free_y",
                       labeller = labeller(opt_var = as_labeller(labels_states_str, 
                                                                 label_parsed))) +
            theme(legend.position = 'right') +
            theme_bw() +
            labs(x = expression(P~load~(gP~m^-2~d^-1))) + 
            labs(dictionary = labels_measures), 
          
          nrow = 2, align = 'hv') |> 
  ggsave(filename = './output/plots/figure_4.png', width = 20, height = 14, unit = 'cm')

multiES_parallel <- lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES, by = join_by(ID)) |> 
  
  ggparcoord(scale = 'std', # rescaling
             columns = 4:7, # measures
             groupColumn = 3, order = 'anyClass',
             showPoints = T, alphaLines = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = labels_measures, name = 'Measure') +
  scale_y_continuous(name = 'Normalised value') + 
  theme(axis.text.x = element_text(vjust = -0.5, hjust = 0.5)) +
  scale_color_manual(labels = c('WQ', 'WQ + RF1', "WQ + RF2"),
                     values = c("#414487FF", "#1E9C89FF", "#BBDF27FF"),
                     name = "Indicators achieved")

ggsave(multiES_parallel, filename = 'output/plots/figure_5.png', 
       height = 10, width = 15, units ='cm')

# what is the general pattern of the three groups?
lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES, by = join_by(ID)) |> 
  reframe(.by = ind_ID,
          across(mPLoadEpi:fManVeg_lag, .fns = list(mean = mean,
                                                    sd = sd,
                                                    min = min,
                                                    max = max),
                 .names = "{.col}__{.fn}")) |> 
  pivot_longer(cols = -ind_ID,
               names_sep = '__',
               names_to = c('var', 'stat')) |> 
  pivot_wider(names_from = stat,
              values_from = value) |> 
  arrange(var)



# PS5 - Prioritised -------------------------------------
summary_prioritisation$iter
nrow(lastpop$lastpop_prioritisation)

lastpop$lastpop_prioritisation |> 
  # filter(fn_out <= 0) |> 
  pivot_longer(mPLoadEpi:fManVeg_lag, names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

# objectively "best pathway"
lastpop$lastpop_prioritisation |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_prioritisation) # what were the associated states


best_prior <- lastpop$lastpop_prioritisation |> 
  slice_min(fn_out) |> # lowest objective function
  select(runID) |> 
  left_join(lastpoppathways$lastpoppathways_multiES, by = join_by(runID == ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.01))  |> 
  
  mutate(fManVeg = fManVeg_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  
  select(all_of(c('ID', 'year', 'aDSubVeg', 'aDFish', 'oChlaEpi', 'mPLoadEpi', 'fManVeg'))) |> #'aDFish'
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  ggplot(aes(x=year, y = value, 
             # size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(name), scales = 'free_y',  nrow = 16,
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black')) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top', 
        strip.background = element_blank(),
        strip.text = element_blank()
        # panel.spacing.y = unit(c( rep( c( rep(0.2,3), 0.6), 3), rep(0.2,3)),"lines")
  ) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, title = '')) +
  facetted_pos_scales(y = list(name == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                      n.breaks = 2),
                               name == "mPLoadEpi" ~ scale_y_continuous(limits = c(0,0.01),
                                                                        n.breaks = 2),
                               name == "oChlaEpi" ~ scale_y_continuous(limits = c(0,200),
                                                                       n.breaks = 2),
                               name == "aDSubVeg" ~ scale_y_continuous(limits = c(0,100),
                                                                       n.breaks = 2),
                               name == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                     n.breaks = 2))) +
  scale_colour_manual(values = c('grey20', 'seagreen','gold', cols_measures), 
                      name = 'Lake state',
                      breaks = c(names(labels_states), names(labels_measures)), 
                      labels = c(labels_states, labels_measures))


ggsave(best_prior, filename = 'output/plots/figure_6.png', height = 18, width = 18, unit = 'cm')

# which targets were met?
prioritisation_ds <- data.frame(opt_var = names(summary_prioritisation$desired_states),
                                lower_range = sapply(summary_prioritisation$desired_states,
                                                     function(x) min(x$target)),
                                upper_range = sapply(summary_prioritisation$desired_states,
                                                     function(x) max(x$target))) |> 
  mutate(opt_var_val = row_number())

lastpopstate$lastpopstate_prioritisation |> 
  full_join(prioritisation_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = sum(opt_var_val)) |> # 1 = 1, 3 = 1+2, 4 = 1+3, 5 = 2+3, 6 = 1+2+3
  reframe(.by = c(total_achieve, ind_ID),
          n = n()) |> 
  mutate(prop = n/nrow(lastpop$lastpop_prioritisation))

ggarrange(lastpopstate$lastpopstate_prioritisation |>
            full_join(multiES_ds, by = join_by(opt_var)) |>
            ggplot(aes(y = out, x = fManVeg, colour = fManVeg_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'plasma', begin = 0.9, end  = 0.5) +
            facet_wrap(~opt_var, 
                       scales   = "free_y",
                       labeller = labeller(opt_var = as_labeller(labels_states_str, 
                                                                 label_parsed))) +
            theme_bw() +
            scale_x_continuous(name = "Vegetation removed (fraction)") +
            labs(dictionary = labels_measures),
          
          lastpopstate$lastpopstate_prioritisation |> 
            full_join(multiES_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = mPLoadEpi, colour = mPLoadEpi_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_x_continuous(n.breaks = 4, labels = scales::label_comma()) +
            scale_colour_viridis_c(option = 'magma', begin = 0.1, end = 0.6) +
            facet_wrap(~opt_var, 
                       scales   = "free_y",
                       labeller = labeller(opt_var = as_labeller(labels_states_str, 
                                                                 label_parsed))) +
            theme(legend.position = 'right') +
            theme_bw() +
            labs(x = expression(P~load~(mg ~ L^-1))) +
            labs(dictionary = labels_measures), 
          
          nrow = 2, align = 'hv') |> 
  ggsave(filename = './output/plots/figure_S6.png', width = 20, height = 14, unit = 'cm')


lastpopstate$lastpopstate_prioritisation |> 
  full_join(prioritisation_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES, by = join_by(ID)) |> 
  
  ggparcoord(scale = 'std', # rescaling
             columns = 4:7, # measures
             groupColumn = 3, order = 'anyClass',
             showPoints = T, alphaLines = 0.5) +
  theme_bw() +
  # scale_x_discrete(labels = labels_measures, name = 'Measure')  +
  scale_y_continuous(name = 'Normalised value') + 
  theme(axis.text.x = element_text(vjust = -0.5, hjust = 0.5)) +
  scale_color_manual(labels = c('WQ', 'WQ + RF1'),
                     values = c("#414487FF", "#1E9C89FF"), name = "Indicators achieved")

# what is the general pattern of the two groups?
lastpopstate$lastpopstate_prioritisation |> 
  full_join(prioritisation_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(row_number()))) |> 
  full_join(lastpopstate$lastpopstate_multiES, by = join_by(ID)) |> 
  reframe(.by = ind_ID,
          across(mPLoadEpi:fManVeg_lag, .fns = list(mean = mean,
                                                    sd = sd,
                                                    min = min,
                                                    max = max),
                 .names = "{.col}__{.fn}")) |> 
  pivot_longer(cols = -ind_ID,
               names_sep = '__',
               names_to = c('var', 'stat')) |> 
  pivot_wider(names_from = stat,
              values_from = value) |> arrange(var)

# Supplementary information ------------------------
# plot all measures as a parallel coordinate plot
single_parallel_SI <- lastpop$lastpop_single |> 
  filter(fn_out <= 0) |> 
  ggparcoord(scale = 'globalminmax', # rescaling
             columns = 1:4, # measures
             showPoints = T, alphaLines = 0.5) +
  theme_bw() +
  scale_x_discrete(labels = labels_measures, name = 'Measure') +
  scale_y_continuous(name = 'Non-scaled value')
ggsave(single_parallel_SI, filename = 'output/plots/figure_S1.png', 
       height = 12, width = 12, units = 'cm')

# PS4 - the one that got the fish but not the macrophytes
fish_multiES <- lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range),
         opt_var == 'aDFish') |> 
  select(ID) |> 
  left_join(lastpoppathways$lastpoppathways_multiES, by = join_by(ID)) |> 
  mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.01))  |> 
  
  mutate(fManVeg = fManVeg_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  
  select(all_of(c('ID', 'year', 'aDSubVeg', 'aDFish', 'oChlaEpi', 'mPLoadEpi', 'fManVeg'))) |> #'aDFish'
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  ggplot(aes(x=year, y = value, 
             # size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(name), scales = 'free_y',  nrow = 16,
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black')) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top', 
        strip.background = element_blank(),
        strip.text = element_blank()
        # panel.spacing.y = unit(c( rep( c( rep(0.2,3), 0.6), 3), rep(0.2,3)),"lines")
  ) +
  guides(color=guide_legend(nrow=2, byrow=TRUE, title = '')) +
  facetted_pos_scales(y = list(name == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                      n.breaks = 2),
                               name == "mPLoadEpi" ~ scale_y_continuous(limits = c(0,0.01),
                                                                        n.breaks = 2),
                               name == "oChlaEpi" ~ scale_y_continuous(limits = c(0,200),
                                                                       n.breaks = 2),
                               name == "aDSubVeg" ~ scale_y_continuous(limits = c(0,100),
                                                                       n.breaks = 2),
                               name == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                     n.breaks = 2))) +
  scale_colour_manual(values = c('grey20', 'seagreen','gold', cols_measures), 
                      name = 'Lake state',
                      breaks = c(names(labels_states), names(labels_measures)), 
                      labels = c(labels_states, labels_measures))
ggsave(fish_multiES, filename = 'output/plots/figure_S5.png',
       height = 18, width = 18, units = 'cm')
