#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Example problem scenario comparisons - read in, plot, and summarise the 
## Date: 2026-01-23
## Author: 
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
labels_measures <- c(expression(paste("P load (mg ", L^-1, ")")), 
                     expression(paste("P load (mg ", L^-1, ")")),#"P load (mg L^-1)",
                     "Fraction\nas marsh area",
                     "Fraction of\nvegetation removed",
                     "Day of year\nfor vegetation removal",
                     "cDredInterval", "cDredStart",
                     "Timing of\nP load (year)",
                     "Timing of\nP load (year)",
                     "Timing of marsh\n area (year)",
                     "Timing of vegetation\n removal (year)")

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


labels_states <- c(oChlaEpi = "Chlorophyll a concnetration",
                   aDSubVeg = "Submerged vegetation biomass",
                   aDFish = "Benthivorous fish biomass")
# labels_states <- c(oChlaEpi =expression(paste("Chlorophyll a concnetration (", mu*g~L^-1, ")")),
#                    aDSubVeg =expression(paste("Submerged vegetation biomass (", gDw~m^-2, ")")),
#                    aDFish = expression(paste("Benthivorous fish biomass (", gDw~m^-2, ")")))
# labels_states <- list(expression(paste("Chlorophyll a concnetration (", mu*g~L^-1, ")")), 
#                    expression(paste("Submerged vegetation biomass (", gDw~m^-2, ")")),
#                    expression(paste("Benthivorous fish biomass (", gDw~m^-2, ")")))

# names(labels_states) <- names(summary_multiES$desired_states)
#
# Single ----------------------------------------------------------------
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
single_scatter <-lastpop$lastpop_single |> 
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

ggarrange(single_scatter, single_parallel,widths = c(1,0.7), align = 'h',
          labels = c("A)", "B)"), label.x = -0.01) |> 
  ggsave(filename = './output/plots/figure_2.png', width = 20, height = 9, unit = 'cm')

# Multi --------------------------------------------
summary_multi$desired_states
summary_multi$iter


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
  ggsave(filename = './output/plots/figure_S1.png', width = 20, height = 9, unit = 'cm')


# Two horizon --------------------------------------------
summary_twohorizon2$desired_states
summary_twohorizon2$iter


lastpop$lastpop_twohorizon2 |> 
  filter(fn_out <= 0) |> 
  mutate(mPLoadEpi_lag_first = ifelse(mPLoadEpi_lag < mPLoadEpi_lag2, mPLoadEpi_lag, mPLoadEpi_lag2),
         mPLoadEpi_lag_second = ifelse(mPLoadEpi_lag > mPLoadEpi_lag2, mPLoadEpi_lag, mPLoadEpi_lag2),
         mPLoadEpi_first = ifelse(mPLoadEpi_lag < mPLoadEpi_lag2, mPLoadEpi, mPLoadEpi2),
         mPLoadEpi_second = ifelse(mPLoadEpi_lag > mPLoadEpi_lag2, mPLoadEpi, mPLoadEpi2)) |> 
  pivot_longer(cols = c("fMarsh", "fMarsh_lag", "mPLoadEpi_lag_first", "mPLoadEpi_lag_second", "mPLoadEpi_first", "mPLoadEpi_second"),
               names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

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
             # size = value,
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
        panel.spacing.y = unit(c( rep( c( rep(0.3,2), 0.7), 2), rep(0.3,2)), "lines")) +
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
                      labels = c(labels_measures, 'oChlaEpi' = 'Chlorophyll a\nconcentration'))

ggsave(pathways_2h, filename = './output/plots/figure_3.png', width = 15, height = 15, unit = 'cm')

# MultiES -----------------------------------------
summary_multiES$iter
nrow(lastpop$lastpop_multiES)

lastpop$lastpop_multiES |> 
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
lastpop$lastpop_multiES |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_multiES) # what were the associated states

lastpop$lastpop_multiES |> 
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
  ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y',  nrow = 16, 
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
        # panel.spacing.y = unit(c( rep( c( rep(0.2,3), 0.6), 3), rep(0.2,3)),"lines")
        ) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
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
  scale_colour_manual(values = c('gold', 'grey20', 'seagreen', cols_measures), 
                      name = 'Lake state',
                      breaks = c('aDFish', 'oChlaEpi', "aDSubVeg", names(labels_measures)), 
                      labels = c('aDFish' = 'Fish biomass',
                                 'oChlaEpi' = 'Chlorophyll a\nconcentration',
                                 "aDSubVeg" = 'Submerge vegetation biomass', labels_measures))


# out of all the pathways how many achieve each of the indicator targets
summary_multiES$desired_states

multiES_ds <- data.frame(opt_var = names(summary_multiES$desired_states),
                         lower_range = sapply(summary_multiES$desired_states,
                                              function(x) min(x$target)),
                         upper_range = sapply(summary_multiES$desired_states,
                                              function(x) max(x$target)))

lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = sum(row_number())) |> # 1 = 1, 3 = 1+2, 4 = 1+3, 5 = 2+3, 6 = 1+2+3
  reframe(.by = c(total_achieve, ind_ID),
          n = n())

ggarrange(lastpopstate$lastpopstate_multiES |> 
            full_join(multiES_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = fManVeg, colour = fManVeg_lag)) +# remove colour because the range is so small?
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'plasma', begin = 0.9, end  = 0.5) +
            facet_wrap(~opt_var, scales = 'free_y', labeller = labeller(opt_var = labels_states)) +
            theme_bw()+
            scale_x_continuous(name = 'Fraction of vegetation removed')+
            labs(dictionary = labels_measures) ,
          
          lastpopstate$lastpopstate_multiES |> 
            full_join(multiES_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = mPLoadEpi, colour = mPLoadEpi_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'magma', begin = 0.1, end = 0.6) +
            facet_wrap(~opt_var, scales = 'free_y') +
            theme(legend.position = 'right') +
            theme_bw() +
            labs(dictionary = labels_measures) , 
          
          nrow = 2, align = 'hv') |> 
  ggsave(filename = './output/plots/figure_4.png', width = 20, height = 10, unit = 'cm')

lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(row_number()))) |> 
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
                     values = c("#43BBADFF", "#3D5296FF"), name = "Indicators achieved")

# what is the general pattern of the two groups?
lastpopstate$lastpopstate_multiES |> 
  full_join(multiES_ds, by = join_by(opt_var)) |> 
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
              values_from = value)

# Prioritised -------------------------------------
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


lastpop$lastpop_prioritisation |> 
  slice_min(fn_out) |> # lowest objective function
  select(runID) |> 
  left_join(lastpoppathways$lastpoppathways_prioritisation, by = join_by(runID == ID)) |> 
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
  ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y',  nrow = 16, 
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
        # panel.spacing.y = unit(c( rep( c( rep(0.2,3), 0.6), 3), rep(0.2,3)),"lines")
  ) +
  guides(color=guide_legend(nrow=2, byrow=TRUE)) +
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
  scale_colour_manual(values = c('gold', 'grey20', 'seagreen', cols_measures), 
                      name = 'Lake state',
                      breaks = c('aDFish', 'oChlaEpi', "aDSubVeg", names(labels_measures)), 
                      labels = c('aDFish' = 'Fish biomass',
                                 'oChlaEpi' = 'Chlorophyll a\nconcentration',
                                 "aDSubVeg" = 'Submerge vegetation biomass', labels_measures))

# which targets were met?
prioritisation_ds <- data.frame(opt_var = names(summary_prioritisation$desired_states),
                         lower_range = sapply(summary_prioritisation$desired_states,
                                              function(x) min(x$target)),
                         upper_range = sapply(summary_prioritisation$desired_states,
                                              function(x) max(x$target)))

lastpopstate$lastpopstate_prioritisation |> 
  full_join(prioritisation_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = sum(row_number())) |> # 1 = 1, 3 = 1+2, 4 = 1+3, 5 = 2+3, 6 = 1+2+3
  reframe(.by = c(total_achieve, ind_ID),
          n = n())

ggarrange(lastpopstate$lastpopstate_prioritisation |> 
            full_join(prioritisation_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = fManVeg, colour = fManVeg_lag)) +# remove colour because the range is so small?
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'plasma', begin = 0.9, end  = 0.5) +
            facet_wrap(~opt_var, scales = 'free_y', labeller = labeller(opt_var = labels_states)) +
            theme_bw()+
            scale_x_continuous(name = 'Fraction of vegetation removed')+
            labs(dictionary = labels_measures) ,
          
          lastpopstate$lastpopstate_prioritisation |> 
            full_join(prioritisation_ds, by = join_by(opt_var)) |> 
            ggplot(aes(y = out, x = mPLoadEpi, colour = mPLoadEpi_lag)) +
            geom_point() +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') +
            scale_colour_viridis_c(option = 'magma', begin = 0.1, end = 0.6) +
            facet_wrap(~opt_var, scales = 'free_y') +
            theme(legend.position = 'right') +
            theme_bw() +
            labs(dictionary = labels_measures) , 
          
          nrow = 2, align = 'hv') 

lastpopstate$lastpopstate_prioritisation |> 
  full_join(prioritisation_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(row_number()))) |> 
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
                     values = c("#43BBADFF", "#3D5296FF"), name = "Indicators achieved")

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
