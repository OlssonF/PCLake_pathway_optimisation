#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Example problem scenario comparisons - read in, plot, and summarise the 
## Date: 2026-01-23; updated 2026-06-07
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

library(tidyverse)
library(ggpubr)
library(GGally)
library(ggh4x)

convert_Pload <- function(before = 0.002, after) {
  
  change <- 100*(before - after)/before
  
  return(change)
}

source('scripts/R/optim_functions.R')
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


## Last iteration population, calculate the state values
lastpopstate <- lapply(list.files('output', 'lastpopstate_', full.names = TRUE),
                       read_delim, show_col_types = F)
names(lastpopstate) <- gsub('.csv', '', 
                            basename(list.files('output', 'lastpopstate_', full.names = TRUE)))


## Last iteration population, state values for the full ts
lastpoppathways <- lapply(list.files('output', 'lastpoppathways_', full.names = TRUE),
                          read_csv, show_col_types = F)

names(lastpoppathways) <- gsub('.csv', '', 
                               basename(list.files('output', 'lastpoppathways_', full.names = TRUE)))

## All populations, obj_function output
allpop <- lapply(list.files('output', 'allpops_', full.names = TRUE),
                 read_csv, show_col_types = F)

names(allpop) <- gsub('.csv', '', 
                      basename(list.files('output', 'allpops_', full.names = TRUE)))


# Aesthetics --------------------------------------------------------------
labels_measures <- c(expression(atop('P load', (gP~m^-2~d^-1))),
                     expression(atop('P load', (gP~m^-2~d^-1))),
                     # expression(paste("P load\n(mg ", L^-1, ")")),
                     # expression(paste("P load\n(mg ", L^-1, ")")), 
                     expression(atop("Marsh area", "(fraction)")),
                     expression(atop("Vegetation removed", "(fraction)")),
                     expression(atop("Day of vegetation", "removal (day of year)")),
                     "cDredInterval", "cDredStart",
                     expression(atop("Start P load", "reduction (year)")),
                     expression(atop("Start P load", "reduction (year)")),
                     expression(atop("Start marsh", " area (year)")),
                     expression(atop("Start vegetation", "removal (year)")),
                     expression(atop("Reduction in P load", "(%)")))



names(labels_measures) <- c(possible_measures$parameter, "mPLoadEpi_change")
cols_measures <- c("#3E4A89FF",
                   "#3E4A89FF",
                   "#6DCD59FF",
                   "#F89441FF",
                   "grey",
                   "grey", "grey",
                   "#3E4A89FF",
                   "#3E4A89FF",
                   "#6DCD59FF",
                   "#F89441FF",
                   "#3E4A89FF")

shapes_measures <- c(16,16,16,16,
                     1,1, 1,
                     2,2,2,2,16)

lines_measures <- c(rep('solid',4),
                    rep('dotted',3),
                    rep('dashed', 4),
                    'solid')

cols_states <- c("#C7EF34FF", "#36AAF9FF", "#7A0403FF")

labels_states <- c(oChlaEpi =expression(atop('Chlorophyll-a concentration', (mu*g~L^-1))),
                   aDSubVeg =expression(atop('Submerged vegetation biomass', (gDW~m^-2))),
                   aDFish   =expression(atop("Benthivorous fish biomass", (gDW~m^-2))))

labels_states_str <- c(oChlaEpi = "atop('Chlorophyll-a concentration', (mu*g~L^-1))",
                       aDSubVeg = "atop('Submerged vegetation biomass', (gDW~m^-2))",
                       aDFish   = "atop('Benthivorous fish biomass', (gDW~m^-2))")

labels_measures_str <- c(mPLoadEpi = "atop('P load', (gP~m^-2~d^-1))",
                         mPLoadEpi2 = "atop('P load', (gP~m^-2~d^-1))",
                         fMarsh = "atop('Marsh area', '(fraction)')",
                         fManVeg = "atop('Vegetation removed', '(fraction)')",
                         cDayManVeg1 = "atop('Day of vegetation', 'removal (day of year)')",
                         cDredInterval = 'cDredInterval', cDredStart =  'cDredStart',
                         mPLoadEpi_lag = "atop('Start P load', 'reduction(year)')",
                         mPLoadEpi_lag2 = "atop('Start P load', 'reduction (year)')",
                         fMarsh_lag = "atop('Start marsh', ' area (year)')",
                         fManVeg_lag = "atop('Start vegetation', 'removal (year)')",
                         mPLoadEpi_change = "atop('Reduction in P load', '(%)')")



#

## Simple example ---------------

simple_ds <- data.frame(opt_var = names(summary_simple$desired_states),
                        lower_range = sapply(summary_simple$desired_states,
                                             function(x) min(x$target)),
                        upper_range = sapply(summary_simple$desired_states,
                                             function(x) max(x$target))) |> 
  mutate(opt_var_val = row_number())

# calculate the mean parameter values of each population
allpop$allpops_simple |> 
  mutate(par1_change = convert_Pload(before = 0.002, after = par1)) |> 
  rename('mPLoadEpi_change' = par1_change, 'mPLoadEpi_lag' = par2) |>
  select(-par1) |> 
  pivot_longer(all_of(c('mPLoadEpi_lag', 'mPLoadEpi_change')),
               names_to = 'measure') |> 
  reframe(.by = c('iteration', 'measure'), 
          n = n(),
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          min = min(value),
          max = max(value),
          mad = mad(value)) |> 
  mutate(cv = sd/mean,
         rcv  = mad/median)

popmean_simple <-
  allpop$allpops_simple |> 
  mutate(par1_change = convert_Pload(before = 0.002, after = par1)) |> 
  select(par1_change, par2, iteration) |> 
  pivot_longer(cols = c(par1_change, par2), names_to = 'parameter') |> 
  reframe(.by = all_of(c('parameter', 'iteration')),
          mean = mean(value)) |> 
  pivot_wider(id_cols = iteration, names_from = parameter, values_from = mean) |> 
  mutate(type = 'mean')

allpop_p1 <- allpop$allpops_simple |> 
  mutate(par1_change = convert_Pload(before = 0.002, after = par1)) |> 
  select(par1_change, par2, iteration) |> 
  rename('mPLoadEpi_change' = par1_change, 'mPLoadEpi_lag' = par2) |>
  pivot_longer(cols = c(mPLoadEpi_change, mPLoadEpi_lag), names_to = 'parameter') |>
  mutate(iteration = as_factor(iteration)) |> 
  ggplot() + 
  geom_boxplot(aes(y=value, x = iteration)) + 
  facet_wrap(~parameter, scales = 'free_y', nrow = 2,
             labeller = labeller(parameter = as_labeller(labels_measures_str, 
                                                         label_parsed))) + 
  theme_bw()  +
  scale_x_discrete(breaks = c('1','5','10','15','21')) +
  labs(y = 'Parameter value')

# panel B
allpop_p2 <- allpop$allpops_simple |> 
  mutate(par1_change = convert_Pload(before = 0.002, after = par1),
         type = 'member') |> 
  full_join(popmean_simple, by = join_by(iteration, par2, par1_change, type)) |> 
  ggplot(aes(x=par1_change, y= par2)) + 
  geom_point(aes(colour = iteration, shape = type, size = type)) +
  scale_colour_viridis_c(begin = 0.2, end = 1, option = 'turbo') + 
  scale_shape_manual(values = c(16, 4)) +
  scale_size_manual(values = c(3,1.5)) +
  theme_bw() +
  coord_cartesian(xlim = c(0,100),
                  ylim = c(0,30)) +
  scale_x_continuous(expand = c(0,0))  +
  guides(shape = "none",
         size = "none") +
  labs(x = "Reduction in P load (%)",
       y = "Start P load reduction (year)")

ggarrange(allpop_p1, allpop_p2, widths = c(0.5,1), labels = c('A)', 'B)'))


# panel C
example_pathways <- lastpoppathways$lastpoppathways_simple |>
  filter(ID %in% seq(1,20,4)) |>
  mutate(.by = year, ID = row_number()) |> # renumber the pathways
  mutate(mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.002)) |>
  mutate(mPLoadEpi = mPLoadEpi_use) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi_use)) |> 
  
  select(all_of(c('ID', 'year', 'oChlaEpi', 'mPLoadEpi_change' ))) |> 
  pivot_longer(cols = !any_of(c('ID','year')), names_to = 'opt_var') |> 
  full_join(simple_ds, by = join_by(opt_var)) |>
  mutate(opt_var = factor(opt_var, levels = c("mPLoadEpi_change", 'oChlaEpi'))) |>
  ggplot(aes(x=year, y = value,
             colour = opt_var)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(ID, opt_var), scales = 'free_y',  nrow = 16,
                           dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'),
                           strip = strip_split(position = c('right', 'top'),
                                               text_y = element_text(),
                                               text_x = element_blank(),
                                               background_y = element_rect(),
                                               background_x = element_blank())) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.title = element_text(hjust = 0.5),
        panel.spacing.y = unit(c( rep( c( rep(0.2,), 0.6), 4), 0.2),"lines")) +
  scale_colour_manual(values = c(cols_measures, cols_states),
                      name = 'Lake state',
                      breaks = c(names(labels_measures), names(labels_states)),
                      labels = c(labels_measures, labels_states)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  guides(colour = guide_legend(direction = 'vertical', ncol = 2)) +
  facetted_pos_scales(y = list(opt_var == "mPLoadEpi_change" ~ scale_y_continuous(limits = c(0,110),
                                                                                  n.breaks = 2,
                                                                                  minor_breaks = c(50,90)),
                               opt_var == "oChlaEpi" ~ scale_y_continuous(limits = c(0,110),
                                                                          n.breaks = 2,
                                                                          minor_breaks = c(20,100)))) +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed') 


ggarrange(ggarrange(lastpopstate$lastpopstate_simple |> 
                      mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
                      ggplot(aes(x=mPLoadEpi_change, y=mPLoadEpi_lag)) +
                      geom_point() +
                      theme_bw() +
                      theme(legend.position="bottom") +
                      labs(dictionary = labels_measures),
                    
                    bestmemit$bestmemit_simple |> 
                      filter(fn_out == 0) |> 
                      mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
                      ggplot(aes(x=mPLoadEpi_change, y=mPLoadEpi_lag)) +
                      geom_point() +
                      theme_bw() +
                      theme(legend.position="bottom")   +
                      labs(dictionary = labels_measures),
                    labels = c('A)', 'B)')),
          example_pathways, 
          nrow = 2, heights = c(1,2), 
          labels = c('', 'C)'))


lastpop$lastpop_simple |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  filter(fn_out <= 0) |> 
  pivot_longer(all_of(c('mPLoadEpi', 'mPLoadEpi_lag', 'mPLoadEpi_change')),
               names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          median = median(value),
          sd = sd(value),
          min = min(value),
          max = max(value),
          mad = mad(value)) |> 
  mutate(cv = sd/mean,
         rcv  = mad/median)

# Figure S3
lastpoppathways$lastpoppathways_simple |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  filter(oChlaEpi < 20) |> 
  group_by(ID) |> 
  slice_min(year) |> 
  mutate(lag = year - mPLoadEpi_lag) |> 
  arrange(lag) |> ggplot() + 
  geom_point(aes(x = mPLoadEpi_lag, colour = lag, 
                 y = mPLoadEpi_change))+ 
  scale_color_viridis_c(name = 'Time from implementation to\ntarget achieved (years)') +
  theme_bw() +
  labs(dictionary = labels_measures) +
  theme(legend.position = 'top')


### -----------------------------------------------------------------#

## Full example ------------------------------------------------------
multiES_mega_ds <- data.frame(opt_var = names(summary_multiES_mega$desired_states),
                              lower_range = sapply(summary_multiES_mega$desired_states,
                                                   function(x) min(x$target)),
                              upper_range = sapply(summary_multiES_mega$desired_states,
                                                   function(x) max(x$target))) |> 
  mutate(opt_var_val = row_number())

colnames(allpop$allpops_multiES_mega) <- c('iteration', summary_multiES_mega$possible_measures$parameter)

# exploration of parameter space
allpop$allpops_multiES_mega |>
  mutate(mPLoadEpi_change = convert_Pload(0.002, mPLoadEpi)) |> select(-mPLoadEpi) |> 
  reframe(.by = iteration,
          across(-any_of(c('iteration')), .fns = list(mean = mean, sd = sd), .names = "{.col}!{.fn}")) |> 
  pivot_longer(-iteration, names_to = c('variable', 'stat'), names_sep = '!' ) |> 
  pivot_wider(names_from = stat, values_from = value) |> 
  mutate(cv = sd/mean,
         variable = as_factor(variable)) |>
  ggplot(aes(x=iteration, y = cv, colour = variable, linetype = variable, shape = variable)) +
  geom_point() + 
  geom_line() + 
  theme_bw() +
  theme(legend.position = 'top') +
  scale_colour_manual(values = cols_measures, 
                      name = 'Measure', 
                      breaks = names(labels_measures), 
                      labels = labels_measures) +
  scale_linetype_manual(values = lines_measures, 
                        name = 'Measure', 
                        breaks = names(labels_measures), 
                        labels = labels_measures)  +
  scale_shape_manual(values = shapes_measures, 
                     name = 'Measure', 
                     breaks = names(labels_measures), 
                     labels = labels_measures)



lastpop$lastpop_multiES_mega |> 
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
lastpop$lastpop_multiES_mega |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_multiES_mega) # what were the associated states

lastpop$lastpop_multiES_mega |> 
  slice_min(fn_out) |> # lowest objective function
  select(runID) |> 
  left_join(lastpoppathways$lastpoppathways_multiES_mega, by = join_by(runID == ID)) |> 
  
  mutate(.by = year, ID = row_number()) |> # renumber the pathways
  mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
         fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.002)) |>
  mutate(fManVeg = fManVeg_use,
         fMarsh = fMarsh_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi_use)) |> 
  
  select(all_of(c('ID', 'year', 'aDSubVeg', 'aDFish', 'oChlaEpi', 'mPLoadEpi_change', 'fManVeg', 'fMarsh'))) |> 
  pivot_longer(cols = !any_of(c('ID','year')), names_to = 'opt_var') |> 
  full_join(multiES_mega_ds, by = join_by(opt_var)) |>
  mutate(opt_var = factor(opt_var, levels = c("mPLoadEpi_change", 'fManVeg', 'fMarsh', 'oChlaEpi', 'aDSubVeg', 'aDFish'))) |>
  ggplot(aes(x=year, y = value,
             colour = opt_var)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(opt_var), scales = 'free_y',  nrow = 16,
                           dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'),
                           strip = strip_split(position = c('top'),
                                               text_y = element_text(),
                                               text_x = element_blank(),
                                               background_y = element_rect(),
                                               background_x = element_blank())) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.title = element_text(hjust = 0.5)) +
  scale_colour_manual(values = c(cols_measures, cols_states),
                      name = 'Lake state',
                      breaks = c(names(labels_measures), names(labels_states)),
                      labels = c(labels_measures, labels_states)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  guides(colour = guide_legend(direction = 'vertical', ncol = 2)) +
  facetted_pos_scales(y = list(opt_var == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                         n.breaks = 4),
                               opt_var == "mPLoadEpi_change" ~ scale_y_continuous(limits = c(0,110),
                                                                                  n.breaks = 2,
                                                                                  minor_breaks = c(50,90)),
                               opt_var == "oChlaEpi" ~ scale_y_continuous(limits = c(0,110),
                                                                          n.breaks = 2,
                                                                          minor_breaks = c(20,100)),
                               opt_var == "aDSubVeg" ~ scale_y_continuous(limits = c(0,110),
                                                                          n.breaks = 2),
                               opt_var == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                        n.breaks = 2),
                               opt_var == "fMarsh" ~ scale_y_continuous(limits = c(0,1),
                                                                        n.breaks = 2))) +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed') 

# Panel D
target_check_multiES_mega <- lastpopstate$lastpopstate_multiES_mega |> 
  full_join(multiES_mega_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES_mega, by = join_by(ID))


lastpopstate$lastpopstate_multiES_mega |> 
  full_join(multiES_mega_ds, by = join_by(opt_var)) |>
  rowwise() |>
  mutate(distance = range_obj(out, target = c(lower_range, upper_range))) |> 
  full_join(target_check_multiES_mega) |> 
  filter(distance != 0) |>  
  mutate(labs = ifelse(ind_ID == 1, 'chla', 
                       ifelse(ind_ID == 3, 'chlaM',
                              ifelse(ind_ID == 4, 'chlaF', NA)))) |> 
  ggplot(aes(x=labs, y= distance)) + 
  geom_boxplot() + 
  facet_wrap(~opt_var, scales= 'free',
             labeller = labeller(opt_var = as_labeller(labels_states_str, label_parsed))) +
  theme_bw()  +
  scale_x_discrete(labels = c('Chla only', 'Chla +\nMacrophytes', "Chla +\nFish"),
                   breaks = c('chla', 'chlaM', 'chlaF'),
                   name = "Attribute target achieved") 

# trade off among targets
lastpopstate$lastpopstate_multiES_mega |> 
  full_join(multiES_mega_ds, by = join_by(opt_var)) |>
  rowwise() |>
  mutate(distance = range_obj(out, target = c(lower_range, upper_range))) |> 
  full_join(target_check_multiES_mega) |> 
  filter(distance != 0) |>  group_by(ind_ID, opt_var) |> summarise(distance = mean(distance))

# how many in each group?
target_check_multiES_mega |> 
  distinct(ID, ind_ID) |> 
  reframe(.by = ind_ID, n = n())

target_check_multiES_mega |> 
  group_by(ind_ID, opt_var) |> 
  summarise(mean = mean(out),
            sd = sd(out),
            cv  = sd(out)/mean(out)) |> 
  arrange(opt_var)


target_check_multiES_mega |>
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  ggparcoord(scale = 'std', # rescaling
             columns = c(5:9, 12), # measures
             groupColumn = 3, 
             order = 'anyClass', #match(c("mPLoadEpi", 'fManVeg', 'fMarsh', "mPLoadEpi_lag", 'fManVeg_lag', 'fMarsh_lag'), colnames(target_check_multiES_mega)),
             showPoints = T, 
             # alphaLines = 0.4, 
             mapping = ggplot2::aes(linewidth = 0.7)) +
  theme_bw() +
  scale_x_discrete(labels = labels_measures, name = 'Measure') +
  scale_y_continuous(name = 'Normalised value') + 
  theme(axis.text.x = element_text(vjust = -0.5, hjust = 0.5)) +
  scale_color_manual(labels = c('Chla', 'Chla + Macrophytes', "Chla + Fish"),
                     values = viridis::turbo(n = 3, begin = 0.3, end = 0.9),
                     name = "Attribute target achieved") +
  scale_linewidth_identity() +
  theme(legend.position = 'top')



# alt figure
ggarrange(lastpopstate$lastpopstate_multiES_mega |> 
            mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
            full_join(multiES_mega_ds, by = join_by(opt_var)) |> 
            ggplot() + 
            geom_point(aes(colour=mPLoadEpi_change, x = mPLoadEpi_lag, y = out)) +
            facet_wrap(~opt_var, scales = 'free', 
                       labeller = labeller(opt_var = as_labeller(labels_states_str, label_parsed))) +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') + 
            scale_colour_viridis_c(option = 'plasma', begin = 0, end  = 0.5) +
            theme_bw(base_size = 12) +
            theme(legend.position = "right",
                  legend.title.position = 'right',
                  legend.justification = "left",
                  legend.box.margin = margin(l = 0.5, unit = "cm"))+
            labs(dictionary = c(labels_measures, 'out' = 'Final state value')),
          
          lastpopstate$lastpopstate_multiES_mega |> 
            full_join(multiES_mega_ds, by = join_by(opt_var)) |> 
            ggplot() + 
            geom_point(aes(colour=fManVeg, x = fManVeg_lag, y = out)) +
            facet_wrap(~opt_var, scales = 'free', 
                       labeller = labeller(opt_var = as_labeller(labels_states_str, label_parsed))) +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') + 
            scale_colour_viridis_c(option = 'viridis', begin = 0.3, end = 0.9) +
            theme_bw(base_size = 12) +
            theme(legend.position = "right",
                  legend.title.position = 'right',
                  legend.justification = "left",
                  legend.box.margin = margin(l = 0.5, unit = "cm"))+
            labs(dictionary = c(labels_measures, 'out' = 'Final state value')),
          
          lastpopstate$lastpopstate_multiES_mega |> 
            full_join(multiES_mega_ds, by = join_by(opt_var)) |> 
            ggplot() + 
            geom_point(aes(colour=fMarsh, x = fMarsh_lag, y = out)) +
            facet_wrap(~opt_var, scales = 'free', 
                       labeller = labeller(opt_var = as_labeller(labels_states_str, label_parsed))) +
            geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
            geom_hline(aes(yintercept = upper_range), linetype = 'dashed') + 
            scale_colour_viridis_c(option = 'plasma', begin = 0.9, end  = 0.5) +
            theme_bw(base_size = 12) +
            theme(legend.position = "right",
                  legend.title.position = 'right',
                  legend.justification = "left",
                  legend.box.margin = margin(l = 0.5, unit = "cm"))+
            labs(dictionary = c(labels_measures, 'out' = 'Final state value')),
          nrow = 3, align = 'hv'
)


# example pathways:pathway to achieve fish
example_IDs <- target_check_multiES_mega |> group_by(ind_ID) |>  slice_head() 
# objectively "best pathway"
best_ID <- lastpop$lastpop_multiES_mega |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_multiES_mega) |> # what were the associated states
  distinct(ID) |>
  inner_join(target_check_multiES_mega) |> 
  slice_head() |> 
  mutate(ind_ID = 'best')


lastpoppathways$lastpoppathways_multiES_mega |> 
  inner_join(bind_rows(example_IDs, best_ID)) |> # get the best and an example of each
  mutate(.by = year, name = ifelse(ind_ID == 1, 'Chla', 
                                   ifelse(ind_ID == 3, 'Chl a + macrophytes',
                                          ifelse(ind_ID == 4, 'Chl a +  fish', 'Best')))) |> # renumber the pathways
  mutate(fManVeg_use = ifelse(fManVeg_lag < year, fManVeg, 0),
         fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.002)) |>
  mutate(fManVeg = fManVeg_use,
         fMarsh = fMarsh_use,
         mPLoadEpi = mPLoadEpi_use) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi_use)) |> 
  
  select(all_of(c('name', 'year', 'aDSubVeg', 'aDFish', 'oChlaEpi', 'mPLoadEpi_change', 'fManVeg', 'fMarsh'))) |> 
  pivot_longer(cols = !any_of(c('name','year')), names_to = 'opt_var') |> 
  full_join(multiES_mega_ds, by = join_by(opt_var)) |>
  mutate(opt_var = factor(opt_var, levels = c("mPLoadEpi_change", 'fManVeg', 'fMarsh', 'oChlaEpi', 'aDSubVeg', 'aDFish'))) |>
  ggplot(aes(x=year, y = value,
             colour = opt_var)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  ggh4x::facet_nested_wrap(vars(name, opt_var), scales = 'free_y',  nrow = 6,
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           strip = strip_split(position = c('right', 'top'),
                                               text_y = element_text(),
                                               text_x = element_blank(),
                                               background_y = element_rect(),
                                               background_x = element_blank())) +
  theme_bw(base_size = 12) +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.title = element_text(hjust = 0.5))  +
  scale_colour_manual(values = c(cols_measures, cols_states),
                      name = 'Lake state',
                      breaks = c(names(labels_measures), names(labels_states)),
                      labels = c(labels_measures, labels_states)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  guides(colour = guide_legend(direction = 'vertical', ncol = 6)) +
  facetted_pos_scales(y = list(opt_var == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                         n.breaks = 4),
                               opt_var == "mPLoadEpi_change" ~ scale_y_continuous(limits = c(0,110),
                                                                                  n.breaks = 2,
                                                                                  minor_breaks = c(50,90)),
                               opt_var == "oChlaEpi" ~ scale_y_continuous(limits = c(0,110),
                                                                          n.breaks = 2,
                                                                          minor_breaks = c(20,100)),
                               opt_var == "aDSubVeg" ~ scale_y_continuous(limits = c(0,115),
                                                                          n.breaks = 2),
                               opt_var == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                        n.breaks = 2),
                               opt_var == "fMarsh" ~ scale_y_continuous(limits = c(0,1),
                                                                        n.breaks = 2)))  +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed') 


#
label_pairs  <- as_labeller(labels_measures_str, 
                                      label_parsed)

lastpopstate$lastpopstate_multiES_mega |>
  slice_min(out, prop = 0.1) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  ggpairs(columns = c(3:7, 10), 
          diag = list('continuous' = 'densityDiag'),
          upper = list('continuous' = 'cor'),
          labeller = label_pairs) + 
  theme_bw() 


# individual pathway timings
target_check_multiES_mega |> 
 select(ID, mPLoadEpi_lag, fMarsh_lag, fManVeg_lag) |> 
  pivot_longer(-ID, names_to = 'parameter', values_to = 'lag') |> 
  ggplot(aes(x=ID, y = lag, colour = parameter)) +
  geom_point()

target_check_multiES_mega |> 
 distinct(ind_ID, ID, mPLoadEpi_lag, fMarsh_lag, fManVeg_lag) |> 
  pivot_longer(all_of(c('mPLoadEpi_lag', 'fMarsh_lag', 'fManVeg_lag')), names_to = 'parameter', values_to = 'lag') |> 
  arrange(lag) |> 
  group_by(ID) |> 
  slice_head(n = 1) |> 
  group_by(parameter, ind_ID) |> 
  summarise(n())
