#--------------------------------------#
## Project: Pathway Optimisation Framework
## Script purpose: Example problem scenario comparisons - read in, plot, and summarise the 
## Date: 2026-01-23; updated 2026-06-07
## Author: Freya Olsson
# Created with R version 4.5.2 (2025-10-31 ucrt)
#--------------------------------------#

library(tidyverse)
library(ggpubr)
library(ggh4x)

source('R/optim_functions.R')
source('scripts/example model initialisation.R')
# Read in data ------------------------------------------------------------
setwd(here::here())
possible_measures <- read_csv('input/possible_measures.csv')


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
{labels_measures <- c(expression(atop('P load', (gP~m^-2~d^-1))),
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
labels_measures <- labels_measures[sort(names(labels_measures))]
labels_measures <- labels_measures[order(grepl("_lag", names(labels_measures)))]
}

cols_measures <- c( "grey","grey", "grey",
                    "#F89441FF",  
                    "#6DCD59FF",
                    "#3E4A89FF","#3E4A89FF","#3E4A89FF",
                    "#F89441FF",  
                    "#6DCD59FF",
                    "#3E4A89FF",
                    "#3E4A89FF")

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
                         mPLoadEpi_lag = "atop('Start P load', 'reduction (year)')",
                         mPLoadEpi_lag2 = "atop('Start P load', 'reduction (year)')",
                         fMarsh_lag = "atop('Start marsh', ' area (year)')",
                         fManVeg_lag = "atop('Start vegetation', 'removal (year)')",
                         mPLoadEpi_change = "atop('Reduction in P load', '(%)')")

#

# Figure 1 - conceptual figure ------------------------#
# Figure 2 - simple example ---------------------------

simple_ds <- data.frame(opt_var = names(summary_simple$desired_states),
                        lower_range = sapply(summary_simple$desired_states,
                                             function(x) min(x$target)),
                        upper_range = sapply(summary_simple$desired_states,
                                             function(x) max(x$target))) |> 
  mutate(opt_var_val = row_number())

# panel A
allpop_p1 <- allpop$allpops_simple |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  select(mPLoadEpi_change, mPLoadEpi_lag, iteration) |> 
  pivot_longer(cols = c(mPLoadEpi_change, mPLoadEpi_lag), names_to = 'parameter') |>
  mutate(iteration = as_factor(iteration)) |> 
  ggplot() + 
  geom_boxplot(aes(y=value, x = iteration)) + 
  facet_wrap(~parameter, scales = 'free_y', nrow = 2,
             labeller = labeller(parameter = as_labeller(labels_measures_str, 
                                                         label_parsed))) + 
  theme_bw()  +
  scale_x_discrete(breaks = c('1','5','10','15','21')) +
  labs(y = 'Parameter value', x = 'Iteration')

# panel B

example_pathways_df <- allpop$allpops_simple  |>
  mutate(success = ifelse(obj == 0, T, F),
         mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  filter(success == T) |> 
  mutate(example = ifelse((mPLoadEpi_change == max(mPLoadEpi_change)) |
                            (mPLoadEpi_change == min(mPLoadEpi_change)),
                          T, F)) |> filter(example == T) |> 
  select(iteration, mPLoadEpi, mPLoadEpi_lag) |> 
  mutate(ID = c('example 1', 'example 2'))

allpop_p2 <- allpop$allpops_simple |>
  mutate(success = ifelse(obj == 0, T, F),
         mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi),
         example = ifelse(mPLoadEpi %in% example_pathways_df$mPLoadEpi,
                          T, F)) |>
  mutate(shape = ifelse(success == T & example == T, 's1',
                        ifelse(success == T & example == F, 's2', 's3'))) |> 
  ggplot(aes(y=mPLoadEpi_change, x = mPLoadEpi_lag, 
             colour = obj, 
             shape = shape, size = shape)) +
  geom_point() +
  geom_label(data = ~ subset(.x, example == TRUE),
             aes(label = c('B-2', 'B-1')),
             label.size = 0.2, nudge_y = 3,
             nudge_x =1, 
             fill = "white")+
  scale_shape_manual(values = c(19, 1, 4)) +
  scale_colour_viridis_c(begin = 0.85, end = 0, option = 'plasma', 
                         name = 'Objective function value\n0=target achieved') + 
  scale_size_manual(values = c(3, 2, 1.5)) +
  scale_x_continuous(expand = c(0.01,0.1))  +
  scale_y_continuous(expand = c(0.01,0.1))  +
  theme_bw() +
  coord_cartesian(xlim = c(0,30),
                  ylim = c(0,100)) +
  guides(shape = "none", size = 'none') +
  labs(y = "Reduction in P load (%)",
       x = "Start P load reduction (year)") +
  theme(legend.position = 'top', 
        legend.title.position = 'left', 
        legend.title = element_text(hjust = 0.5))



# panel C

example_pathways_ls <- example_pathways_df |> 
  select(-ID) |> 
  group_split(iteration, .keep = F) |> 
  map(~ unlist(.x)) |> 
  set_names(c('example 1','example 2'))

# run the pathways with the example parameter sets
example_runs <- map(.x = example_pathways_ls,
                    ~run_pathway(val_pars = .x, name_pars = names(.x), 
                                 current_val = c(0.002, 0), initial_conditions = equilibrium_states)) |> 
  set_names(c('example 1','example 2')) |> 
  list_rbind(names_to = 'ID') |> 
  mutate(year = floor((time-1)/365) + 1,
         doy = yday(as_date(time - (year * 365) + 364, origin = '2025-01-01'))) |> 
  filter(doy %in% 50:300, year %in% 1:30) |> 
  select(c('year', 'ID', names(summary_simple$desired_states))) |> 
  group_by(year, ID) |> 
  summarise(across(any_of(names(summary_simple$desired_states)), max))

setwd(here::here()) # weird things happen when you run the pathways

# plot the example pathways
example_p3 <- example_runs |>
  full_join(example_pathways_df, by = 'ID') |> 
  mutate(mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.002)) |>
  mutate(mPLoadEpi = mPLoadEpi_use) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi_use)) |> 
  
  select(all_of(c('ID', 'year', 'oChlaEpi', 'mPLoadEpi_change'))) |> 
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
        panel.spacing.y = unit(c( rep( c( rep(0.2,), 0.8), 1), 0.2),"lines")) +
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
                                                                          minor_breaks = c(50,100)))) +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed') 

# arrange and save
fig2 <- ggarrange(ggarrange(allpop_p1, allpop_p2, widths = c(0.6,1), labels = c('A)', 'B)')),
          ggarrange(NULL, example_p3, NULL, widths = c(0.1,1, 0.1), ncol = 3, 
                    labels = c('', 'C)', '')),
          nrow = 2, heights = c(1,1))
ggsave(fig2, filename = 'output/plots/ms/Figure2.jpg', height = 20, width = 20, units = 'cm')


# Figure 3 - correlations among measures ---------------
cor_vars <- lastpopstate$lastpopstate_multiES_constrained |> 
  select(any_of(c(possible_measures$parameter))) |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  select(-mPLoadEpi)

cor_varnames <- names(cor_vars)
cor_varnames <- names(labels_measures)[which(names(labels_measures) %in% cor_varnames)]

cor_results <- expand_grid(var1 = cor_varnames,
                           var2 = cor_varnames)  |> 
  mutate(var1 = factor(var1, levels = names(labels_measures)),
         var2 = factor(var2, levels = names(labels_measures))) |> 
  filter(var1 != var2) |> 
  filter(as.numeric(var1) <= as.numeric(var2)) |> # keep unique pairs only
  mutate(cor = map2_dbl(var1, var2, ~ cor(cor_vars[[.x]], cor_vars[[.y]], method = "spearman", )),
         p_value = format(round(map2_dbl(var1, var2, ~ cor.test(cor_vars[[.x]], cor_vars[[.y]], method = "spearman")$p.value), 
                                digits = 3),
                          scientific = F)) |> 
  mutate(sig_symbol = case_when(p_value < 0.001 ~ "***",
                                p_value < 0.01  ~ "**",
                                p_value < 0.05  ~ "*",
                                TRUE ~ "")) |>
  arrange(var1, var2)


corplot <- cor_results |> 
  ggplot(aes(y= var2, x = var1)) + 
  geom_tile(aes(fill = cor)) + 
  geom_text(aes(label = paste(round(cor, 3), sig_symbol))) +
  scale_fill_continuous(palette = 'RdBu', limits = c(-1, 1), name = 'Spearman rank correlation') + 
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank()) +
  scale_y_discrete(labels = labels_measures,
                   breaks = names(labels_measures),
                   expand = expansion(mult = c(0.1, 0.1))) +
  theme(panel.background = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        legend.position = 'top',
        legend.title.position = 'top',
        legend.title = element_text(hjust = 0.5)) +
  scale_x_discrete(position = 'top', 
                   labels = labels_measures,
                   expand = expansion(mult = c(0.1, 0.1)))

ggsave(corplot, filename = 'output/plots/ms/Figure3.png',
       height = 18, width = 18, units = 'cm')

# Figure 4 - boxplot comparisons ------------------------------
multiES_constrained_ds <- data.frame(opt_var = names(summary_multiES_constrained$desired_states),
                                     lower_range = sapply(summary_multiES_constrained$desired_states,
                                                          function(x) min(x$target)),
                                     upper_range = sapply(summary_multiES_constrained$desired_states,
                                                          function(x) max(x$target))) |> 
  mutate(opt_var_val = c(1,2,4))

# which pathways achieved which targets?
target_check_multiES_constrained <- lastpopstate$lastpopstate_multiES_constrained |> 
  full_join(multiES_constrained_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES_constrained, by = join_by(ID))

figure4 <- 
  target_check_multiES_constrained |>
  select(any_of(c('ID', 'ind_ID', possible_measures$parameter))) |> 
  distinct() |> 
  mutate(mPLoadEpi_change = convert_Pload(before = 0.002, after = mPLoadEpi)) |> 
  pivot_longer(any_of(c('mPLoadEpi_change', "fMarsh", "fManVeg", 
                        "mPLoadEpi_lag","fMarsh_lag", "fManVeg_lag")),
               names_to = 'parameter', values_to = 'value') |> 
  ggplot(aes(x=ind_ID, y = value)) +
  geom_boxplot() +
  facet_wrap(~parameter, nrow = 3, scales = 'free_y', 
             labeller = labeller(parameter = as_labeller(labels_measures_str, 
                                                         label_parsed))) +
  scale_x_discrete(labels = c('Chla only', 'Macrophytes\nonly', 'Fish\nonly', 'Chla +\nMacrophytes', "Chla +\nFish", 'All'),
                   breaks = c(1, 2, 4, 3, 5, 7),
                   name = "Attribute target achieved") +
  theme_bw() +
  facetted_pos_scales(y = list(scale_y_continuous(limits = c(0,1)), scale_y_continuous(limits = c(0,30)),
                               scale_y_continuous(limits = c(0,1)), scale_y_continuous(limits = c(0,30)),
                               scale_y_continuous(limits = c(0,100)), scale_y_continuous(limits = c(0,30)))) +
  labs(y = 'Parameter value')

ggsave(figure4, filename = 'output/plots/ms/Figure4.jpg', height = 15, width = 15, units = 'cm')

# Figure 5 - boxplot of distance from target ---------------
# per target_achieved group, add point for the best_fit variant in the figure. 

f5_p1 <-
  lastpopstate$lastpopstate_multiES_constrained |> 
  full_join(multiES_constrained_ds, by = join_by(opt_var)) |>
  rowwise() |>
  mutate(distance = range_obj(out, target = c(lower_range, upper_range))) |> 
  full_join(target_check_multiES_constrained) |> 
  filter(distance != 0) |>
  ggplot(aes(y=ind_ID, x= distance)) + 
  geom_boxplot() + 
  facet_wrap(~opt_var, scales= 'free', nrow = 3,
             labeller = labeller(opt_var = as_labeller(labels_states_str, label_parsed))) +
  theme_bw()  +
  scale_x_continuous(name = 'Absolute distance to target') +
  scale_y_discrete(labels = c('Chla only', 'Macrophytes\nonly', 'Fish\nonly', 'Chla +\nMacrophytes', "Chla +\nFish", 'All'),
                   breaks = c(1, 2, 4, 3, 5, 7),
                   name = "Attribute target achieved") +
  ggh4x::force_panelsizes(rows = c(3,1, 1))

#Panel B shows "best_fit" pathway
f5_p2 <- lastpop$lastpop_multiES_constrained |> 
  slice_min(fn_out) |> # lowest objective function
  select(runID) |> 
  left_join(lastpoppathways$lastpoppathways_multiES_constrained, by = join_by(runID == ID)) |> 
  filter(year %in% 1:30) |> 
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
  full_join(multiES_constrained_ds, by = join_by(opt_var)) |>
  mutate(opt_var = factor(opt_var, levels = c("mPLoadEpi_change", 'fManVeg', 'fMarsh', 'oChlaEpi', 'aDSubVeg', 'aDFish'))) |>
  ggplot(aes(x=year, y = value,
             colour = opt_var)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 0.7) +
  ggh4x::facet_nested_wrap(vars(ID, opt_var), scales = 'free_y',  nrow = 6,
                           dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'),
                           strip = strip_split(position = c('top'),
                                               text_y = element_text(),
                                               text_x = element_blank(),
                                               background_y = element_rect(),
                                               background_x = element_blank())) +
  theme_bw()  +
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
                               opt_var == "aDSubVeg" ~ scale_y_continuous(limits = c(0,100),
                                                                          n.breaks = 2),
                               opt_var == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                        n.breaks = 2),
                               opt_var == "fMarsh" ~ scale_y_continuous(limits = c(0,1),
                                                                        n.breaks = 2))) +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed')+
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.title = element_text(hjust = 0.5),
        panel.spacing.y = unit(c( rep( c( rep(0.2, 2), 0.8), ), rep(0.2, 2)),"lines"))

figure5 <- cowplot::plot_grid(f5_p1,                  
                              f5_p2, 
                              align = 'h', axis = 'b',
                              labels = c('A)', 'B)'), rel_widths = c(0.9,1))
ggsave(figure5, filename = 'output/plots/ms/Figure5.jpg', height = 15, width = 20, units = 'cm')


# Supplementary figures ------------------------------

## Figure S1 - equilibrium states spin-up -------------
# run during source('scripts/example model initialisation.R')

## Figure S2- objective function ----------------------
# show the output of the objective function(s)
target_ex <- c(6,10)
result_ex <- 2:14

obj_output <- map_dfr(result_ex, ~ tibble(range = range_obj(.x, target = target_ex),
                                          below = below_obj(.x, target = mean(target_ex)),
                                          above = above_obj(.x, target = mean(target_ex)), 
                                          exact = exact_obj(.x, target = mean(target_ex)))) |> 
  mutate(input = result_ex) |> 
  pivot_longer(-input, names_to = 'method', values_to = 'output') |> 
  full_join(data.frame(method = c('range', 'below', 'above', 'exact'),
                       lower = c(6,NA, NA, NA),
                       upper = c(10, NA, NA,NA),
                       val = c(NA, 8,8,8))) |> 
  ggplot(aes(x=input, y=output)) +
  geom_point() +
  facet_wrap(~method, scales = 'free') +
  theme_bw() +
  geom_vline(aes(xintercept = val), linetype = 'dashed') +
  geom_vline(aes(xintercept = lower), linetype = 'dotted') +
  geom_vline(aes(xintercept = upper), linetype = 'dotted') +
  labs(x = 'Evaluated result', y = 'Output from objective function')

ggsave(plot = obj_output, filename ='output/plots/ms/FigureS2.jpg',
       height = 10, width = 10, units = 'cm')

## Figure S3 - convergence of parameters, MoMm -----------
shapes_measures <- c(1,1,1,
                     16,16,16,16,16,
                     2,2,2,2)

lines_measures <- c(rep('dotted',3),
                    rep('solid',5),
                    rep('dashed', 4))
# exploration of parameter space
cv_evo <- ggarrange(allpop$allpops_multiES_constrained |>
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
                      guides(colour = guide_legend(direction = 'vertical', ncol = 2)) + 
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
                                         labels = labels_measures)+
                      scale_x_continuous(limits = c(0,50), name = 'Iteration') +
                      scale_y_continuous(limits = c(0.15,0.8), name = 'CV'),
                    allpop$allpops_multiES_compromise |>
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
                      guides(colour = guide_legend(direction = 'vertical', ncol = 2)) + 
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
                                         labels = labels_measures) +
                      scale_x_continuous(limits = c(0,50), name = 'Iteration') +
                      scale_y_continuous(limits = c(0.15,0.8), name = 'CV'),
                    nrow = 2, common.legend = T, labels = c('A)', 'B)'), vjust = 0, hjust = -0.5
)

ggsave(cv_evo,
       filename = 'output/plots/ms/FigureS3.jpg', height = 15, width = 15, units = 'cm')


## Figure S4 - successful pathways compromise ------------
labels_measures_str_oneline <- c(mPLoadEpi = "P load (gP~m^-2~d^-1)",
                                 mPLoadEpi2 = "P load (gP~m^-2~d^-1)",
                                 fMarsh = "Marsh area (fraction)",
                                 fManVeg = "Vegetation removed (fraction)",
                                 cDayManVeg1 = "atop('Day of vegetation', 'removal (day of year)')",
                                 cDredInterval = 'cDredInterval', cDredStart =  'cDredStart',
                                 mPLoadEpi_lag = "Start P load reduction (year)",
                                 mPLoadEpi_lag2 = "Start P load reduction (year)",
                                 fMarsh_lag = "Start marsh area (year)",
                                 fManVeg_lag = "Start vegetation removal (year)",
                                 mPLoadEpi_change = "Reduction in P load (%)")

multiES_compromise_ds <- data.frame(opt_var = names(summary_multiES_compromise$desired_states),
                                    lower_range = sapply(summary_multiES_compromise$desired_states,
                                                         function(x) min(x$target)),
                                    upper_range = sapply(summary_multiES_compromise$desired_states,
                                                         function(x) max(x$target))) |> 
  mutate(opt_var_val = c(1,2,4))


# which pathways achieved with targets
target_check_multiES_compromise <- lastpopstate$lastpopstate_multiES_compromise |> 
  full_join(multiES_compromise_ds, by = join_by(opt_var)) |> 
  filter(between(out, lower_range, upper_range)) |> 
  reframe(.by = ID,
          total_achieve = n(), 
          ind_ID = as_factor(sum(opt_var_val))) |> 
  full_join(lastpopstate$lastpopstate_multiES_compromise, by = join_by(ID))

success_pathways <- 
  lastpop$lastpop_multiES_compromise |> filter(fn_out == 0) |> 
  pull(runID)

figure_s4 <- lastpoppathways$lastpoppathways_multiES_compromise |> 
  filter(ID %in% success_pathways) |> 
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
  full_join(multiES_compromise_ds, by = join_by(opt_var)) |>
  mutate(opt_var = factor(opt_var, levels = c("mPLoadEpi_change", 'fManVeg', 'fMarsh', 'oChlaEpi', 'aDSubVeg', 'aDFish'))) |> 
  ggplot(aes(x=year, y = value,
             colour = opt_var)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 0.6) +
  ggh4x::facet_nested_wrap(vars(ID, opt_var), scales = 'free_y',  nrow = 6,
                           dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'),
                           strip = strip_split(position = c('top'),
                                               text_y = element_text(),
                                               text_x = element_blank(),
                                               background_y = element_rect(),
                                               background_x = element_blank())) +
  theme_bw()  +
  scale_colour_manual(values = c(cols_measures, cols_states),
                      name = 'Lake state',
                      breaks = c(names(labels_measures), names(labels_states)),
                      labels = c(labels_measures_str_oneline, labels_states)) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  guides(colour = guide_legend(direction = 'vertical', ncol = 2)) +
  facetted_pos_scales(y = list(opt_var == "fManVeg" ~ scale_y_continuous(limits = c(0,1),
                                                                         n.breaks = 4),
                               opt_var == "mPLoadEpi_change" ~ scale_y_continuous(limits = c(0,100),
                                                                                  n.breaks = 2,
                                                                                  minor_breaks = c(50,90)),
                               opt_var == "oChlaEpi" ~ scale_y_continuous(limits = c(0,100),
                                                                          n.breaks = 2,
                                                                          minor_breaks = c(20,100)),
                               opt_var == "aDSubVeg" ~ scale_y_continuous(limits = c(0,130),
                                                                          n.breaks = 2),
                               opt_var == "aDFish" ~ scale_y_continuous(limits = c(0,10),
                                                                        n.breaks = 2),
                               opt_var == "fMarsh" ~ scale_y_continuous(limits = c(0,1),
                                                                        n.breaks = 2))) +
  labs(x='Year', y = '') +
  geom_hline(aes(yintercept = lower_range), linetype = 'dashed') +
  geom_hline(aes(yintercept = upper_range), linetype = 'dashed')+
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top',
        legend.key.spacing.y = unit(0.4, "lines"), 
        legend.title = element_text(hjust = 0.5),
        panel.spacing.y = unit(rep(c( c( rep(0.2, 2), 0.8), rep(0.2, 2))),"lines")) 

ggsave(figure_s4, create.dir = T, filename = 'output/plots/ms/FigureS4.png',
       height = 12, width = 15, unit ='cm')

# Stats and summaries--------------------
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
          max = max(value)) |> 
  mutate(cv = sd/mean)


lastpop$lastpop_multiES_constrained |> 
  pivot_longer(mPLoadEpi:fManVeg_lag, names_to = 'measure') |> 
  reframe(.by = 'measure', 
          n = n(),
          mean = mean(value),
          sd = sd(value),
          min = min(value),
          max = max(value)) |> 
  mutate(cv = sd/mean)

# how many in each group?
target_check_multiES_constrained |> 
  distinct(ID, ind_ID) |> 
  reframe(.by = ind_ID, n = n())

target_check_multiES_compromise |> 
  distinct(ID, ind_ID) |> 
  reframe(.by = ind_ID,
          n())

target_check_multiES_constrained |> 
  group_by(ind_ID, opt_var) |> 
  summarise(mean = mean(out),
            sd = sd(out),
            cv  = sd(out)/mean(out)) |> 
  arrange(opt_var)

# objectively "best pathway"
lastpop$lastpop_multiES_constrained |> 
  slice_min(fn_out) |> # lowest objective function
  inner_join(lastpopstate$lastpopstate_multiES_constrained) # what were the associated states

# trade off among targets
lastpopstate$lastpopstate_multiES_constrained |> 
  full_join(multiES_constrained_ds, by = join_by(opt_var)) |>
  rowwise() |>
  mutate(distance = range_obj(out, target = c(lower_range, upper_range))) |> 
  full_join(target_check_multiES_constrained) |> 
  filter(distance != 0) |>  group_by(ind_ID, opt_var) |> summarise(distance = mean(distance))

## timing of measures
latest_measures <- target_check_multiES_constrained |> 
  distinct(ind_ID, ID, mPLoadEpi_lag, fMarsh_lag, fManVeg_lag) |> 
  pivot_longer(all_of(c('mPLoadEpi_lag', 'fMarsh_lag', 'fManVeg_lag')), names_to = 'parameter', values_to = 'lag') |> 
  arrange(desc(lag)) |># want the last
  group_by(ID) |> 
  slice_head(n = 1) |> 
  group_by(parameter) |> 
  summarise(n = n()) |> 
  pivot_wider(names_from = parameter, values_from = n) 

earliest_measures <- target_check_multiES_constrained |> 
  distinct(ind_ID, ID, mPLoadEpi_lag, fMarsh_lag, fManVeg_lag) |> 
  pivot_longer(all_of(c('mPLoadEpi_lag', 'fMarsh_lag', 'fManVeg_lag')), names_to = 'parameter', values_to = 'lag') |> 
  arrange(lag) |># want the first
  group_by(ID) |> 
  slice_head(n = 1) |> 
  group_by(parameter) |> 
  summarise(n = n()) |> 
  pivot_wider(names_from = parameter, values_from = n) 

earliest_measures |> 
  chisq.test()

latest_measures |> 
  chisq.test()
