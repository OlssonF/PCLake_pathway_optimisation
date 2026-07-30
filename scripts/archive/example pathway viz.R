
library(tidyverse)
library(ggpubr)
library(ggh4x)

# Read in data ------------------------------------------------------------
setwd(here::here())
possible_measures <- read_csv('possible_measures.csv')

## Rdata output of the optimisation
rds_files <-  list.files('output', '.RData', full.names = T)[1]

for (input_file in rds_files) {
  assign(gsub('.RData', '', basename(input_file)),
         read_rds(input_file))
}

## Best member of each iteration
bestmemit <- read_csv(list.files('output', 'bestmemit', full.names = TRUE)[1], show_col_types = F)

## Last iteration population, obj_function output
lastpop <- read_csv(list.files('output', 'lastpop_', full.names = TRUE)[1], show_col_types = F)

## Last iteration population, obj_function output
lastpopstate <- read_delim(list.files('output', 'lastpopstate_', full.names = TRUE)[1], show_col_types = F)

## Last iteration population, obj_function output
lastpoppathways <- read_csv(list.files('output', 'lastpoppathways_', full.names = TRUE)[1], show_col_types = F)


# Plotting ----------------------------------------------------------------

lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  ggplot(aes(x=year, y = oChlaEpi, group = ID)) +
  geom_line()


lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways
  ggplot(aes(x=year, y = ID, linewidth = oChlaEpi, group = ID)) +
  geom_line()

lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01, 1)) |> 
  ggplot(aes(x=year, colour = mPLoadEpi_use, y = ID, group = ID)) +
  geom_line()


lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0)) |> 
  ggplot(aes(x=year, colour = fMarsh_use, y = ID, group = ID)) +
  geom_line()



lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01 * 50, 0),
         oChlaEpi_use = oChlaEpi/200) |> 
  select(all_of(c('ID','year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID', 'year'))) |> 
  ggplot(aes(x=year, y = name, 
             size = value, 
             group = name)) +
  # geom_line() +
  geom_point()

lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01 * 50, 0),
         oChlaEpi_use = oChlaEpi/200) |> 
  select(all_of(c('ID','year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID', 'year'))) |> 
  ggplot(aes(y=year, x = as_factor(ID), 
             size = value, 
             colour = name,
             group = interaction(ID, name))) +
  # geom_line() +
  geom_point(position = position_dodge(width = 1))


lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01 * 50, 0),
         oChlaEpi_use = oChlaEpi/200)  |> 
  mutate(ID2 = ifelse(ID < 13, 'p1', 'p2')) |> 
  mutate(.by = any_of(c("ID2","year")), ID = row_number()) |> # renumber the pathways  
  select(all_of(c('ID', 'ID2','year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID','ID2', 'year'))) |> 
  ggplot(aes(y=year, x = name, 
             linewidth = value, 
             colour = name,
             group = name)) +
  geom_line(position = position_dodge(width = 1)) +
  facet_grid(ID~ID2) +
  coord_flip()


lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01 * 50, 0),
         oChlaEpi_use = oChlaEpi/200)  |> 
  select(all_of(c('ID', 'year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  ggplot(aes(y=year, x = name, 
             size = value, 
             colour = name,
             group = name)) +
  geom_line() +
  facet_wrap(~ID) +
  coord_flip()


lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi/0.01 * 50, 0),
         oChlaEpi_use = oChlaEpi/200)  |> 
  select(all_of(c('ID', 'year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  # filter(ID %in% 1:6) |> 
  ggplot(aes(x=year, y = name, 
             size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1) +
  # ggh4x::facet_nested(vars(ID,name),  scales = 'free_y', axes = 'y', switch = 'y') +
  ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y', axes = 'y', nrow = 18, 
                           strip.position = 'left', dir = 'v', remove_labels = 'all',
                           nest_line = element_line(colour = 'black'), 
                           strip = strip_nested(text_y = list(element_text(), 
                                                              element_text(colour = 'white')),
                                                background_y = list(element_rect(),element_blank()), 
                                                by_layer_y = TRUE)) +
  theme_bw() +
  theme(axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.border = element_rect(colour = 'black'))


example_pthway_plot <- lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  select(ID) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.01),
         oChlaEpi_use = oChlaEpi)  |> 
  select(all_of(c('ID', 'year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('ID','year'))) |> 
  # filter(ID %in% 7:9) |>
  ggplot(aes(x=year, y = value, 
             # size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  # ggh4x::facet_nested(vars(ID,name),  scales = 'free_y', axes = 'y', switch = 'y') +
  ggh4x::facet_nested_wrap(vars(ID, name), scales = 'free_y',  nrow = 18, 
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'), 
                           strip = strip_nested(text_y = list(element_text(), 
                                                              element_text(colour = 'white', 
                                                                           size = 1)),
                                                background_y = list(element_rect(),
                                                                    element_blank()), 
                                                by_layer_y = TRUE)) +
  theme_bw(base_size = 12) +
  scale_y_continuous(n.breaks = 3) +
  scale_colour_discrete(name = '') +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top') +
  facetted_pos_scales(y = list(name == "fMarsh_use" ~ scale_y_continuous(limits = c(0,1),
                                                                         n.breaks = 2),
                               name == "mPLoadEpi_use" ~ scale_y_continuous(limits = c(0,0.01),
                                                                            n.breaks = 2),
                               name == "oChlaEpi_use" ~ scale_y_continuous(limits = c(0,200),
                                                                           n.breaks = 2)))

ggsave(example_pthway_plot, filename = 'output/plots/example_pathways.png',
       height = 20, width = 30, units = 'cm')


## select a few extremes -------
extreme_pathways <- lastpoppathways |> 
  filter(year == max(lastpoppathways$year),
         oChlaEpi <= 20) |> 
  pivot_longer(cols = all_of(names(summary_1$bestmem)), names_to = 'measure') |> 
  arrange(measure, value) |> # arrange in order of value
  slice(.by = measure, c(1, n())) |>  # take first and last in each group
  mutate(labels = ifelse(value == min(value), 
                         paste0('Minimum ', measure),
                         paste0('Maximum ', measure)),
         .by = measure) |> 
  select(ID, labels) |> 
  left_join(lastpoppathways, by = join_by(ID)) |> 
  mutate(.by = year, ID = row_number()) |> # renumber the pathways 
  mutate(fMarsh_use = ifelse(fMarsh_lag < year, fMarsh, 0),
         mPLoadEpi_use = ifelse(mPLoadEpi_lag < year, mPLoadEpi, 0.01),
         oChlaEpi_use = oChlaEpi)  |> 
  select(all_of(c('labels', 'year', 'oChlaEpi_use', 'mPLoadEpi_use', 'fMarsh_use'))) |> 
  pivot_longer(cols = !any_of(c('labels','year'))) |> 
  ggplot(aes(x=year, y = value, 
             # size = value, 
             colour = name)) +
  geom_line(lineend = 'round', linejoin = 'round', linemitre = 1, linewidth = 1) +
  # ggh4x::facet_nested(vars(ID,name),  scales = 'free_y', axes = 'y', switch = 'y') +
  ggh4x::facet_nested_wrap(vars(labels, name), scales = 'free_y',  nrow = 6, 
                           strip.position = 'right', dir = 'v', remove_labels = 'y',
                           nest_line = element_line(colour = 'black'), 
                           strip = strip_nested(text_y = list(element_text(), 
                                                              element_text(colour = 'white', 
                                                                           size = 1)),
                                                background_y = list(element_rect(),
                                                                    element_blank()), 
                                                by_layer_y = TRUE)) +
  theme_bw(base_size = 12) +
  scale_y_continuous(n.breaks = 3) +
  scale_colour_discrete(name = '') +
  theme(panel.border = element_rect(colour = 'black'),
        legend.position = 'top') +
  facetted_pos_scales(y = list(name == "fMarsh_use" ~ scale_y_continuous(limits = c(0,1),
                                                                         n.breaks = 2),
                               name == "mPloadEpi_use" ~ scale_y_continuous(limits = c(0,0.01),
                                                                            n.breaks = 2),
                               name == "oChlaEpi_use" ~ scale_y_continuous(limits = c(0,200),
                                                                           n.breaks = 2)))

ggsave(extreme_pathways, filename = 'output/plots/extreme_pathways.png',
       height = 12, width = 30, units = 'cm')
