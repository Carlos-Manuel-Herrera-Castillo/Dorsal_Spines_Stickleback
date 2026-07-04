#!/usr/bin/env Rscript

# Step 08d: Plot developmental sequence data.
# Original section: V* Plotting of developmental sequence
#
# Run with:
#   Rscript scripts/08_phenotypes/development_sequence.R

library(tidyverse)
library(patchwork)

data <- read.csv("/my/folder.path/DS_Developmental_graph_data.csv")
names(data) <- str_replace_all(names(data), " ", "_")
names(data) <- str_trim(names(data))

legend_levels <- c("Nothing", "Cartilage only", "Dorsal plate ossification", "Dorsal spine")

prepare_bar_data <- function(df, prefix) {
  df %>%
    select(Age, all_of(c(
      paste0(prefix, "_Nothing"),
      paste0(prefix, "_Cartilage_only"),
      paste0(prefix, "_Dorsal_plate_ossification")
    ))) %>%
    pivot_longer(cols = -Age, names_to = "Category", values_to = "Count") %>%
    mutate(
      Category = str_remove(Category, paste0(prefix, "_")),
      Category = str_replace_all(Category, "_", " "),
      Category = factor(Category, levels = legend_levels)
    ) %>%
    group_by(Age) %>%
    mutate(Proportion = Count / sum(Count)) %>%
    ungroup()
}

prepare_spine_data <- function(df, prefix) {
  spine_col <- paste0(prefix, "_Dorsal_spine")
  total_cols <- df %>%
    select(starts_with(prefix)) %>%
    select(-all_of(spine_col)) %>%
    names()

  df %>%
    select(Age, spine = all_of(spine_col)) %>%
    left_join(
      df %>% transmute(Age, Total = rowSums(select(., all_of(total_cols)))),
      by = "Age"
    ) %>%
    distinct() %>%
    mutate(Proportion = spine / Total)
}

duin_bar <- prepare_bar_data(data, "DUIN")
scad_bar <- prepare_bar_data(data, "SCAD")

duin_spine <- prepare_spine_data(data, "DUIN")
scad_spine <- prepare_spine_data(data, "SCAD")

custom_palette <- c(
  "Nothing" = "lightgray",
  "Cartilage only" = "#1E88E5",
  "Dorsal plate ossification" = "#D81B60",
  "Dorsal spine" = "#004D40"
)

dummy_legend <- data.frame(
  Age = 0,
  Proportion = 0,
  Category = factor("Dorsal spine", levels = legend_levels)
)

duin_plot <- ggplot() +
  geom_bar(
    data = duin_bar,
    aes(x = as.factor(Age), y = Proportion, fill = Category),
    stat = "identity",
    color = "white",
    width = 0.8
  ) +
  geom_line(
    data = duin_spine,
    aes(x = as.factor(Age), y = Proportion, group = 1),
    color = "#004D40",
    linewidth = 1.2
  ) +
  geom_point(
    data = duin_spine,
    aes(x = as.factor(Age), y = Proportion),
    color = "#004D40",
    size = 2
  ) +
  geom_point(
    data = dummy_legend,
    aes(x = Age, y = Proportion, fill = Category),
    color = "#004D40",
    shape = 21,
    size = 0,
    alpha = 0,
    show.legend = TRUE,
    inherit.aes = FALSE
  ) +
  scale_fill_manual(values = custom_palette, drop = FALSE) +
  labs(
    title = "Spined G. aculeatus - Pelvic Developmental Stages",
    y = "Proportion",
    x = "Age (days post hatching)"
  ) +
  theme_light(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    legend.title = element_blank(),
    plot.title = element_text(face = "bold", size = 15)
  )

scad_plot <- ggplot() +
  geom_bar(
    data = scad_bar,
    aes(x = as.factor(Age), y = Proportion, fill = Category),
    stat = "identity",
    color = "white",
    width = 0.8
  ) +
  geom_line(
    data = scad_spine,
    aes(x = as.factor(Age), y = Proportion, group = 1),
    color = "#004D40",
    linewidth = 1.2
  ) +
  geom_point(
    data = scad_spine,
    aes(x = as.factor(Age), y = Proportion),
    color = "#004D40",
    size = 2
  ) +
  scale_fill_manual(values = custom_palette, drop = FALSE) +
  labs(
    title = "Spineless G. aculeatus - Pelvic Developmental Stages",
    y = "Proportion",
    x = "Age (days post hatching)"
  ) +
  theme_light(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "none",
    plot.title = element_text(face = "bold", size = 15)
  )

final_plot <- duin_plot / scad_plot +
  plot_annotation(
    theme = theme(plot.title = element_text(size = 18, face = "bold", hjust = 0.5))
  )

final_plot
