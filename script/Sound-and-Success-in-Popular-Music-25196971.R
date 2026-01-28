# Module: IJC437 Coursework (2025–26)
# Student Registration Number: 25196971
# MSc Data Science
# University of Sheffield

#Project Title: Sound and Success in Popular Music

# Research Questions:
#   EDA : How have average acoustic characteristics changed over time (2000–2018)?
#   RQ1 : Which acoustic features changed most between 2000-2018?
#   RQ2 : Which acoustic features are linked to being a hit (top 25%)?
#   RQ3 : How do acoustic feature distributions differ between early (2000–2003) and late (2015–2018) eras?

# Dataset: MusicOSet (2000–2018)  |  Files (TSV):
# - song_chart.csv
# - song_pop.csv
# - acoustic_features.csv


graphics.off()
rm(list = ls())

# Packages
pkgs <- c("tidyverse", "lubridate", "broom", "ggrepel", "patchwork",
          "DiagrammeR", "DiagrammeRsvg", "rsvg")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install) > 0) install.packages(to_install)

library(tidyverse)
library(lubridate)
library(broom)
library(ggrepel)
library(patchwork)
library(DiagrammeR)

set.seed(2026)

# Folders
dir.create("plots_final", showWarnings = FALSE)
dir.create("outputs", showWarnings = FALSE)



run_info <- tibble(
  run_time = as.character(Sys.time()),
  r_version    = paste(R.Version()$major, R.Version()$minor, sep = "."),
  seed         = 2026
)
message("Run information")
print(run_info)



# Methodology diagram
g <- DiagrammeR::grViz("
digraph methodology {

  graph [layout = dot, rankdir = TB]
  node [shape = box, style = filled, fillcolor = '#E8EEF7']

  A [label = 'MusicOSet TSV files\\n(song_chart, song_pop,\\nacoustic_features)']
  B [label = 'Load data']
  C [label = 'Derive Debut_Year\\n(first chart appearance)']
  D [label = 'Inner join tables\\n(by song_id)']
  E [label = 'Clean data\\n(select vars, drop NA)']
  F [label = 'Min–max scaling\\n(acoustic features)']
  G [label = 'EDA:\\nYearly means (2000–2018)']
  H [label = 'RQ1:\\nLinear trends per feature']
  I [label = 'RQ2:\\nLogistic regression\\n(Hit vs Non-Hit)']
  J [label = 'RQ3:\\nEra comparison\\n(Early vs Late)']
  K [label = 'Outputs:\\nFigures + tables']

  A -> B -> C -> D -> E -> F -> G
  G -> H
  G -> I
  G -> J
  H -> K
  I -> K
  J -> K
}
")

# Show in Viewer (htmlwidget)
print(g)

# Save as PNG + open file (system viewer)
if (requireNamespace("DiagrammeRsvg", quietly = TRUE) &&
    requireNamespace("rsvg", quietly = TRUE)) {
  
  svg_txt <- DiagrammeRsvg::export_svg(g)
  
  rsvg::rsvg_png(
    charToRaw(svg_txt),
    file  = "plots_final/Methodology_Flowchart.png",
    width = 1400
  )
  
  message("Saved methodology diagram to plots_final/Methodology_Flowchart.png")
  
  if (interactive()) {
    browseURL(normalizePath("plots_final/Methodology_Flowchart.png"))
  }
  
} else {
  message("Diagram shown in Viewer; PNG not saved.")
}



# Loader (TSV)
load_data <- function(path) readr::read_delim(path, delim = "\t", show_col_types = FALSE)

# 1) Load files
df_chart <- load_data("data_raw/song_chart.csv") %>% mutate(week = as.Date(week))
df_pop   <- load_data("data_raw/song_pop.csv")
df_audio <- load_data("data_raw/acoustic_features.csv")

# Data audit
audit_raw <- list(
  song_chart_rows = nrow(df_chart),
  song_pop_rows   = nrow(df_pop),
  audio_rows      = nrow(df_audio),
  chart_na_week   = sum(is.na(df_chart$week))
)

# NA counts per column
na_chart <- sort(colSums(is.na(df_chart)), decreasing = TRUE)
na_pop   <- sort(colSums(is.na(df_pop)), decreasing = TRUE)
na_audio <- sort(colSums(is.na(df_audio)), decreasing = TRUE)

message("Raw data row counts")
print(enframe(audit_raw))

message("NA counts per column: song_chart")
print(enframe(na_chart, name = "column", value = "na_count") %>% filter(na_count > 0))

message("NA counts per column: song_pop")
print(enframe(na_pop, name = "column", value = "na_count") %>% filter(na_count > 0))

message("NA counts per column: acoustic_features")
print(enframe(na_audio, name = "column", value = "na_count") %>% filter(na_count > 0))

# 2) Create one Debut_Year per song using its first chart appearance (prevents double-counting songs across years)
df_debuts <- df_chart %>%
  filter(!is.na(week)) %>%
  group_by(song_id) %>%
  arrange(week) %>%
  slice(1) %>%
  ungroup() %>%
  transmute(song_id, Debut_Year = year(week)) %>%
  filter(between(Debut_Year, 2000, 2018))

message("Unique songs with Debut_Year (2000–2018): ", nrow(df_debuts))
message("Songs debuting per year (audit)")
print(
  df_debuts %>%
    count(Debut_Year, name = "songs_debuting"),
  n = Inf
)

# Join audit
n_audio <- df_audio %>% distinct(song_id) %>% nrow()
n_pop   <- df_pop   %>% distinct(song_id) %>% nrow()
n_debut <- df_debuts %>% distinct(song_id) %>% nrow()

message("Unique song_id counts | audio: ", n_audio,
        " | pop: ", n_pop,
        " | debuts: ", n_debut)


# 3) Join table and clean
df_raw <- df_audio %>%
  inner_join(df_pop, by = "song_id") %>%
  inner_join(df_debuts, by = "song_id") %>%
  select(
    song_id, Debut_Year,
    year_end_score,
    danceability, energy, valence, acousticness, loudness, tempo
  ) %>%
  drop_na()


# Cleaning audit
df_pre_na <- df_audio %>%
  inner_join(df_pop, by = "song_id") %>%
  inner_join(df_debuts, by = "song_id") %>%
  select(song_id, Debut_Year, year_end_score,
         danceability, energy, valence, acousticness, loudness, tempo)

dropped_rows <- nrow(df_pre_na) - nrow(df_raw)

# Summary audit
message("Rows after join (before drop_na): ", nrow(df_pre_na))
message("Rows after drop_na: ", nrow(df_raw))
message("Rows dropped due to missing values: ", dropped_rows)

# NA counts per column before drop_na
message("NA counts per column (after join, before drop_na)")
print(
  enframe(sort(colSums(is.na(df_pre_na)), decreasing = TRUE),
          name = "column", value = "na_count") %>%
    filter(na_count > 0),
  n = Inf
)



# Min–max scaling (0–1) so features are comparable on the same scale.
# makes interpretation relative (higher vs lower), not in original units such as dB (loudness) or BPM (tempo).
norm_safe <- function(x) {
  r <- max(x, na.rm = TRUE) - min(x, na.rm = TRUE)
  if (is.na(r) || r == 0) return(rep(0.5, length(x)))
  (x - min(x, na.rm = TRUE)) / r
}

df <- df_raw %>%
  mutate(across(danceability:tempo, norm_safe))

message("Years covered: ", min(df$Debut_Year), "–", max(df$Debut_Year))
message("Final analysis rows: ", nrow(df))
message("Cleaned and scaled data:")
print(head(df, 6))


# DATA DICTIONARY
# Describes variables used in the final analysis dataset (df)
data_dictionary <- tibble(
  variable = c(
    "song_id",
    "Debut_Year",
    "year_end_score",
    "danceability",
    "energy",
    "valence",
    "acousticness",
    "loudness",
    "tempo"
  ),
  description = c(
    "Unique song identifier used to join datasets",
    "Year of first appearance in the chart (earliest chart week)",
    "Year-end popularity score (higher values indicate greater chart success)",
    "Danceability of the track (relative, scaled 0–1)",
    "Perceived intensity and activity of the track (relative, scaled 0–1)",
    "Musical positiveness conveyed by the track (relative, scaled 0–1)",
    "Degree to which the track is acoustic (relative, scaled 0–1)",
    "Overall loudness of the track (originally in dB, now relative 0–1)",
    "Estimated tempo of the track (originally BPM, now relative 0–1)"
  ),
  source = c(
    "All tables",
    "Derived from song_chart",
    "song_pop",
    "acoustic_features",
    "acoustic_features",
    "acoustic_features",
    "acoustic_features",
    "acoustic_features",
    "acoustic_features"
  ),
  notes = c(
    "Primary join key",
    "Derived to avoid counting the same song multiple times across years",
    "Used to define a 'Hit' as top 25%",
    "Min–max scaled across full dataset",
    "Min–max scaled across full dataset",
    "Min–max scaled across full dataset",
    "Min–max scaled across full dataset",
    "Min–max scaled across full dataset; unitless after scaling",
    "Min–max scaled across full dataset; unitless after scaling"
  )
)

message("DATA DICTIONARY")
print(data_dictionary, n = Inf)




feature_colors <- c(
  acousticness  = "#1b9e77",
  energy        = "#66a61e",
  tempo         = "#7570b3",
  danceability  = "#e6ab02",
  loudness      = "#a6761d",
  valence       = "#e7298a"
)

feature_labels <- c(
  acousticness = "Acousticness",
  danceability = "Danceability",
  energy       = "Energy",
  loudness     = "Loudness",
  tempo        = "Tempo",
  valence      = "Valence"
)

feature_order <- names(feature_colors)

# EDA: How have average acoustic characteristics changed over time (2000–2018)?
eda_yearly <- df %>%
  group_by(Debut_Year) %>%
  summarise(across(danceability:tempo, mean), .groups = "drop") %>%
  pivot_longer(-Debut_Year, names_to = "Feature", values_to = "MeanValue") %>%
  mutate(
    Feature = factor(Feature, levels = feature_order),
    FeatureLabel = recode(Feature, !!!feature_labels)
  )

p_eda1 <- ggplot(eda_yearly, aes(x = Debut_Year, y = MeanValue, color = Feature)) +
  geom_line(linewidth = 1.15, alpha = 0.95) +
  geom_point(size = 1.9) +
  geom_text_repel(
    data = eda_yearly %>% filter(Debut_Year == 2018),
    aes(label = FeatureLabel),
    nudge_x = 0.6,
    segment.color = NA,
    size = 4
  ) +
  scale_color_manual(values = feature_colors, breaks = feature_order, labels = feature_labels) +
  scale_x_continuous(breaks = seq(2000, 2018, by = 2), limits = c(2000, 2019)) +
  scale_y_continuous(limits = c(0, 1), labels = scales::number_format(accuracy = 0.1)) +
  labs(
    title = "EDA1. How have the average acoustic characteristics changed over time (2000–2018)?",
    subtitle = "Each line shows the mean feature value for songs debuting in that year (scaled 0–1).",
    x = "Song debut year (first appearance in the chart)",
    y = "Mean scaled feature value (0 = low, 1 = high)",
    color = "Feature"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave("plots_final/Fig1_EDA1_Acoustic_Trends_2000_2018.png", p_eda1, width = 12, height = 7, dpi = 200)
print(p_eda1)


# RQ1 : Which acoustic features changed most between 2000-2018?
# Linear regression
# Uses yearly means and fits MeanValue ~ Debut_Year per feature.
line_data <- eda_yearly %>% transmute(Debut_Year, Feature, Value = MeanValue)

rq1_trends <- line_data %>%
  group_by(Feature) %>%
  do({
    model <- lm(Value ~ Debut_Year, data = .)
    tibble(
      slope_per_year = unname(coef(model)[2]),
      change_18yrs   = unname(coef(model)[2]) * (2018 - 2000)
    )
  }) %>%
  ungroup() %>%
  mutate(
    Feature = factor(Feature, levels = feature_order),
    FeatureLabel = recode(Feature, !!!feature_labels),
    Direction = ifelse(change_18yrs > 0, "Increased", "Decreased")
  )

write_csv(rq1_trends, "outputs/RQ1_linear_total_change_2000_2018.csv")
message("RQ1 summary")
print(rq1_trends)

max_abs <- max(abs(rq1_trends$change_18yrs), na.rm = TRUE)
p_rq1 <- ggplot(
  rq1_trends,
  aes(x = reorder(FeatureLabel, change_18yrs), y = change_18yrs, fill = Feature)
) +
  geom_col(width = 0.72) +
  coord_flip() +
  scale_fill_manual(values = feature_colors, guide = "none") +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey35") +
  scale_y_continuous(
    limits = c(-max_abs, max_abs),
    breaks = scales::pretty_breaks(n = 5),
    labels = scales::number_format(accuracy = 0.01)
  ) +
  labs(
    title = "RQ1. Which acoustic features changed most between 2000 and 2018?",
    subtitle = "Total change estimated from linear trends in yearly average feature values.",
    x = NULL,
    y = "Estimated change in scaled yearly mean (2018–2000)"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave("plots_final/Fig2_RQ1_TotalChange_2000_2018.png", p_rq1, width = 10, height = 6, dpi = 200)
print(p_rq1)


# RQ2 : Which acoustic features are linked to being a hit (top 25%)?
# Logistic regression — Hit vs Non-Hit (Top 25% year_end_score)
# Forest plot (point + 95% CI)
hit_thresh <- quantile(df$year_end_score, 0.75, na.rm = TRUE)

df_rq2 <- df %>%
  mutate(
    is_hit = factor(ifelse(year_end_score >= hit_thresh, "Hit", "Non-Hit"),
                    levels = c("Non-Hit", "Hit"))
  )

logit_model <- glm(
  is_hit ~ danceability + energy + valence + acousticness + loudness + tempo,
  data = df_rq2,
  family = binomial()
)

rq2_or <- broom::tidy(logit_model) %>%
  filter(term != "(Intercept)") %>%
  mutate(
    Feature = term,
    FeatureLabel = recode(Feature, !!!feature_labels),
    OR      = exp(estimate),
    low_95  = exp(estimate - 1.96 * std.error),
    high_95 = exp(estimate + 1.96 * std.error)
  ) %>%
  select(Feature, FeatureLabel, OR, low_95, high_95, p.value) %>%
  arrange(desc(OR))

write_csv(rq2_or, "outputs/RQ2_logistic_odds_ratios.csv")
message("RQ2 summary")
print(rq2_or)

p_rq2 <- ggplot(
  rq2_or,
  aes(x = OR, y = fct_reorder(FeatureLabel, OR), color = Feature)
) +
  geom_vline(xintercept = 1, linetype = "dashed", colour = "grey35") +
  geom_errorbarh(aes(xmin = low_95, xmax = high_95), height = 0.22, linewidth = 0.9) +
  geom_point(size = 3.2) +
  scale_color_manual(values = feature_colors, guide = "none") +
  scale_x_log10(
    breaks = c(0.25, 0.5, 1, 2, 4, 8),
    labels = c("0.25", "0.5", "1", "2", "4", "8")
  ) +
  labs(
    title = "RQ2. Which acoustic features are linked to being a hit (top 25%)?",
    subtitle = "Logistic regression odds ratios (95% confidence intervals).",
    x = "Odds ratio (log scale) \n>1 higher hit odds; <1 lower hit odds.",
    y = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 11),
    axis.title.x = element_text(size = 12),
    panel.grid.minor = element_blank()
  )

ggsave("plots_final/Fig3_RQ2_ForestPlot_OddsRatios.png", p_rq2, width = 10, height = 6, dpi = 200)
print(p_rq2)


# RQ3 : How do acoustic feature distributions differ between early (2000–2003) and late (2015–2018) eras?
# Era comparison — distributions early vs late
df_rq3 <- df %>%
  mutate(
    Era = case_when(
      between(Debut_Year, 2000, 2003) ~ "Early era (2000–2003)",
      between(Debut_Year, 2015, 2018) ~ "Late era (2015–2018)",
      TRUE ~ NA_character_
    )
  ) %>%
  filter(!is.na(Era)) %>%
  mutate(Era = factor(Era, levels = c("Early era (2000–2003)", "Late era (2015–2018)")))

rq3_long <- df_rq3 %>%
  select(Era, danceability, energy, valence, acousticness, loudness, tempo) %>%
  pivot_longer(-Era, names_to = "Feature", values_to = "Value") %>%
  mutate(
    Feature = factor(Feature, levels = feature_order),
    FeatureLabel = recode(Feature, !!!feature_labels)
  )


p_rq3 <- ggplot(
  rq3_long,
  aes(x = Era, y = Value, fill = Era)
) +
  geom_violin(
    trim = TRUE,
    alpha = 0.5,
    linewidth = 0.35
  ) +
  geom_boxplot(
    width = 0.18,
    outlier.alpha = 0.1,
    linewidth = 0.4,
    fill = "white"
  ) +
  facet_wrap(~FeatureLabel, ncol = 3) +
  
  scale_y_continuous(
    limits = c(0, 1),
    breaks = seq(0, 1, by = 0.25),
    labels = scales::number_format(accuracy = 0.1)
  ) +
  
  scale_fill_manual(
    values = c(
      "Early era (2000–2003)" = "#4C72B0",
      "Late era (2015–2018)"  = "#DD8452"
    )
  ) +
  
  labs(
    title = "RQ3. How do acoustic feature distributions differ between two eras?",
    subtitle = "Comparison of early 2000s (2000–2003) and late 2010s (2015–2018) chart music.",
    x = NULL,
    y = "Normalised feature value (0–1)",
    fill = "Era"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 12),
    plot.subtitle = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      size = 10,
      angle = 10,
      hjust = 0.5,
      face = "bold"
    ),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

ggsave("plots_final/Fig4_RQ3_Era_Distributions.png", p_rq3, width = 12, height = 7, dpi = 200)
print(p_rq3)