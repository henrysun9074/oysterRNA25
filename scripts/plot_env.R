library(cowplot)
library(tidyverse)
library(ggpubr)

## load data
env <- read.csv("raw/env_bio/env2025.csv") %>%
  select(-X) %>%
  mutate(DateTime = as.POSIXct(DateTime, format = "%Y-%m-%d %H:%M:%S", tz = "UTC"))

## (darkened ~20% from "#7FC97F", "#BEAED4", "#FDC086")
loc_colors <- c("CMT" = "#66A166", "DAF" = "#988BAA", "WC" = "#CA9A6B")

## date window
start_date <- as.POSIXct("2025-06-21", tz = "UTC")
end_date   <- as.POSIXct("2025-08-01", tz = "UTC")

env <- env %>%
  filter(DateTime >= start_date, DateTime <= end_date)

## shared styling
bold_theme <- theme(
  text         = element_text(face = "bold"),
  axis.title   = element_text(face = "bold"),
  axis.text    = element_text(face = "bold"),
  strip.text   = element_text(face = "bold", size = 16),
  legend.title = element_text(face = "bold", size = 16),
  legend.text  = element_text(face = "bold", size = 14),
  legend.key.size = unit(1, "cm"),
  plot.margin  = margin(t = 10, r = 10, b = 5, l = 5)
)

# week ticks
week_scale <- scale_x_datetime(date_labels = "%b %d", date_breaks = "1 week")

## shaded date windows for RNA sampling
shade_windows <- data.frame(
  xmin = as.POSIXct(c("2025-06-25", "2025-07-09", "2025-07-23"), tz = "UTC"),
  xmax = as.POSIXct(c("2025-06-27", "2025-07-11", "2025-07-25"), tz = "UTC")
)
shade_layer <- geom_rect(data = shade_windows,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "gray50", alpha = 0.4)

## environmental variables
env_long <- env %>%
  pivot_longer(c(Temp_C, Salinity_ppt, DO_mgL, pH), names_to = "Variable", values_to = "value") %>%
  mutate(Variable = factor(Variable, levels = c("Temp_C", "Salinity_ppt", "DO_mgL", "pH"),
                           labels = c("Temperature (°C)", "Salinity (ppt)", "DO (mg/L)", "pH")))

p_env <- ggplot(env_long, aes(DateTime, value, color = Location, group = Location)) +
  shade_layer +
  geom_line(linewidth = 0.4, alpha = 0.8, na.rm = TRUE) +
  facet_wrap(~ Variable, ncol = 1, scales = "free_y", strip.position = "left") +
  scale_color_manual(values = loc_colors) +
  week_scale +
  coord_cartesian(xlim = c(start_date, end_date)) +
  guides(color = guide_legend(override.aes = list(linewidth = 2, alpha = 1))) +
  labs(x = NULL, y = NULL, color = "Site") +
  theme_pubr(legend = "right") +
  theme(strip.background = element_blank(), strip.placement = "outside") +
  bold_theme
p_env

## save
dir.create("raw/env_bio/plots", showWarnings = FALSE)
ggsave("raw/env_bio/plots/env_jun21_aug1.png", p_env, width = 12, height = 8)
