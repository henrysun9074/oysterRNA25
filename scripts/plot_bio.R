library(cowplot)
library(tidyverse)
library(ggpubr)
library(readxl)

growth <- read.csv("raw/env_bio/growth2025.csv")
growth <- growth %>%
  mutate(Site = case_when(
    Site == "Duke Aquafarm" ~ "DAF",
    Site == "Ward Creek" ~ "WC",
    Site == "CMAST" ~ "CMT",
    TRUE ~ Site
  ))

## mortality data
duml_mort <- read_excel("/work/hs325/bass25/raw/env_bio/mortality2025.xlsx", sheet = "Duke Aquafarm")
duml_mort$location <- "DAF"
wc_mort <- read_excel("/work/hs325/bass25/raw/env_bio/mortality2025.xlsx", sheet = "Ward Creek")
wc_mort$location <- "WC"
cmt_mort <- read_excel("/work/hs325/bass25/raw/env_bio/mortality2025.xlsx", sheet = "CMAST")
cmt_mort$location <- "CMT"
mort <- rbind(duml_mort, wc_mort, cmt_mort)

## use these colors from the previous plots for consistency
loc_colors <- c("CMT" = "#7FC97F", "DAF" = "#BEAED4", "WC" = "#FDC086")

## date windows: growth was sampled through August, mortality through October 1
start_date      <- as.Date("2025-05-01")
growth_end_date <- as.Date("2025-08-31")
mort_end_date   <- as.Date("2025-10-01")

# plot only 3n
growth <- growth %>% filter(Strain == "Triploid")
mort <- mort %>%
  mutate(Date = as.Date(Date)) %>%
  filter(Strain == "Triploid", Date >= start_date, Date <= mort_end_date)

growth <- growth %>%
  mutate(Date = as.Date(Date, format = "%m/%d/%y")) %>%
  filter(Date >= start_date, Date <= growth_end_date)

## shared styling
bold_theme <- theme(
  text         = element_text(face = "bold"),
  axis.title   = element_text(face = "bold"),
  axis.text    = element_text(face = "bold"),
  strip.text   = element_text(face = "bold", size = 16),  # facet titles in panel A
  legend.title = element_text(face = "bold", size = 16),
  legend.text  = element_text(face = "bold", size = 14),
  legend.key.size = unit(1, "cm"),
  plot.margin  = margin(t = 10, r = 10, b = 5, l = 5)
)

# month ticks
month_scale <- scale_x_date(date_labels = "%b", date_breaks = "1 month")

## shaded date windows for RNA sampling
shade_windows <- data.frame(
  xmin = as.Date(c("2025-06-25", "2025-07-09", "2025-07-23")),
  xmax = as.Date(c("2025-06-27", "2025-07-11", "2025-07-25"))
)
shade_layer <- geom_rect(data = shade_windows,
                         aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
                         inherit.aes = FALSE, fill = "gray50", alpha = 0.4)

## growth
growth_long <- growth %>%
  pivot_longer(c(L_mm, W_mm, D_mm), names_to = "Dimension", values_to = "mm") %>%
  mutate(Dimension = factor(Dimension, levels = c("L_mm", "W_mm", "D_mm"),
                            labels = c("Length (mm)", "Width (mm)", "Depth (mm)")))

growth_sum <- growth_long %>%
  group_by(Site, Date, Dimension) %>%
  summarise(mean = mean(mm, na.rm = TRUE),
            se = sd(mm, na.rm = TRUE) / sqrt(sum(!is.na(mm))),
            .groups = "drop")

p_growth <- ggplot(growth_sum, aes(Date, mean, color = Site, group = Site)) +
  shade_layer + 
  geom_jitter(data = growth_long, aes(Date, mm, color = Site), inherit.aes = FALSE,
              width = 2, height = 0, alpha = 0.5, size = 1) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 2) +
  facet_wrap(~ Dimension, scales = "free_y") +
  scale_color_manual(values = loc_colors) +
  month_scale +
  coord_cartesian(xlim = c(start_date, growth_end_date)) +
  labs(x = NULL, y = "Mean +/- SE", color = "Site") +
  theme_pubr() +
  theme(strip.background = element_blank()) +
  bold_theme

## mortality
mort_sum <- mort %>%
  group_by(location, Date) %>%
  summarise(mean = mean(Survivorship, na.rm = TRUE),
            se = sd(Survivorship, na.rm = TRUE) / sqrt(sum(!is.na(Survivorship))),
            .groups = "drop")

p_mort <- ggplot(mort_sum, aes(Date, mean, color = location, group = location)) +
  shade_layer + 
  geom_jitter(data = mort, aes(Date, Survivorship, color = location), inherit.aes = FALSE,
              width = 2, height = 0, alpha = 0.5, size = 1.5) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 2) +
  scale_color_manual(values = loc_colors) +
  month_scale +
  coord_cartesian(xlim = c(start_date, mort_end_date)) +
  labs(x = NULL, y = "Survivorship (%)", color = "Site") +
  theme_pubr() +
  bold_theme

## COMBINE
shared_legend <- get_legend(p_mort + theme(legend.position = "right",
                                           legend.box.margin = margin(0, 20, 0, 10)))
panels <- plot_grid(
  p_growth + theme(legend.position = "none"),
  p_mort   + theme(legend.position = "none"),
  ncol = 1, align = "v", axis = "l",
  labels = c("A", "B"), label_size = 22, label_fontface = "bold"
)

p_both <- plot_grid(panels, shared_legend, nrow = 1, rel_widths = c(1, 0.15))
p_both

## save
dir.create("raw/env_bio/plots", showWarnings = FALSE)
ggsave("raw/env_bio/plots/growth_mortality_triploid.png", p_both, width = 12, height = 8)