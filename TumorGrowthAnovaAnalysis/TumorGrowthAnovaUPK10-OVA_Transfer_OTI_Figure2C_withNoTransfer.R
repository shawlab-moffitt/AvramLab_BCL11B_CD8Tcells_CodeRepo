############################################################
# Two-way ANOVA:
# Group: NO TRANSFER, WT, KO
# Time: 8, 10, 15
#
# Model:
# Response ~ Group * Time
############################################################


############################################################
# 1. Install required packages
############################################################

required_packages <- c(
  "dplyr",
  "tidyr",
  "ggplot2",
  "car",
  "emmeans"
)

packages_to_install <- required_packages[
  !required_packages %in% rownames(installed.packages())
]

if (length(packages_to_install) > 0) {
  install.packages(packages_to_install)
}


############################################################
# 2. Load packages
############################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(car)
library(emmeans)

############################################################
# 4. Convert from wide format to long format
############################################################

# Create data frame
data_long <- data.frame(
  Response = c(
    # NO TRANSFER
    509.7, 637.9, 474.1, 473.2,
    548.3, 649.4, 499.1, 583.6, 611.3,
    900, 925, 700, 950,
    
    # WT
    426, 410, 391, 373.8,
    519, 419, 397, 293, 504,
    410, 591.1, 573.8, 633.3,
    950, 998, 1000, 1100,
    
    # KO
    240, 291.7, 0, 127, 299.5, 267, 226.5, 253.3,
    247, 259, 178.9, 245, 119,
    299.5, 267, 226.5, 253.3,
    328.1, 339.9, 358.5, 400
  ),
  
  Time = factor(c(
    # NO TRANSFER
    rep(8, 4),
    rep(10, 5),
    rep(15, 4),
    
    # WT
    rep(8, 4),
    rep(10, 5),
    rep(15, 4),
    rep(20, 4),
    
    # KO
    rep(8, 8),
    rep(10, 5),
    rep(15, 4),
    rep(20, 4)
  )),
  
  Group = factor(c(
    # NO TRANSFER
    rep("NO TRANSFER", 13),
    
    # WT
    rep("WT", 17),
    
    # KO
    rep("KO", 21)
  ), levels = c("NO TRANSFER", "WT", "KO"))
)


aov1=aov(Response ~ Group+Time+Group * Time, data = data_long)
summary(aov1)
TukeyHSD(aov1,'Group:Time')

# ---- exclude the following analysis ----
# Inspect number of observations per group and time
print(
  dplyr::count(data_long, Time, Group)
)

summary_statistics <- data_long %>%
  group_by(Time, Group) %>%
  summarise(
    N = n(),
    Mean = mean(Response, na.rm = TRUE),
    SD = sd(Response, na.rm = TRUE),
    SEM = SD / sqrt(N),
    Median = median(Response, na.rm = TRUE),
    Minimum = min(Response, na.rm = TRUE),
    Maximum = max(Response, na.rm = TRUE),
    .groups = "drop"
  )

print(summary_statistics)


raw_data_plot <- ggplot(
  data_long,
  aes(
    x = Time,
    y = Response,
    color = Group,
    group = Group
  )
) +
  geom_jitter(
    aes(shape = Group),
    width = 0.08,
    size = 3,
    alpha = 0.8
  ) +
  stat_summary(
    fun = mean,
    geom = "point",
    size = 4,
    position = position_dodge(width = 0.15)
  ) +
  stat_summary(
    fun = mean,
    geom = "line",
    linewidth = 1,
    position = position_dodge(width = 0.15)
  ) +
  theme_classic() +
  labs(
    title = "Response by Group and Time",
    x = "Time",
    y = "Response",
    color = "Group",
    shape = "Group"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(raw_data_plot)

############################################################
# 8. Plot mean +/- SEM
############################################################

summary_plot <- ggplot(
  summary_statistics,
  aes(
    x = Time,
    y = Mean,
    color = Group,
    group = Group
  )
) +
  geom_line(linewidth = 1) +
  geom_point(size = 4) +
  geom_errorbar(
    aes(
      ymin = Mean - SEM,
      ymax = Mean + SEM
    ),
    width = 0.15,
    linewidth = 0.8
  ) +
  theme_classic() +
  labs(
    title = "Mean Response by Group and Time",
    subtitle = "Error bars represent SEM",
    x = "Time",
    y = "Mean response",
    color = "Group"
  ) +
  theme(
    plot.title = element_text(
      hjust = 0.5,
      face = "bold"
    ),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 12),
    legend.title = element_text(face = "bold")
  )

print(summary_plot)

anova_model <- lm(
  Response ~ Group * Time,
  data = data_long
)

# Standard model summary
print(summary(anova_model))

# Standard sequential ANOVA table
print(anova(anova_model))
