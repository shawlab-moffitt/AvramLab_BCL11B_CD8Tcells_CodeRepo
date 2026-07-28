library(car)
library(emmeans)
library(effectsize)
library(ggplot2)

# Create data frame
data_long <- data.frame(
  Response = c(
    # WT Day 8
    448.3, 449.4, 499.1, 583.6,
    
    # KO Day 8
    80.1, 53.4, 95.3, 110.2,
    
    # WT Day 10
    628.1, 752.5, 622.2, 630.2, 605.4, 635, 476.3,
    
    # KO Day 10
    286.2, 290.5, 283.1, 250.6, 248.6, 213.9, 278,
    
    # WT Day 15
    589.7, 664.2, 1020, 834.8, 998.6, 727, 740.7,
    
    # KO Day 15
    259.5, 304.1, 342.3, 276.9, 411.1, 289.7, 395.7
  ),
  
  Time = factor(c(
    rep(8, 4), rep(8, 4),
    rep(10, 7), rep(10, 7),
    rep(15, 7), rep(15, 7)
  )),
  
  Group = factor(c(
    rep("WT", 4), rep("KO", 4),
    rep("WT", 7), rep("KO", 7),
    rep("WT", 7), rep("KO", 7)
  ))
)

# Verify sample sizes
table(data_long$Time, data_long$Group)

# Fit two-way ANOVA
anova_model <- lm(Response ~ Group * Time, data = data_long)

# Type II ANOVA
Anova(anova_model, type = 2)

# Type III ANOVA (optional)
Anova(anova_model, type = 3)

# Standard ANOVA table
anova(anova_model)

# Estimated marginal means
emm <- emmeans(anova_model, ~ Group | Time)

# WT vs KO at each time point
pairs(emm, adjust = "holm")

# Compare time points within each genotype
emm_time <- emmeans(anova_model, ~ Time | Group)
pairs(emm_time, adjust = "tukey")

# Model diagnostics
par(mfrow = c(2, 2))
plot(anova_model)

# Normality
shapiro.test(residuals(anova_model))

# Equal variance
leveneTest(Response ~ Group * Time, data = data_long)

# Interaction plot
ggplot(data_long,
       aes(Time, Response,
           color = Group,
           group = Group)) +
  stat_summary(fun = mean,
               geom = "line",
               linewidth = 1) +
  stat_summary(fun = mean,
               geom = "point",
               size = 3) +
  stat_summary(fun.data = mean_se,
               geom = "errorbar",
               width = 0.15) +
  theme_classic(base_size = 14) +
  labs(
    x = "Time",
    y = "Response"
  )
