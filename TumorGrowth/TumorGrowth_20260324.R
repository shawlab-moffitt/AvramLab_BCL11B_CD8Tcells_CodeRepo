library(lme4)
library(lmerTest)
library(emmeans)
file = "PMEL_TumorGrowthData_20260330.tsv"
data <- read.delim(file, header = TRUE, sep = "\t")
data$Group=factor(data$Group,level=c('WT','NT','KO'))
model_pmel <- lmer(log2(TumorSize + 1) ~ Days + Group + Days * Group + (1  | ID), data = data)
summary(model_pmel)

data$Group=factor(data$Group,level=c('NT','WT','KO'))
model_pmel <- lmer(log2(TumorSize + 1) ~ Days + Group + Days * Group + (1  | ID), data = data)
summary(model_pmel)


emtrends(model_pmel, pairwise ~ Days * Group, var = "Days")

library(ggplot2)

aov_tab <- anova(model_pmel, type = 3)
p_interaction <- aov_tab["Days:Group", "Pr(>F)"]

ggplot(data, aes(x = Days, y = TumorSize, color = Group)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  annotate(
    "text",
    x = max(data$Days) * 0.7,
    y = max(data$TumorSize, na.rm = TRUE),
    label = paste0("Mixed-effects model\nInteraction p = ", signif(p_interaction, 3)),
    hjust = 0,
    size = 5
  ) +
  theme_classic()


file = "OVA_OTI_TumorGrowthData_20260330.tsv"
data_ova <- read.delim(file, header = TRUE, sep = "\t")
data_ova$Group=factor(data_ova$Group,level=c('WT','NT','KO'))
model_ova <- lmer(log2(TumorSize + 1) ~ Days + Group + Days * Group + (1  | ID), data = data_ova)
summary(model_ova)

emtrends(model_ova, pairwise ~ Group, var = "Days")

data_ova$Group=factor(data_ova$Group,level=c('NT','WT','KO'))
model_ova <- lmer(log2(TumorSize + 1) ~ Days + Group + Days * Group + (1  | ID), data = data_ova)
summary(model_ova)

#emtrends(model_ova, pairwise ~ Group, var = "Days")

emtrends(model_ova, pairwise ~ Days * Group, var = "Days")
