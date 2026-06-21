##################################################
# PHASE 1: LOAD PACKAGES
##################################################

library(tidyverse)
library(skimr)
library(janitor)
library(gtsummary)
library(corrplot)
library(caret)
library(pROC)
library(car)
library(broom)
library(gt)


##################################################
# PHASE 2: DATA INSPECTION
##################################################

dim(heart)
names(heart)
str(heart)
summary(heart)
skim(heart)


##################################################
# TABLE A1: MISSING VALUES
##################################################

missing_table <- data.frame(
  Variable = names(heart),
  Missing = colSums(is.na(heart)),
  Percent = round(
    colSums(is.na(heart))/nrow(heart)*100,
    2
  )
)

missing_table

write.csv(
  missing_table,
  "Table_A1_MissingValues.csv",
  row.names = FALSE
)

##################################################
# TABLE A2: DUPLICATES
##################################################

sum(duplicated(heart))


##################################################
# PHASE 3: DATA PREPARATION
##################################################

heart$Sex <- factor(heart$Sex)

heart$ChestPainType <- factor(
  heart$ChestPainType
)

heart$FastingBS <- factor(
  heart$FastingBS
)

heart$RestingECG <- factor(
  heart$RestingECG
)

heart$ExerciseAngina <- factor(
  heart$ExerciseAngina
)

heart$ST_Slope <- factor(
  heart$ST_Slope
)

heart$HeartDisease <- factor(
  heart$HeartDisease
)

str(heart)

##################################################
# FIGURE 1
# AGE DISTRIBUTION
##################################################

fig1 <- ggplot(
  heart,
  aes(x = Age)
) +
  geom_histogram(
    bins = 20
  ) +
  labs(
    title = "Distribution of Age"
  ) +
  theme_minimal()

fig1

ggsave(
  "Figure1_AgeDistribution.png",
  fig1,
  width = 7,
  height = 5,
  dpi = 300
)

##################################################
# FIGURE 2
##################################################

fig2 <- ggplot(
  heart,
  aes(x = Sex)
) +
  geom_bar() +
  labs(
    title = "Sex Distribution"
  ) +
  theme_minimal()

fig2

ggsave(
  "Figure2_SexDistribution.png",
  fig2,
  width = 7,
  height = 5,
  dpi = 300
)

##################################################
# FIGURE 2: SEX DISTRIBUTION
##################################################

fig2 <- ggplot(
  heart,
  aes(x = Sex)
) +
  geom_bar() +
  labs(
    title = "Distribution of Participants by Sex",
    x = "Sex",
    y = "Frequency"
  ) +
  theme_minimal()

fig2

ggsave(
  filename = "Figure2_SexDistribution.png",
  plot = fig2,
  width = 7,
  height = 5,
  dpi = 300
)



##################################################
# FIGURE 3
##################################################

fig3 <- ggplot(
  heart,
  aes(x = HeartDisease)
) +
  geom_bar() +
  labs(
    title = "Heart Disease Distribution"
  ) +
  theme_minimal()

fig3

ggsave(
  "Figure3_HeartDiseaseDistribution.png",
  fig3,
  width = 7,
  height = 5,
  dpi = 300
)

##################################################
# FIGURE 4
##################################################

fig4 <- ggplot(
  heart,
  aes(
    x = HeartDisease,
    y = Age
  )
) +
  geom_boxplot() +
  labs(
    title = "Age by Heart Disease Status"
  ) +
  theme_minimal()

fig4

ggsave(
  "Figure4_AgeByDisease.png",
  fig4,
  width = 7,
  height = 5,
  dpi = 300
)


##################################################
# FIGURE 5
##################################################

fig5 <- ggplot(
  heart,
  aes(x = Cholesterol)
) +
  geom_histogram(
    bins = 20
  ) +
  labs(
    title = "Cholesterol Distribution"
  ) +
  theme_minimal()

fig5

ggsave(
  "Figure5_CholesterolDistribution.png",
  fig5,
  width = 7,
  height = 5,
  dpi = 300
)


##################################################
# FIGURE 6
##################################################

numeric_data <- heart %>%
  select(
    Age,
    RestingBP,
    Cholesterol,
    MaxHR,
    Oldpeak
  )

cor_matrix <- cor(
  numeric_data
)

png(
  "Figure6_CorrelationHeatmap.png",
  width = 1000,
  height = 800
)

corrplot(
  cor_matrix,
  method = "color",
  type = "upper"
)

dev.off()


##################################################
# TABLE 1
##################################################

table1 <- heart %>%
  tbl_summary(
    by = HeartDisease,
    statistic = list(
      all_continuous() ~
        "{mean} ({sd})"
    )
  ) %>%
  add_p()

table1

as_gt(table1)


chisq.test(
  table(
    heart$Sex,
    heart$HeartDisease
  )
)

chisq.test(
  table(
    heart$ExerciseAngina,
    heart$HeartDisease
  )
)


chisq.test(table(heart$Sex, heart$HeartDisease))

chisq.test(table(heart$ExerciseAngina, heart$HeartDisease))

t.test(Age ~ HeartDisease, data = heart)

t.test(Cholesterol ~ HeartDisease, data = heart)



##################################################
# MULTIVARIABLE MODEL
##################################################

model <- glm(
  HeartDisease ~
    Age +
    Sex +
    ChestPainType +
    RestingBP +
    Cholesterol +
    FastingBS +
    RestingECG +
    MaxHR +
    ExerciseAngina +
    Oldpeak +
    ST_Slope,
  data = heart,
  family = binomial
)

summary(model)

##################################################
# TABLE 3
##################################################
exp(confint(model))


table3 <- tbl_regression(
  model,
  exponentiate = TRUE
)

table3

as_gt(table3)

##################################################
# TABLE 4
##################################################

vif_table <- data.frame(
  Variable = names(vif(model)),
  VIF = vif(model)
)

vif_table

write.csv(
  vif_table,
  "Table4_VIF.csv",
  row.names = FALSE
)

##################################################
# FIGURE 7
##################################################

or_data <- tidy(
  model,
  exponentiate = TRUE,
  conf.int = TRUE
)

fig7 <- ggplot(
  or_data[-1,],
  aes(
    x = reorder(
      term,
      estimate
    ),
    y = estimate
  )
) +
  geom_point() +
  geom_errorbar(
    aes(
      ymin = conf.low,
      ymax = conf.high
    ),
    width = 0.2
  ) +
  coord_flip() +
  theme_minimal() +
  labs(
    title = "Adjusted Odds Ratios"
  )

fig7

ggsave(
  "Figure7_ForestPlot.png",
  fig7,
  width = 8,
  height = 6,
  dpi = 300
)

##################################################
# TABLE 5
##################################################

pred <- predict(
  model,
  type = "response"
)

pred_class <- ifelse(
  pred > 0.5,
  1,
  0
)

cm <- confusionMatrix(
  factor(pred_class),
  heart$HeartDisease
)

cm



















