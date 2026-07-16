##############################################################
# R code to generate Figure 4 in Koch et al. 2026
# Author: Abigail Meyer
# Last edit: 04/23/2026
##############################################################

library(ggplot2)
library(dplyr)
library(tidyr)
library(glmmTMB)
library(DHARMa)
library(patchwork)
library(knitr)

setwd("~/Desktop/Phd/Projects/Urban_Photobionts")


# 1. Load and prepare data

photo <- read.csv(
  "photobiont_correct.csv",
  stringsAsFactors = TRUE,
  header = TRUE,
  row.names = 1,
  na.strings = "NA"
)

physio <- read.csv(
  "physio_land_all_correct.csv",
  stringsAsFactors = TRUE,
  header = TRUE,
  row.names = 1,
  na.strings = "NA"
)

# Pull A46/I05 (dominant clades) counts from photobiont data
counts <- photo[, c("Sample_code", "A46", "I05")]

# Join physiological data (which contains site data) with count data 
merged <- left_join(physio, counts, by = "Sample_code")

# Remove samples that do not have key variables
merged <- merged %>%
  filter(!is.na(ndvi10),
         !is.na(lnd_sur10),
         !is.na(A46),
         !is.na(I05),
         (A46 + I05) > 0)

# Only looking at the LCCMR transplants 
merged <- merged %>%
  filter(Project %in% c("Urban_Lindsey", "Urban_LCCMR"))

# Create column of summed A46 and I05 reads
merged$total_AI <- merged$A46 + merged$I05

# 2. Rarefaction

# Create column of summed A46 and I05 reads
merged$total_AI <- merged$A46 + merged$I05

# Quantiles of A46/I05 read numbers
quantile(merged$total_AI, probs = c(0, 0.05, 0.10, 0.25, 0.50, 0.75, 1))

# Rarefy at 10% 
rarefy_depth <- as.integer(
  quantile(merged$total_AI, probs = 0.10)
)

before <- nrow(merged)
merged <- merged %>%
  filter(total_AI >= rarefy_depth)

# Samples lost at this depth:
before - nrow(merged)

set.seed(123)

# Function to rarefy A46 and I05 reads
rarefy_counts <- function(A46, I05, depth) {
  reads <- c(rep(1L, A46), rep(0L, I05))
  subsampled <- sample(reads, depth)
  A46_rare <- sum(subsampled)
  c(A46_rare = A46_rare,
    I05_rare = depth - A46_rare)
}

# Apply to dataframe
rare_results <- mapply(
  rarefy_counts,
  A46 = merged$A46,
  I05 = merged$I05,
  depth = rarefy_depth
)

# Add rarefied reads to dataframe
merged$A46_rare <- rare_results["A46_rare", ]
merged$I05_rare <- rare_results["I05_rare", ]

# Calculate proportions
merged$prop_A46_rare <- merged$A46_rare / rarefy_depth
merged$prop_A46 <- merged$A46 / merged$total_AI

merged$prop_I05_rare <- merged$I05_rare / rarefy_depth
merged$prop_I05 <- merged$I05 / merged$total_AI


# 3. Rarefied vs original proportions

# Check rarefaction didn't change proportions 
p_rarecheck <- ggplot(merged,
                      aes(x = prop_A46, y = prop_A46_rare)) +
  geom_point(alpha = 0.6) +
  geom_abline(slope = 1, intercept = 0,
              linetype = "dashed", color = "grey50") +
  labs(x = "Original proportion A46",
       y = "Rarefied proportion A46",
       title = "Rarefaction comparison") +
  theme_minimal()

p_rarecheck

# How are proportions dispersed?
p_hist <- ggplot(merged, aes(x = prop_A46_rare)) +
  geom_histogram(binwidth = 0.05, fill = "#588157", color = "white") +
  labs(x = "Rarefied A46 proportion",
       y = "Samples",
       title = "Distribution of rarefied proportions") +
  theme_minimal()

print(p_hist)


# 4. Beta-binomial model

# Fit null model
m_bb_null <- glmmTMB(
  cbind(A46_rare, I05_rare) ~ 1,
  family = betabinomial(link = "logit"),
  data = merged
)

# Fit alternative model
m_bb <- glmmTMB(
  cbind(A46_rare, I05_rare) ~ ndvi10,
  family = betabinomial(link = "logit"),
  data = merged
)

summary(m_bb)

# Compare models
AIC(m_bb_null, m_bb)

# Check model fit
sim_bb <- simulateResiduals(m_bb, n = 1000)
plot(sim_bb)


# 5. Logistic regression

# Add "dominant" column for logistic regression (must be binary)
merged <- merged %>% 
  mutate(dominant = factor(case_when(
    prop_A46_rare > 0.7 ~ "A46",
    prop_I05_rare > 0.7 ~ "I05",
  ), levels = c("I05", "A46")))

# Inspect, NA = neither A46 nor I05 dominant
table(merged$dominant, useNA = "always")

# Fit null model
m_logistic_null <- glm(
  dominant ~ 1,
  family = binomial(link = "logit"),
  data = merged
)

# Fit alternative model 
m_logistic <- glm(
  dominant ~ ndvi10,
  family = binomial(link = "logit"),
  data = merged
)
summary(m_logistic)

# Compare models
AIC(m_logistic_null, m_logistic)

# Check model fit
sim_log <- simulateResiduals(m_logistic, n = 1000)
plot(sim_log)


# 6. Model interpretation 

# For every 0.1 unit increase in NDVI, the odds of A46 dominance are increased by a factor of 1.822
round(exp(coef(m_logistic)["ndvi10"] * 0.1), 3)

# For every 0.1 unit increase in NDVI, the odds that a given read is A46 increase by a factor of 1.347
round(exp(fixef(m_bb)$cond["ndvi10"] * 0.1), 3)

# At NDVI 0.4, predicted A46 proportion is 0.30
pred_low  <- plogis(predict(m_bb, newdata = data.frame(ndvi10 = 0.4), type = "link"))

# At NDVI 0.8, predicted A46 proportion is 0.59
pred_high <- plogis(predict(m_bb, newdata = data.frame(ndvi10 = 0.8), type = "link"))


# 7. Model visualization

# Create vector of NDVI values
nd <- data.frame(ndvi10 = seq(min(merged$ndvi10),
                              max(merged$ndvi10),
                              length.out = 200))

# Generate predicted data based on beta binomial model for each NDVI value 
pred_bb <- predict(m_bb, newdata = nd, type = "link", se.fit = TRUE)
# Predicted mean at each NDVI
nd$bb_mean  <- plogis(pred_bb$fit)
# Confidence intervals
nd$bb_lower <- plogis(pred_bb$fit - 1.96 * pred_bb$se.fit)
nd$bb_upper <- plogis(pred_bb$fit + 1.96 * pred_bb$se.fit)

# Generate predicted data based on logisteic regression model for each NDVI value 
pred_log <- predict(m_logistic, newdata = nd, type = "link", se.fit = TRUE)
# Predicted probability at each NDVI
nd$log_prob  <- plogis(pred_log$fit)
# Confidence intervals
nd$log_lower <- plogis(pred_log$fit - 1.96 * pred_log$se.fit)
nd$log_upper <- plogis(pred_log$fit + 1.96 * pred_log$se.fit)

#Create numeric column for plotting logistic regression
merged$dominant_numeric <- ifelse(merged$dominant == "A46", 1, 0)


# Figure 3
p_bb <- ggplot() +
  geom_point(data = merged,
             aes(x = ndvi10, y = prop_A46_rare),
             alpha = 0.5, size = 2, color = "#01a7ad") +
  geom_ribbon(data = nd,
              aes(x = ndvi10, ymin = bb_lower, ymax = bb_upper),
              fill = "#01a7ad", alpha = 0.2) +
  geom_line(data = nd, aes(x = ndvi10, y = bb_mean),
            color = "#01a7ad", linewidth = 1.2) +
  labs(x = "NDVI", y = "A46 proportion",
       title = "Beta-binomial (continuous proportion)") +
  theme_classic()

p_log <- ggplot() +
  geom_jitter(data = merged,
              aes(x = ndvi10, y = dominant_numeric, color = dominant),
              width = 0, height = 0.03, size = 2.5, alpha = 0.6) +
  geom_ribbon(data = nd,
              aes(x = ndvi10, ymin = log_lower, ymax = log_upper),
              fill = "grey50", alpha = 0.2) +
  geom_line(data = nd, aes(x = ndvi10, y = log_prob),
            linewidth = 1.2, color = "grey30") +
  scale_color_manual(values = c("A46" = "#01a7ad", "I05" = "#c9ca6d"), 
                     na.translate = FALSE) +
  scale_y_continuous(breaks = c(0, 0.5, 1),
                     labels = c("I05", "50%", "A46")) +
  labs(x = "NDVI", y = "P(A46 dominance)",
       title = "Logistic Regression (which clade dominates)",
       color = NULL) +
  theme_classic()

fig <- (p_bb / p_log) +
  plot_annotation(tag_levels = "A")

print(fig)

ggsave("Beta:Log.png", fig, width = 15, height = 11, dpi = 600)
ggsave("Beta:Log.pdf", fig, width = 15, height = 11)


# 8. Who are the low NDVI values? Check they do not all come from one site

ndvi_cutoff <- 0.45

low_ndvi <- merged %>%
  filter(ndvi10 < ndvi_cutoff)

# Low NDVI samples are dispersed among 4 different sites. 
low_ndvi %>%
    select(Project, Location, ndvi10, lnd_sur10,
           A46_rare, I05_rare, prop_A46_rare, dominant) %>%
    arrange(ndvi10)
