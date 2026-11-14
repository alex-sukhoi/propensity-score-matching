# =============================================================================
# 02_matching.R
#   1. Data cleaning: harmonise yes/no coding (1/2 with stray 0s)
#   2. Propensity score (logistic regression on age, performance status,
#      neoadjuvant chemotherapy, uni/bilobar disease)
#   3. 1:1 nearest-neighbour matching without replacement, caliper 0.2 SD of
#      the logit PS (Austin, 2011); sensitivity with caliper 0.3 and with
#      exact matching on disease extent
#   4. Balance diagnostics: standardised mean differences (SMD < 0.1),
#      variance ratios, love plot, PS overlap
#   5. Outcomes in the matched sample: paired Wilcoxon (length of stay),
#      McNemar (90-day complications)
# =============================================================================

suppressPackageStartupMessages({ library(MatchIt); library(ggplot2) })
dir.create("results", showWarnings = FALSE); dir.create("figures", showWarnings = FALSE)
raw <- read.csv("data/liver_resection_registry_raw.csv")

# --- 1. Cleaning -------------------------------------------------------------
cat("neoadjuvant_chemo codes in raw export:\n"); print(table(raw$neoadjuvant_chemo))
# Codebook: 1 = Yes, 2 = No. Stray 0s are treated as "No" (natural coding);
# this assumption is logged and should be confirmed with the data owner.
d <- raw
d$neoadjuvant_chemo <- as.integer(raw$neoadjuvant_chemo == 1)
d$complication_90d  <- as.integer(raw$complication_90d == 1)
d$laparoscopic <- as.integer(d$technique == "Laparoscopic")
d$bilobar <- as.integer(d$disease_extent == "Bilobar")
write.csv(data.frame(field = "neoadjuvant_chemo", raw_code = c(0, 1, 2),
                     cleaned = c("No (assumed)", "Yes", "No"),
                     n = as.vector(table(factor(raw$neoadjuvant_chemo, levels = 0:2)))),
          "results/01_recoding_log.csv", row.names = FALSE)
write.csv(d, "data/liver_resection_registry_clean.csv", row.names = FALSE)

covs <- c("age", "performance_status", "neoadjuvant_chemo", "bilobar")
f <- reformulate(covs, response = "laparoscopic")

# --- 2-3. Propensity score and matching ---------------------------------------
m02 <- matchit(f, data = d, method = "nearest", distance = "glm", link = "linear.logit",
               ratio = 1, replace = FALSE, caliper = 0.2, std.caliper = TRUE)
m03 <- update(m02, caliper = 0.3)
m02x <- update(m02, exact = ~ bilobar)

smd_table <- function(m, label) {
  s <- summary(m, standardize = TRUE)
  rbind(data.frame(sample = "Before matching", variable = rownames(s$sum.all),
                   smd = s$sum.all[, "Std. Mean Diff."], var_ratio = s$sum.all[, "Var. Ratio"], spec = label),
        data.frame(sample = "After matching", variable = rownames(s$sum.matched),
                   smd = s$sum.matched[, "Std. Mean Diff."], var_ratio = s$sum.matched[, "Var. Ratio"], spec = label))
}
bal <- rbind(smd_table(m02, "NN 1:1, caliper 0.2"), smd_table(m03, "NN 1:1, caliper 0.3"),
             smd_table(m02x, "NN 1:1, caliper 0.2 + exact on bilobar"))
bal$smd <- round(bal$smd, 3); bal$var_ratio <- round(bal$var_ratio, 2)
write.csv(bal, "results/02_balance_smd.csv", row.names = FALSE)

pairs_tab <- function(m, label) {
  nn <- summary(m)$nn
  data.frame(specification = label, treated_total = nn["All", "Treated"],
             matched_pairs = nn["Matched", "Treated"],
             pct_treated_matched = round(100 * nn["Matched", "Treated"] / nn["All", "Treated"], 1),
             max_abs_smd_after = round(max(abs(summary(m)$sum.matched[covs, "Std. Mean Diff."])), 3))
}
spec <- rbind(pairs_tab(m02, "NN 1:1, caliper 0.2"), pairs_tab(m03, "NN 1:1, caliper 0.3"),
              pairs_tab(m02x, "NN 1:1, caliper 0.2 + exact on bilobar"))
write.csv(spec, "results/03_matching_specifications.csv", row.names = FALSE)
cat("\nMatching specifications:\n"); print(spec, row.names = FALSE)

# Primary specification = the best-balanced one among those that keep >= 80% of treated
ok <- spec[spec$pct_treated_matched >= 80, ]
primary_label <- ok$specification[which.min(ok$max_abs_smd_after)]
primary <- list("NN 1:1, caliper 0.2" = m02, "NN 1:1, caliper 0.3" = m03,
                "NN 1:1, caliper 0.2 + exact on bilobar" = m02x)[[primary_label]]
cat("Primary specification:", primary_label, "\n")

# --- 4. Balance figures --------------------------------------------------------
lp <- subset(bal, spec == primary_label & variable != "distance")
lp$variable <- factor(lp$variable, levels = rev(covs),
                      labels = rev(c("Age", "Performance status", "Neoadjuvant chemotherapy", "Bilobar disease")))
g1 <- ggplot(lp, aes(abs(smd), variable, colour = sample, shape = sample)) +
  geom_vline(xintercept = 0.1, linetype = "dashed", colour = "grey50") +
  geom_point(size = 3.5) +
  scale_colour_manual(values = c("Before matching" = "#c53030", "After matching" = "#2b6cb0")) +
  labs(title = "Covariate balance (love plot)", subtitle = paste0(primary_label, "; dashed line = |SMD| 0.1"),
       x = "Absolute standardised mean difference", y = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("figures/love_plot.png", g1, width = 7, height = 4, dpi = 110)

d$ps <- primary$distance
d$matched <- primary$weights > 0
psd <- rbind(data.frame(ps = plogis(d$ps), technique = d$technique, sample = "All patients"),
             data.frame(ps = plogis(d$ps[d$matched]), technique = d$technique[d$matched], sample = "Matched sample"))
g2 <- ggplot(psd, aes(ps, fill = technique)) +
  geom_density(alpha = 0.5, colour = NA) + facet_wrap(~ sample) +
  scale_fill_manual(values = c(Laparoscopic = "#2b6cb0", Open = "#dd6b20")) +
  labs(title = "Propensity score overlap", x = "Propensity score (probability of laparoscopy)",
       y = "Density", fill = NULL) +
  theme_minimal(base_size = 12) + theme(legend.position = "top")
ggsave("figures/ps_overlap.png", g2, width = 8, height = 3.8, dpi = 110)

# --- 5. Outcomes ---------------------------------------------------------------
md <- match.data(primary)
md <- md[order(md$subclass, -md$laparoscopic), ]
lap  <- md[md$laparoscopic == 1, ]; open <- md[md$laparoscopic == 0, ]
stopifnot(all(lap$subclass == open$subclass))
w <- wilcox.test(lap$los_days, open$los_days, paired = TRUE, exact = FALSE)
mc <- mcnemar.test(table(lap$complication_90d, open$complication_90d))
crude_los <- wilcox.test(los_days ~ technique, data = d, exact = FALSE)
outcomes <- data.frame(
  outcome = c("Length of stay, median days", "90-day complications, %"),
  laparoscopic_unmatched = c(median(d$los_days[d$laparoscopic == 1]), round(100 * mean(d$complication_90d[d$laparoscopic == 1]), 1)),
  open_unmatched = c(median(d$los_days[d$laparoscopic == 0]), round(100 * mean(d$complication_90d[d$laparoscopic == 0]), 1)),
  laparoscopic_matched = c(median(lap$los_days), round(100 * mean(lap$complication_90d), 1)),
  open_matched = c(median(open$los_days), round(100 * mean(open$complication_90d), 1)),
  test_matched = c("Paired Wilcoxon", "McNemar"),
  p_matched = signif(c(w$p.value, mc$p.value), 3))
write.csv(outcomes, "results/04_outcomes.csv", row.names = FALSE)
cat("\nOutcomes:\n"); print(outcomes, row.names = FALSE)
cat("\nBalance (primary):\n"); print(subset(bal, spec == primary_label)[, 1:4], row.names = FALSE)
