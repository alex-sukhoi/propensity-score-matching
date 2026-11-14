# =============================================================================
# 01_simulate_data.R
# Synthetic surgical registry: laparoscopic vs open liver resection (n = 160).
# Treatment assignment depends on the covariates (confounding by indication):
# laparoscopy is chosen more often for younger, fitter patients with
# unilobar disease. The raw export also reproduces a common data problem:
# yes/no fields coded 1 = Yes, 2 = No, but with stray 0s.
# =============================================================================

set.seed(77)
RNGkind(sample.kind = "Rejection")
n <- 160

age <- round(rnorm(n, 63, 11))
performance_status <- sample(0:2, n, TRUE, c(0.50, 0.38, 0.12))       # ECOG 0-2
neoadjuvant_chemo <- rbinom(n, 1, 0.55)
bilobar <- rbinom(n, 1, 0.40)
bmi <- round(rnorm(n, 27, 4), 1)

lin <- -0.2 - 0.035 * (age - 63) - 0.5 * performance_status - 1.3 * bilobar +
       0.2 * neoadjuvant_chemo - 0.03 * (bmi - 27)
laparoscopic <- rbinom(n, 1, plogis(lin))

# Outcomes: laparoscopy shortens stay (true effect -2 days) and lowers complications
los_days <- round(pmax(2, rlnorm(n, log(7.5) + 0.03 * performance_status + 0.15 * bilobar +
                                   0.01 * (age - 63) - 0.25 * laparoscopic, 0.35)))
complication_90d <- rbinom(n, 1, plogis(-1.2 + 0.4 * performance_status + 0.5 * bilobar +
                                         0.02 * (age - 63) - 0.6 * laparoscopic))

# Export as the registry would deliver it: 1 = Yes, 2 = No, plus stray zeros
yn <- function(x) ifelse(x == 1, 1, 2)
raw <- data.frame(
  patient_id = sprintf("H%03d", 1:n),
  technique = ifelse(laparoscopic == 1, "Laparoscopic", "Open"),
  age, performance_status, bmi,
  neoadjuvant_chemo = yn(neoadjuvant_chemo),
  disease_extent = ifelse(bilobar == 1, "Bilobar", "Unilobar"),
  los_days, complication_90d = yn(complication_90d)
)
zero_rows <- sample(which(neoadjuvant_chemo == 0), 25)      # "No" exported as 0
raw$neoadjuvant_chemo[zero_rows] <- 0

dir.create("data", showWarnings = FALSE)
write.csv(raw, "data/liver_resection_registry_raw.csv", row.names = FALSE)
cat(sprintf("Saved raw registry: %d patients (%d laparoscopic, %d open)\n",
            n, sum(laparoscopic), sum(1 - laparoscopic)))
