# =============================================================================
# NHANES 2003-2006 Data Loading Script
# Uses rnhanesdata built-in data + exclude_accel() for correct API usage
# =============================================================================

# --- 1. Fix library path and load packages -----------------------------------

# Must be first line — ensures updated rlang is found before locked version
.libPaths(c("C:/Users/zilef/R_libs", .libPaths()))

cat("rlang version:", as.character(packageVersion("rlang")), "\n")

cran_packages <- c("nhanesA", "dplyr", "readr", "haven")
for (pkg in cran_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, type = "binary")
  }
}

library(rnhanesdata)
library(nhanesA)
library(dplyr)
library(readr)

cat("All packages loaded.\n")

# --- 2. Load and process accelerometry ---------------------------------------

cat("Loading accelerometry data...\n")

data("PAXINTEN_C")   # 2003-2004 minute-level activity (1440 cols per day)
data("PAXINTEN_D")   # 2005-2006 minute-level activity
data("Flags_C")      # 2003-2004 wear/non-wear flags
data("Flags_D")      # 2005-2006 wear/non-wear flags

cat("Applying wear-time exclusions...\n")

# exclude_accel() returns row indices of GOOD days
# threshold_lower = 600 means >= 600 mins (10 hrs) valid wear required
keep_C <- exclude_accel(act = PAXINTEN_C, flags = Flags_C, threshold_lower = 600)
keep_D <- exclude_accel(act = PAXINTEN_D, flags = Flags_D, threshold_lower = 600)

accel_good_C <- PAXINTEN_C[keep_C, ]
accel_good_D <- PAXINTEN_D[keep_D, ]

cat(paste("Good days - Wave C:", nrow(accel_good_C),
          "| Wave D:", nrow(accel_good_D), "\n"))

# Compute per-participant summary features from minute-level data
# MIN1:MIN1440 = activity counts for each minute of the day
compute_accel_summary <- function(accel_good, wave_label) {
  
  min_cols <- paste0("MIN", 1:1440)
  
  accel_good %>%
    mutate(
      # Total activity counts for this day
      TAC_day  = rowSums(select(., all_of(min_cols)), na.rm = TRUE),
      # Sedentary minutes (counts < 100)
      ST_day   = rowSums(select(., all_of(min_cols)) < 100, na.rm = TRUE),
      # MVPA minutes (counts >= 2020)
      MVPA_day = rowSums(select(., all_of(min_cols)) >= 2020, na.rm = TRUE),
      # Light PA minutes (100 <= counts < 2020)
      LIPA_day = rowSums(select(., all_of(min_cols)) >= 100 &
                           select(., all_of(min_cols)) < 2020, na.rm = TRUE),
      # Wear time = sum of non-NA valid minutes
      WC_day   = rowSums(!is.na(select(., all_of(min_cols))))
    ) %>%
    select(SEQN, PAXCAL, PAXSTAT, SDDSRVYR,
           TAC_day, ST_day, MVPA_day, LIPA_day, WC_day) %>%
    group_by(SEQN) %>%
    summarise(
      TAC      = mean(TAC_day,  na.rm = TRUE),
      ST       = mean(ST_day,   na.rm = TRUE),
      MVPA     = mean(MVPA_day, na.rm = TRUE),
      LIPA     = mean(LIPA_day, na.rm = TRUE),
      WC       = mean(WC_day,   na.rm = TRUE),
      ndays    = n(),
      PAXCAL   = first(PAXCAL),
      PAXSTAT  = first(PAXSTAT),
      SDDSRVYR = first(SDDSRVYR),
      .groups  = "drop"
    ) %>%
    filter(ndays >= 4) %>%   # require at least 4 valid days
    mutate(
      TLAC  = log(TAC + 1),  # log-transform TAC
      wave  = wave_label
    )
}

cat("Computing activity summaries (this takes a few minutes)...\n")

accel_C <- compute_accel_summary(accel_good_C, "C")
accel_D <- compute_accel_summary(accel_good_D, "D")

accel_all <- bind_rows(accel_C, accel_D)

cat(paste("Accelerometry summary:", nrow(accel_all), "participants.\n"))

# --- 3. Helper function for safe downloads -----------------------------------

safe_nhanes <- function(name, vars) {
  tryCatch({
    df <- nhanes(name)
    df[, intersect(vars, names(df)), drop = FALSE]
  }, error = function(e) {
    cat(paste("Warning: could not load", name, "\n"))
    NULL
  })
}

# --- 4. Load demographic data ------------------------------------------------

cat("Loading demographic data...\n")

demo_vars <- c("SEQN", "RIDAGEYR", "RIAGENDR", "RIDRETH1",
               "INDFMPIR", "WTMEC2YR", "SDMVPSU", "SDMVSTRA", "SDDSRVYR")

demo_c <- safe_nhanes("DEMO_C", demo_vars)
demo_d <- safe_nhanes("DEMO_D", demo_vars)

demo_all <- bind_rows(demo_c, demo_d) %>%
  mutate(WTMEC4YR = WTMEC2YR / 2)

cat(paste("Demographics:", nrow(demo_all), "participants.\n"))

# --- 5. Load lab data --------------------------------------------------------

cat("Loading lab data...\n")
# HbA1c - 2003-2004 uses L10_C, variable is still LBXGH
ghb_c  <- safe_nhanes("L10_C",   c("SEQN", "LBXGH"))

# Triglycerides - 2003-2004 uses L13AM_C, variable is still LBXTR
# Fasting time is in a separate file PH_C for 2003-2004
trig_c <- safe_nhanes("L13AM_C", c("SEQN", "LBXTR"))
fast_c <- safe_nhanes("PH_C",    c("SEQN", "PHAFSTHR"))
trig_c <- left_join(trig_c, fast_c, by = "SEQN")

# CRP - 2003-2004 uses L11_C, variable is still LBXCRP
crp_c  <- safe_nhanes("L11_C",   c("SEQN", "LBXCRP"))

# Glucose - 2003-2004 uses L10AM_C, variable is still LBXGLU
glu_c  <- safe_nhanes("L10AM_C", c("SEQN", "LBXGLU"))

# HDL - 2003-2004 uses L13_C, variable is LBXHDD
hdl_c  <- safe_nhanes("L13_C",   c("SEQN", "LBXHDD"))

cat("Labs loaded.\n")

# --- 6. Load dietary data ----------------------------------------------------

cat("Loading dietary data...\n")

diet_vars <- c("SEQN", "DR1TKCAL", "DR1TCARB", "DR1TTFAT",
               "DR1TPROT", "DR1TFIBE", "DR1TSFAT")

diet_c <- safe_nhanes("DR1TOT_C", diet_vars)
diet_d <- safe_nhanes("DR1TOT_D", diet_vars)

diet_all <- bind_rows(diet_c, diet_d)

cat("Diet loaded.\n")

# --- 7. Merge all data -------------------------------------------------------

cat("Merging datasets...\n")

labs_c <- demo_c %>%
  left_join(ghb_c,  by = "SEQN") %>%
  left_join(trig_c, by = "SEQN") %>%
  left_join(crp_c,  by = "SEQN") %>%
  left_join(glu_c,  by = "SEQN") %>%
  left_join(hdl_c,  by = "SEQN") %>%
  left_join(bmx_c,  by = "SEQN") %>%
  left_join(bpx_c,  by = "SEQN") %>%
  left_join(smq_c,  by = "SEQN")

labs_d <- demo_d %>%
  left_join(ghb_d,  by = "SEQN") %>%
  left_join(trig_d, by = "SEQN") %>%
  left_join(crp_d,  by = "SEQN") %>%
  left_join(glu_d,  by = "SEQN") %>%
  left_join(hdl_d,  by = "SEQN") %>%
  left_join(bmx_d,  by = "SEQN") %>%
  left_join(bpx_d,  by = "SEQN") %>%
  left_join(smq_d,  by = "SEQN")

labs_all <- bind_rows(labs_c, labs_d) %>%
  mutate(WTMEC4YR = WTMEC2YR / 2)

final_df <- accel_all %>%
  left_join(labs_all, by = "SEQN") %>%
  left_join(diet_all, by = "SEQN")

cat(paste("Merged:", nrow(final_df), "participants,",
          ncol(final_df), "variables.\n"))

# --- 8. Inclusion criteria ---------------------------------------------------

cat("Applying inclusion criteria...\n")

analytic_df <- final_df %>%
  filter(
    RIDAGEYR >= 20,
    RIDAGEYR <= 85,
    !is.na(LBXGH),
    !is.na(LBXTR),
    !is.na(LBXCRP),
    !is.na(PHAFSTHR),
    PHAFSTHR >= 8
  )

cat(paste("Analytic sample:", nrow(analytic_df), "participants.\n"))

# --- 9. Summary --------------------------------------------------------------

cat("\n--- Outcomes ---\n")
print(summary(analytic_df[, c("LBXGH", "LBXTR", "LBXCRP")]))

cat("\n--- Activity ---\n")
print(summary(analytic_df[, c("TAC", "TLAC", "ST", "MVPA", "LIPA")]))

cat("\n--- Demographics ---\n")
cat("Age (mean):  ", round(mean(analytic_df$RIDAGEYR, na.rm=TRUE), 1), "\n")
cat("Female (%):  ", round(mean(analytic_df$RIAGENDR == 2, na.rm=TRUE)*100, 1), "\n")
cat("Wave C (n):  ", sum(analytic_df$wave == "C", na.rm=TRUE), "\n")
cat("Wave D (n):  ", sum(analytic_df$wave == "D", na.rm=TRUE), "\n")

# --- 10. Export --------------------------------------------------------------

if (!dir.exists("data")) dir.create("data")
write_csv(analytic_df, "data/nhanes_analytic_23062026.csv")

cat("\nSaved: data/nhanes_analytic.csv\n")
cat("Columns:", paste(names(analytic_df), collapse = ", "), "\n")
cat("Done.\n")
