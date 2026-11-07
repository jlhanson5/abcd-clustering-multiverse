# =============================================================================
# 01_data.R: build the analysis files.
#   data_derived/long.rds : one row per child x wave (BL-Y3), observed values
#   data_derived/wide.rds : one row per child, observed + single-imputed CBCL,
#                           income, covariates, family/pilot/inclusion flags
# Sample flow is written to output/logs/sample_flow.txt.
# =============================================================================

source("00_setup.R")
suppressPackageStartupMessages(library(mice))
dir.create(DER, showWarnings = FALSE); dir.create(LOGS, showWarnings = FALSE)
unlink(file.path(LOGS, "sample_flow.txt"))

# ---- CBCL, age, site, family --------------------------------------------------
cbcl <- read_abcd("mental-health/mh_p_cbcl.csv",
                  c("cbcl_scr_syn_totprob_r", "cbcl_scr_syn_totprob_t",
                    "cbcl_scr_syn_internal_r", "cbcl_scr_syn_external_r")) |>
  filter(eventname %in% EVENTS) |>
  mutate(across(starts_with("cbcl_"), num))

lt <- read_abcd("abcd-general/abcd_y_lt.csv",
                c("interview_age", "site_id_l", "rel_family_id")) |>
  filter(eventname %in% EVENTS) |>
  mutate(age = as.numeric(interview_age) / 12)

long <- lt |>
  select(src_subject_id, eventname, age) |>
  left_join(cbcl, by = c("src_subject_id", "eventname")) |>
  mutate(wave = factor(names(EVENTS)[match(eventname, EVENTS)], levels = WAVES),
         wave_num = as.integer(wave) - 1L,
         tot_r = cbcl_scr_syn_totprob_r, tot_t = cbcl_scr_syn_totprob_t) |>
  select(src_subject_id, wave, wave_num, age, tot_r, tot_t,
         int_r = cbcl_scr_syn_internal_r, ext_r = cbcl_scr_syn_external_r)

bl <- lt |> filter(eventname == EV_BL) |>
  transmute(src_subject_id, site = site_id_l, family = rel_family_id, age_bl = age)

# ---- Baseline income and sex ---------------------------------------------------
demo <- read_abcd("abcd-general/abcd_p_demo.csv", c("demo_comb_income_v2", "demo_sex_v2")) |>
  filter(eventname == EV_BL) |>
  transmute(src_subject_id,
            income_bin = num(demo_comb_income_v2),
            income = unname(INCOME_MID[as.character(income_bin)]),
            sex = case_when(demo_sex_v2 == "1" ~ "M", demo_sex_v2 == "2" ~ "F"))

# ---- Sample flow ----------------------------------------------------------------
nwave <- long |> group_by(src_subject_id) |> summarise(n_waves = sum(!is.na(tot_r)), .groups = "drop")
base <- bl |> left_join(demo, by = "src_subject_id") |> left_join(nwave, by = "src_subject_id")

log_line("Sample flow (ABCD 5.1, CBCL total problems BL-Y3), run ", format(Sys.time()))
log_line("Baseline participants: ", nrow(base))
s1 <- filter(base, !is.na(income))
log_line("Excluded, income refuse/don't know (777/999) or missing: ", nrow(base) - nrow(s1),
         "  -> ", nrow(s1))
s2 <- filter(s1, !is.na(sex))
log_line("Excluded, sex not 1/2 (intersex-male or missing): ", nrow(s1) - nrow(s2), "  -> ", nrow(s2))
s3 <- filter(s2, n_waves >= 3)
log_line("Excluded, fewer than 3 of 4 CBCL waves: ", nrow(s2) - nrow(s3), "  -> ", nrow(s3))
log_line("  of which all 4 waves: ", sum(s3$n_waves == 4))

# ---- Wide file + single seeded imputation (raw and T jointly) -------------------
long <- semi_join(long, s3, by = "src_subject_id")
wide <- long |>
  select(src_subject_id, wave, age, tot_r, tot_t) |>
  pivot_wider(names_from = wave, values_from = c(age, tot_r, tot_t)) |>
  right_join(s3, by = "src_subject_id") |>
  mutate(log_income = log(income))

# Missing ages at missing waves: baseline age + median wave offset (only used for
# imputed values under the age time metric).
offs <- long |> left_join(select(bl, src_subject_id, age_bl), by = "src_subject_id") |>
  group_by(wave) |> summarise(off = median(age - age_bl, na.rm = TRUE))
for (w in WAVES) {
  a <- paste0("age_", w)
  wide[[a]] <- coalesce(wide[[a]], wide$age_bl + offs$off[offs$wave == w])
}

imp_vars <- c(paste0("tot_r_", WAVES), paste0("tot_t_", WAVES), "age_bl", "sex", "log_income")
imp_dat <- wide |> select(all_of(imp_vars)) |> mutate(sex = factor(sex))
imp <- mice(imp_dat, m = 1, method = "pmm", seed = SEED, printFlag = FALSE)
done <- complete(imp, 1)
for (w in WAVES) {
  wide[[paste0("imp_r_", w)]] <- done[[paste0("tot_r_", w)]]
  wide[[paste0("imp_t_", w)]] <- done[[paste0("tot_t_", w)]]
}
log_line("Imputed CBCL values (raw total, single pmm imputation, seed ", SEED, "): ",
         sum(is.na(as.matrix(wide[paste0("tot_r_", WAVES)]))))

# ---- Family and pilot flags -----------------------------------------------------
set.seed(SEED)
wide <- wide |>
  group_by(family) |>
  mutate(one_per_family = row_number() == sample.int(n(), 1)) |>
  ungroup()
log_line("One child per family (seeded random pick): ", sum(wide$one_per_family),
         " of ", nrow(wide), " (", n_distinct(wide$family), " families)")

set.seed(SEED + 1)
fams <- sample(unique(wide$family))
csum <- cumsum(table(wide$family)[fams])
pilot_fams <- fams[seq_len(which(csum >= PILOT_N)[1])]
wide$pilot <- wide$family %in% pilot_fams
log_line("Pilot (whole families sampled until >= ", PILOT_N, " children): ", sum(wide$pilot),
         " children, ", length(pilot_fams), " families")

wide$all4 <- wide$n_waves == 4

# Children with no interview_age at all (1 in 5.1): fill baseline and wave ages
# with the sample median at each wave. Done after imputation and flags so it
# changes nothing else. Age is needed for the age time metric and as a
# poverty-model covariate.
no_age <- is.na(wide$age_bl)
if (any(no_age)) {
  wide$age_bl[no_age] <- median(wide$age_bl, na.rm = TRUE)
  for (w in WAVES) {
    a <- paste0("age_", w)
    wide[[a]][no_age] <- median(wide[[a]], na.rm = TRUE)
  }
  log_line("Ages filled with wave medians (no interview_age at any wave): ", sum(no_age), " child(ren)")
  fill <- select(wide, src_subject_id, all_of(paste0("age_", WAVES)))
  long <- long |>
    left_join(fill |> pivot_longer(-src_subject_id, names_to = "wave", names_prefix = "age_",
                                   values_to = "age_fill"),
              by = c("src_subject_id", "wave")) |>
    mutate(age = coalesce(age, age_fill), wave = factor(wave, levels = WAVES)) |>
    select(-age_fill)
}
long <- long |> left_join(select(wide, src_subject_id, one_per_family, pilot, all4), by = "src_subject_id")

saveRDS(long, file.path(DER, "long.rds"))
saveRDS(wide, file.path(DER, "wide.rds"))
log_line("Wrote ", file.path(DER, "long.rds"), " and ", file.path(DER, "wide.rds"))
