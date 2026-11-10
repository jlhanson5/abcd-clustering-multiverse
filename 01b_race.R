# =============================================================================
# 01b_race.R: race/ethnicity for the race-adjusted poverty models (added 2026-11-09).
#   Rscript 01b_race.R
# Reads race_ethnicity from abcd-general/abcd_p_demo.csv (ABCD 5.1; ABCD's derived
# 5-level variable: 1 White, 2 Black, 3 Hispanic, 4 Asian, 5 Other). Baseline value;
# a child blank at baseline takes a later wave's code (codes do not differ across
# waves in 5.1). Still blank -> NA, left out of the adjusted model only.
# Output: data_derived/race.rds (src_subject_id, race) for every child in wide.rds,
# and output/logs/race.txt. wide.rds and long.rds are not changed.
# =============================================================================

source("00_setup.R")
logf <- file.path(LOGS, "race.txt"); unlink(logf)

wide <- readRDS(file.path(DER, "wide.rds"))
d <- read_abcd("abcd-general/abcd_p_demo.csv", "race_ethnicity") |>
  filter(!is.na(race_ethnicity), race_ethnicity != "")
stopifnot(all(d$race_ethnicity %in% names(RACE_LEV)))
n_conflict <- d |> group_by(src_subject_id) |> summarise(n = n_distinct(race_ethnicity)) |>
  filter(n > 1) |> nrow()
first <- d |> arrange(src_subject_id, eventname != EV_BL) |> distinct(src_subject_id, .keep_all = TRUE)

race <- tibble(src_subject_id = wide$src_subject_id) |>
  left_join(first, by = "src_subject_id") |>
  transmute(src_subject_id, race = unname(RACE_LEV[race_ethnicity]),
            from_baseline = eventname == EV_BL)
saveRDS(race, file.path(DER, "race.rds"))

log_line("Race/ethnicity (race_ethnicity, ABCD 5.1), run ", format(Sys.time()), file = logf)
log_line("Children in wide.rds: ", nrow(race), file = logf)
log_line("Ids with conflicting codes across waves (whole table): ", n_conflict, file = logf)
log_line("Code taken from a later wave (blank at baseline): ", sum(race$from_baseline %in% FALSE), file = logf)
log_line("No code at any wave (left out of the adjusted model): ", sum(is.na(race$race)), file = logf)
for (s in list(list("All children", TRUE), list("One per family", wide$one_per_family))) {
  tab <- table(factor(race$race[s[[2]]], levels = RACE_LEV), useNA = "always")
  log_line(s[[1]], ": ", paste(sprintf("%s %d (%.1f%%)", replace(names(tab), is.na(names(tab)), "Missing"),
                                       tab, 100 * tab / sum(tab)), collapse = ", "), file = logf)
}
