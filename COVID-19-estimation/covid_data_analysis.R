# covid_data_analysis.R
# Karnataka COVID-19 data analysis -- produces all Section 4 numbers.
# Files: ann2.csv, contacts.csv, 54_serial_interval_data.xlsx,
#        Delays_3_for_histogram.xlsx

cat('Current working directory', getwd())
setwd(getwd())

library(readxl)

ann2     <- read.csv("covid-19-data/ann2.csv", stringsAsFactors = FALSE)
contacts <- read.csv("covid-19-data/contacts.csv", stringsAsFactors = FALSE)
si_data  <- read_xlsx("covid-19-data/54_serial_interval_data.xlsx")
delays   <- read_xlsx("covid-19-data/Delays_3 for histogram.xlsx")

# Month filter
get_month <- function(s) {
  m <- c(Mar=3, Apr=4, May=5)
  for (nm in names(m))
    if (grepl(nm, as.character(s), ignore.case=TRUE)) return(m[nm])
  NA_integer_
}
ann2$month <- sapply(ann2$dt_conf, get_month)
sub <- ann2[!is.na(ann2$month) & ann2$month %in% 3:5, ]
cat(sprintf("March-May total: %d\n", nrow(sub)))
print(sort(table(sub$cat_4), decreasing=TRUE))

# Eligibility
traced   <- sub[sub$cat_4 == "Local Traced", ]
untraced <- sub[sub$cat_4 != "Local Traced", ]
pobs <- nrow(untraced) / nrow(sub)
cat(sprintf("Excluded: %d | Included: %d | p_obs=%.4f\n",
            nrow(traced), nrow(untraced), pobs))

# Save untrace data for parameter estimation
write.csv(untraced, "covid-19-data/untraced.csv", row.names = FALSE)

# j_obs
cat("\nj_obs distribution:\n")
print(table(untraced$children_pri))
cat(sprintf("mean=%.4f, max=%d, j_obs=0: %.1f%%\n",
            mean(untraced$children_pri), max(untraced$children_pri),
            100*mean(untraced$children_pri==0)))

# Cluster
cat(sprintf("cluster_size: mean=%.2f, median=%d, max=%d\n",
            mean(untraced$cluster_size),
            as.integer(median(untraced$cluster_size)),
            max(untraced$cluster_size)))

# Serial interval
si_vals <- si_data$serial_interval
cat(sprintf("SI (n=%d): mean=%.2f, sd=%.2f, range=[%d,%d]\n",
            length(si_vals), mean(si_vals), sd(si_vals),
            min(si_vals), max(si_vals)))

# Delays
d_conf <- "Symptom onset to lab confirmation"
cat(sprintf("Onset->confirm (n=%d): mean=%.2f, sd=%.2f\n",
            sum(!is.na(delays[[d_conf]])),
            mean(delays[[d_conf]], na.rm=TRUE),
            sd(delays[[d_conf]], na.rm=TRUE)))

# Timing availability
delays_ids <- as.integer(delays$Case)
un_onset   <- sum(untraced$ID %in% delays_ids)
cat(sprintf("With onset: %d/%d (%.1f%%) | Fully latent: %d/%d (%.1f%%)\n",
            un_onset, nrow(untraced), 100*un_onset/nrow(untraced),
            nrow(untraced)-un_onset, nrow(untraced),
            100*(nrow(untraced)-un_onset)/nrow(untraced)))
cat(sprintf("Index cases (parentid=0): %d | Infector known: %d\n",
            sum(untraced$parentid==0), sum(untraced$parentid>0)))

# Summary
cat("\n=== FINAL COHORT ===\n")
cat(sprintf("Period: 9 Mar - 31 May 2020\n"))
cat(sprintf("Total: %d | Excluded: %d | Included: %d\n",
            nrow(sub), nrow(traced), nrow(untraced)))
cat(sprintf("p_obs: %.4f | j_obs max: %d | K_MAX: %d\n",
            pobs, max(untraced$children_pri), max(untraced$cluster_size)))
cat(sprintf("SI mean: %.2f d (n=54) | Onset->confirm: %.2f d (n=261)\n",
            mean(si_vals), mean(delays[[d_conf]], na.rm=TRUE)))

