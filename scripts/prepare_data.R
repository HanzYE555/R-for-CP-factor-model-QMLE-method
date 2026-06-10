```r
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(lubridate)
  library(stringr)
})
# -------------------------
# Config
# -------------------------
zip_path_ff3   <- file.path("data", "raw", "F-F_Research_Data_Factors_CSV.zip")
zip_path_10x10 <- file.path("data", "raw", "100_Portfolios_10x10_CSV.zip")

out_dir <- file.path("data", "processed")
out_rds <- file.path(out_dir, "ten_by_ten_prepared.rds")

start_date <- ymd("1964-01-01")
end_date   <- ymd("2021-12-01")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path) {
  if (!file.exists(path)) stop("Missing file: ", path, call. = FALSE)
}
stop_if_missing(zip_path_ff3)
stop_if_missing(zip_path_10x10)

# -------------------------
# Helpers
# -------------------------
ensure_num <- function(x) {
  readr::parse_number(as.character(x),
                      locale = readr::locale(decimal_mark = ".", grouping_mark = ","))
}

is_yyyymm <- function(x) grepl("^\\s*[12][0-9]{5}\\s*$", as.character(x))

# unzip and return extracted csv path
extract_first_csv <- function(zip_path) {
  zlist <- unzip(zip_path, list = TRUE)
  csvs <- zlist$Name[grepl("\\.csv$", tolower(zlist$Name))]
  if (length(csvs) == 0) stop("No CSV found inside zip: ", zip_path, call. = FALSE)
  inner_csv <- csvs[1]
  tmp_dir <- tempdir()
  unzip(zip_path, files = inner_csv, exdir = tmp_dir, overwrite = TRUE)
  file.path(tmp_dir, inner_csv)
}

# Read a FF-style CSV that has text headers and then a table where first col is YYYYMM.
# Strategy:
# - find first line that starts with YYYYMM
# - treat the line above it as the header row
# - read from that header row onward with readr::read_csv
read_monthly_block_from_csv <- function(csv_path) {
  raw_lines <- readLines(csv_path, warn = FALSE, encoding = "UTF-8")
  data_pat  <- "^\\s*[12][0-9]{5}\\b"
  first_data_idx <- which(grepl(data_pat, raw_lines))[1]
  if (is.na(first_data_idx) || first_data_idx <= 1) {
    stop("Failed to locate monthly data block in: ", csv_path, call. = FALSE)
  }
  header_idx <- first_data_idx - 1
  block <- raw_lines[header_idx:length(raw_lines)]
  
  tf <- tempfile(fileext = ".csv")
  writeLines(block, tf)
  
  suppressWarnings(readr::read_csv(
    tf,
    col_types = readr::cols(.default = readr::col_character()),
    show_col_types = FALSE,
    progress = FALSE
  ))
}

# CAPM residuals: regress excess returns on MktRF column-by-column
capm_resid_excess <- function(R, rf, mktrf) {
  stopifnot(is.matrix(R), is.numeric(R), length(rf) == nrow(R), length(mktrf) == nrow(R))
  Rex <- sweep(R, 1, rf, FUN = "-")
  Tn <- nrow(Rex); Nn <- ncol(Rex)
  E <- matrix(NA_real_, Tn, Nn)
  
  for (j in seq_len(Nn)) {
    fit <- lm(Rex[, j] ~ mktrf)
    E[, j] <- resid(fit)
  }
  colnames(E) <- colnames(R)
  rownames(E) <- rownames(R)
  E
}

# -------------------------
# A) Read FF3 factors (monthly)
# -------------------------
csv_path_ff3 <- extract_first_csv(zip_path_ff3)
ff_raw <- read_monthly_block_from_csv(csv_path_ff3)

names(ff_raw)[1] <- "Date"
names(ff_raw) <- trimws(names(ff_raw))

ff_monthly <- ff_raw %>%
  filter(is_yyyymm(Date)) %>%
  mutate(
    Date  = ensure_num(Date),
    date  = ymd(paste0(as.character(Date), "01"))
  )

# Normalize column names
if ("Mkt-RF" %in% names(ff_monthly)) {
  names(ff_monthly)[names(ff_monthly) == "Mkt-RF"] <- "MktRF"
}
required <- c("MktRF", "SMB", "HML", "RF")
miss <- setdiff(required, names(ff_monthly))
if (length(miss) > 0) stop("FF3 is missing columns: ", paste(miss, collapse = ", "), call. = FALSE)

ff3 <- ff_monthly %>%
  transmute(
    date  = date,
    MktRF = ensure_num(MktRF) / 100,
    SMB   = ensure_num(SMB)   / 100,
    HML   = ensure_num(HML)   / 100,
    RF    = ensure_num(RF)    / 100
  ) %>%
  filter(date >= start_date, date <= end_date) %>%
  arrange(date)

fac_monthly <- ff3 %>% select(date, MktRF, RF)

cat(sprintf("FF3: %d rows, sample %s — %s\n",
            nrow(ff3), format(min(ff3$date)), format(max(ff3$date))))

# -------------------------
# -------------------------
# B) Read 10x10 portfolios (monthly) - robust mode
# -------------------------
csv_path_10x10 <- extract_first_csv(zip_path_10x10)

raw_lines <- readLines(csv_path_10x10, warn = FALSE, encoding = "UTF-8")
data_pat  <- "^\\s*[12][0-9]{5}\\b"
first_data_idx <- which(grepl(data_pat, raw_lines))[1]
if (is.na(first_data_idx)) stop("Failed to locate YYYYMM rows in 10x10 CSV.", call. = FALSE)

# take only the contiguous monthly block: YYYYMM... until first non-YYYYMM after start
tail_lines <- raw_lines[first_data_idx:length(raw_lines)]
is_data_line <- grepl(data_pat, tail_lines)
end_idx <- which(!is_data_line)[1]
if (is.na(end_idx)) {
  block_lines <- tail_lines
} else {
  block_lines <- tail_lines[seq_len(end_idx - 1)]
}

tf <- tempfile(fileext = ".csv")
writeLines(block_lines, tf)

# Read without headers to avoid header/format issues
ten_tab <- suppressWarnings(readr::read_csv(
  tf,
  col_names = FALSE,
  col_types = readr::cols(.default = readr::col_character()),
  show_col_types = FALSE,
  progress = FALSE
))

# We expect: col1 = YYYYMM, next 100 cols = returns
if (ncol(ten_tab) < 101) {
  stop(sprintf("10x10 block has only %d columns (<101). File format may differ.", ncol(ten_tab)), call. = FALSE)
}
ten_tab <- ten_tab[, 1:101]

names(ten_tab)[1] <- "Date"
port_cols <- sprintf("P%03d", 1:100)
names(ten_tab)[2:101] <- port_cols

ten_monthly <- ten_tab %>%
  filter(is_yyyymm(Date)) %>%
  mutate(
    Date_num = ensure_num(Date),
    date     = ymd(paste0(as.character(Date_num), "01"))
  ) %>%
  filter(!is.na(date)) %>%
  filter(date >= start_date, date <= end_date) %>%
  arrange(date)

# Convert to numeric decimal returns; set missing to 0
ten_monthly_num <- ten_monthly
for (cl in port_cols) {
  x <- ensure_num(ten_monthly_num[[cl]])
  
  # common FF missing codes
  x[!is.na(x) & x <= -90] <- NA_real_
  
  # percent -> decimal
  x <- x / 100
  
  # set missing to 0
  x[is.na(x)] <- 0
  
  ten_monthly_num[[cl]] <- x
}

cat(sprintf("10x10: %d rows after date filter (%s — %s)\n",
            nrow(ten_monthly_num),
            format(min(ten_monthly_num$date)),
            format(max(ten_monthly_num$date))))
# -------------------------
# C) Align and build matrices
# -------------------------
data_merged <- ten_monthly_num %>%
  select(date, all_of(port_cols)) %>%
  inner_join(fac_monthly, by = "date") %>%
  arrange(date)

R_mat <- as.matrix(data_merged[, port_cols, drop = FALSE])
rownames(R_mat) <- format(data_merged$date, "%Y-%m")

MktRF <- as.numeric(data_merged$MktRF)
RF    <- as.numeric(data_merged$RF)

cat(sprintf("Merged: %d months, %d portfolios\n", nrow(data_merged), ncol(R_mat)))


col_na <- colSums(is.na(R_mat))
summary(col_na)
max(col_na)
# -------------------------
# D) CAPM residuals
# -------------------------
# Note: lm() will drop rows with NA automatically per column.
# If you have many NAs, you may want a stricter handling policy.
E_capm <- capm_resid_excess(R_mat, RF, MktRF)

# -------------------------
# E) Save
# -------------------------
saveRDS(list(
  date = data_merged$date,
  R_mat = R_mat,
  MktRF = MktRF,
  RF = RF,
  E_capm = E_capm,
  colnames_10x10 = colnames(R_mat)
), file = out_rds)

cat("Saved: ", out_rds, "\n", sep = "")
cat(sprintf("mean(RF)=%.6f, sd(RF)=%.6f\n", mean(RF, na.rm = TRUE), sd(RF, na.rm = TRUE)))
cat(sprintf("E_capm: mean=%.6f, sd=%.6f\n",
            mean(E_capm, na.rm = TRUE), sd(as.vector(E_capm), na.rm = TRUE)))
