# Explicit acquisition from original public sources; no fund data are bundled.
# The private reader is deliberately limited to the source archive's MATLAB v5
# real double arrays R and F. It is not a general MATLAB-file interface.
.gpy_mat_numeric <- function(path) {
  bytes <- readBin(path, "raw", n = file.info(path)$size)
  if (length(bytes) < 128L || rawToChar(bytes[127:128]) != "IM")
    stop("Expected the author's little-endian MATLAB v5 file.")
  element <- function(buf, pos) {
    if (pos + 7L > length(buf)) stop("Truncated MATLAB header.")
    tag <- readBin(buf[pos:(pos + 3L)], "integer", n = 1L, size = 4L, endian = "little")
    small <- bitwShiftR(tag, 16L)
    if (small > 0L) {
      kind <- bitwAnd(tag, 65535L); size <- small; begin <- pos + 4L; next_pos <- pos + 8L
      if (size > 4L) stop("Invalid MATLAB small element.")
    } else {
      kind <- tag
      size <- readBin(buf[(pos + 4L):(pos + 7L)], "integer", n = 1L, size = 4L, endian = "little")
      if (is.na(size) || size < 0L) stop("Invalid MATLAB element size.")
      begin <- pos + 8L
      next_pos <- if (kind == 15L) begin + size else begin + ceiling(size / 8) * 8
    }
    if (begin + size - 1 > length(buf)) stop("Truncated MATLAB payload.")
    list(kind = kind, data = if (size) buf[seq.int(begin, length.out = size)] else raw(0), next_pos = next_pos)
  }
  parse <- function(buf) {
    out <- list(); pos <- 1L
    while (pos + 7L <= length(buf)) {
      x <- element(buf, pos); pos <- x$next_pos
      if (x$kind == 15L) {
        out <- c(out, parse(memDecompress(x$data, type = "gzip")))
      } else if (x$kind == 14L) {
        flags <- element(x$data, 1L)
        dims <- element(x$data, flags$next_pos)
        shape <- readBin(dims$data, "integer", n = length(dims$data) / 4L, size = 4L, endian = "little")
        name <- element(x$data, dims$next_pos); label <- rawToChar(name$data)
        values <- element(x$data, name$next_pos)
        if (label %in% c("R", "F")) {
          flag <- readBin(flags$data, "integer", n = 1L, size = 4L, endian = "little")
          if (bitwAnd(flag, 255L) != 6L || bitwAnd(flag, 2048L) != 0L ||
              values$kind != 9L || length(shape) != 2L || any(shape <= 0L) ||
              prod(shape) * 8 != length(values$data)) stop("R/F must be ordinary real double matrices.")
          out[[label]] <- matrix(readBin(values$data, "double", n = prod(shape),
            size = 8L, endian = "little"), nrow = shape[1L], ncol = shape[2L])
        }
      }
    }
    out
  }
  out <- parse(bytes[-seq_len(128L)])
  if (!all(c("R", "F") %in% names(out))) stop("MAT file must contain R and F.")
  out
}

.gpy_rf_csv <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[grepl("^[[:space:]]*[0-9]{6}[[:space:]]*,", lines)]
  if (!length(lines)) stop("No monthly YYYYMM records in the French archive.")
  values <- utils::read.csv(text = paste(lines, collapse = "\n"), header = FALSE,
    strip.white = TRUE, stringsAsFactors = FALSE)
  if (ncol(values) < 5L) stop("French factor archive lacks RF.")
  dates <- sprintf("%06d", as.integer(values[[1L]]))
  wanted <- format(seq(as.Date("1987-02-01"), by = "month", length.out = 368L), "%Y%m")
  index <- match(wanted, dates)
  if (anyNA(index) || anyDuplicated(dates)) stop("Missing or duplicated article dates.")
  rf <- as.numeric(values[[5L]][index])
  if (any(!is.finite(rf))) stop("Non-finite French RF.")
  rf
}

.gpy_numeric_md5 <- function(x) {
  path <- tempfile("fund_numeric_checksum_")
  on.exit(unlink(path), add = TRUE)
  con <- file(path, "wb")
  tryCatch(writeBin(as.double(x), con, size = 8L, endian = "little"), finally = close(con))
  unname(tools::md5sum(path))
}

.gpy_write_csv <- function(x, path) {
  if (is.null(dim(x))) x <- matrix(x, ncol = 1L)
  text <- matrix(formatC(as.double(x), digits = 17L, format = "g", decimal.mark = "."),
                 nrow = nrow(x), ncol = ncol(x))
  utils::write.table(text, path, sep = ",", quote = FALSE, row.names = FALSE, col.names = FALSE)
}

#' Acquire raw article inputs directly from their original public sources
#' @export
fetch_empirical_data <- function(out_dir, author_archive = NULL, rf_archive = NULL,
    author_url = paste0("https://www.dropbox.com/scl/fo/v9qgrwn5xgx7ahwdyqhki/",
      "AFJ9_QuPz2BrHEXfl29BQlU?rlkey=sn2za70nnwojcw8jnogw34uhs&dl=1"),
    rf_url = paste0("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/",
      "Historical_Archives/08%202018%20Update/ftp/F-F_Research_Data_Factors_CSV.zip"),
    verify = TRUE, quiet = FALSE) {
  restore_locale <- .tf_utf8_scope(); on.exit(restore_locale(), add = TRUE)
  if (!is.logical(verify) || length(verify) != 1L || is.na(verify) ||
      !is.logical(quiet) || length(quiet) != 1L || is.na(quiet)) stop("Invalid verify/quiet flags.")
  out_dir <- .emp_output_dir(out_dir)
  if (is.null(out_dir)) stop("An empty output directory is required.")
  old_timeout <- getOption("timeout")
  options(timeout = max(60, old_timeout))
  on.exit(options(timeout = old_timeout), add = TRUE)
  temporary <- tempfile("author_input_acquisition_"); dir.create(temporary)
  on.exit(unlink(temporary, recursive = TRUE), add = TRUE)
  obtain <- function(local, url, basename) {
    if (!is.null(local)) {
      if (!is.character(local) || length(local) != 1L || is.na(local) || !file.exists(local))
        stop("Supplied source archive does not exist.")
      return(local)
    }
    if (!is.character(url) || length(url) != 1L || !grepl("^https://", url)) stop("Source URLs must be HTTPS.")
    dest <- file.path(temporary, basename)
    status <- tryCatch(utils::download.file(url, dest, mode = "wb", quiet = quiet),
      error = function(e) stop("Original-source download failed. Download the two source ZIPs in a browser and supply author_archive and rf_archive for offline processing. Details: ",
                               conditionMessage(e), call. = FALSE))
    if (!identical(as.integer(status), 0L)) stop("Source download failed.")
    dest
  }
  author <- obtain(author_archive, author_url, "author_code.zip")
  french <- obtain(rf_archive, rf_url, "french_2018.zip")
  entries <- utils::unzip(author, list = TRUE)$Name
  entry <- entries[grepl("(^|/)CaseStudy/mutual_fund[.]mat$", entries)]
  unsafe <- function(x) grepl("(^/|^[A-Za-z]:|(^|/)\\.\\.(/|$)|\\\\)", x)
  if (length(entry) != 1L || unsafe(entry)) stop("ZIP must contain one safe CaseStudy/mutual_fund.mat entry.")
  mat_file <- utils::unzip(author, files = entry, exdir = temporary)
  parsed <- .gpy_mat_numeric(mat_file)
  if (!identical(dim(parsed$R), c(668L, 726L)) || !identical(dim(parsed$F), c(668L, 4L)))
    stop("Public archive dimensions differ from the article snapshot.")
  # All 726 original fund series are retained: do not run the separate author's
  # less-than-20-percent-missing case-study screen in RunRegression.m.
  fund <- t(parsed$R[301:668, , drop = FALSE]); factors <- parsed$F[301:668, , drop = FALSE]
  entries <- utils::unzip(french, list = TRUE)$Name
  rf_entry <- entries[grepl("(^|/)F-F_Research_Data_Factors[.]csv$", entries, ignore.case = TRUE)]
  if (length(rf_entry) != 1L || unsafe(rf_entry)) stop("French ZIP must contain one safe monthly factor CSV.")
  rf_file <- utils::unzip(french, files = rf_entry, exdir = temporary)
  rf <- .gpy_rf_csv(rf_file)
  checksums <- c(fund = .gpy_numeric_md5(fund), factor = .gpy_numeric_md5(factors), rf = .gpy_numeric_md5(rf))
  expected <- c(fund = "e4be121319019e5061e5c82611816cd7", factor = "1d2bc45c111f1dad9436d5b6ecbee7c3",
                rf = "207a14725af16951c69cb3815f1a9037")
  if (verify && !identical(checksums, expected))
    stop("Source numeric checksums differ from the article inputs. Inspect the archive vintage; no changed input is silently substituted.")
  files <- file.path(out_dir, c("fund.csv", "factor.csv", "rf.csv"))
  .gpy_write_csv(fund, files[1L]); .gpy_write_csv(factors, files[2L]); .gpy_write_csv(rf, files[3L])
  panel <- read_fund_panel(files[1L], files[2L], files[3L])
  panel$acquisition <- list(author_url = author_url, rf_url = rf_url,
    source_archive_md5 = c(author = unname(tools::md5sum(author)), french = unname(tools::md5sum(french))),
    numeric_md5 = checksums, numeric_verified = identical(checksums, expected),
    rf_vintage = "French July 2018 data cut (August 2018 update)",
    transformation = "R[301:668,] transposed; F[301:668,]; all726 funds; RF198702--201709",
    rights = "Downloaded directly from original public sources. Users remain responsible for source permissions; the package does not redistribute fund returns.")
  writeLines(capture.output(dput(panel$acquisition)), file.path(out_dir, "acquisition_provenance.txt"))
  if (!quiet) message("Raw inputs acquired; no model or result cache was read. Numeric checksums verified: ",
                       panel$acquisition$numeric_verified)
  panel
}
