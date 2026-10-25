library(TVFundSelect)
reader <- TVFundSelect:::.gpy_mat_numeric
must_error <- function(x) stopifnot(inherits(tryCatch(force(x), error = identity), "error"))
int_raw <- function(x) writeBin(as.integer(x), raw(), size = 4L, endian = "little")
element <- function(kind, value) c(int_raw(c(kind, length(value))), value,
                                  raw((8L - length(value) %% 8L) %% 8L))
matrix_element <- function(name, x, compressed = FALSE, flags = 6L) {
  stopifnot(nchar(name) == 1L)
  # Small data-element name tag, as used by the actual one-letter F/R names.
  name_tag <- c(int_raw(1L + 65536L), charToRaw(name), raw(3L))
  body <- c(element(6L, int_raw(c(flags, 0L))), element(5L, int_raw(dim(x))),
            name_tag, element(9L, writeBin(as.double(x), raw(), size = 8L, endian = "little")))
  block <- element(14L, body)
  if (compressed) {
    z <- memCompress(block, type = "gzip")
    block <- c(int_raw(c(15L, length(z))), z) # Compressed elements have no added padding.
  }
  block
}
header <- raw(128L); title <- charToRaw("MATLAB 5.0 MAT-file synthetic fixture")
header[seq_along(title)] <- title; header[125:126] <- as.raw(c(0, 1)); header[127:128] <- charToRaw("IM")
R <- matrix(seq_len(12L) / 10, 4L, 3L); F <- matrix(seq_len(8L) / 7, 4L, 2L)
path <- tempfile(fileext = ".mat")
writeBin(c(header, matrix_element("R", R), matrix_element("F", F)), path)
plain <- reader(path)
stopifnot(identical(plain$R, R), identical(plain$F, F))
writeBin(c(header, matrix_element("R", R, TRUE), matrix_element("F", F, TRUE)), path)
compressed <- reader(path)
stopifnot(identical(compressed, plain))
bad_header <- header; bad_header[127:128] <- charToRaw("MI")
writeBin(c(bad_header, matrix_element("R", R), matrix_element("F", F)), path)
must_error(reader(path))
writeBin(c(header, matrix_element("R", R, flags = 1L), matrix_element("F", F)), path)
must_error(reader(path))
writeBin(c(header, matrix_element("R", R, flags = 2054L), matrix_element("F", F)), path)
must_error(reader(path))
writeBin(c(header, matrix_element("R", R)), path)
must_error(reader(path))
writeBin(c(header, as.raw(1:3)), path)
must_error(reader(path))

# A self-generated complete French-style monthly CSV verifies date parsing,
# exclusion of annual rows, and use of the fifth (RF) field, with no network.
dates <- format(seq(as.Date("1987-02-01"), by = "month", length.out = 368L), "%Y%m")
rf <- seq_len(368L) / 100
lines <- c("Synthetic factor description", ",Mkt-RF,SMB,HML,RF",
           paste(dates, 1, 2, 3, rf, sep = ","), "Annual Factors: January-December", "1987,1,2,3,4")
csv <- tempfile(fileext = ".csv"); writeLines(lines, csv)
stopifnot(identical(TVFundSelect:::.gpy_rf_csv(csv), rf))
writeLines(lines[-5L], csv); must_error(TVFundSelect:::.gpy_rf_csv(csv))
stopifnot(identical(TVFundSelect:::.gpy_numeric_md5(c(.31, .32)),
                   TVFundSelect:::.gpy_numeric_md5(matrix(c(.31, .32), 2L, 1L))))
cat("Offline data acquisition tests passed: plain/compressed MATLAB v5 numeric arrays, small tags, unsupported-format rejection and dated RF parsing.\n")
