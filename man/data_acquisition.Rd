\name{fetch_empirical_data}
\alias{fetch_empirical_data}
\title{Acquire raw empirical inputs from the original public sources}
\description{On an explicit user call, downloads the original authors' public
code/data archive and a pinned historical French RF archive, extracts the article
sample, and writes input CSVs. No raw fund returns are distributed by this package.}
\usage{
fetch_empirical_data(out_dir, author_archive = NULL, rf_archive = NULL,
  author_url = paste0("https://www.dropbox.com/scl/fo/v9qgrwn5xgx7ahwdyqhki/",
    "AFJ9_QuPz2BrHEXfl29BQlU?rlkey=sn2za70nnwojcw8jnogw34uhs&dl=1"),
  rf_url = paste0("https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/",
    "Historical_Archives/08\%202018\%20Update/ftp/F-F_Research_Data_Factors_CSV.zip"),
  verify = TRUE, quiet = FALSE)
}
\arguments{
\item{out_dir}{New or empty directory for \code{fund.csv}, \code{factor.csv}, \code{rf.csv} and source provenance. Nonempty directories are not overwritten.}
\item{author_archive}{Optional already downloaded author ZIP. Supply with \code{rf_archive} for fully offline processing.}
\item{rf_archive}{Optional already downloaded July 2018 French factor ZIP.}
\item{author_url}{Public archive linked from original coauthor Bin Peng's research page. Nothing from the archive is executed.}
\item{rf_url}{Pinned official historical July 2018 data cut, released in August 2018.}
\item{verify}{Require canonical numeric checksums to match the article's inputs. Default TRUE. Set FALSE only for an explicitly investigated source revision; changed data will not reproduce the article exactly.}
\item{quiet}{Suppress download and completion messages.}
}
\details{The source is Gao, Peng and Yan's (2024) author archive linked at
\url{https://sites.google.com/view/binpengau/research-related}. Its
\code{CaseStudy/mutual_fund.mat} contains \code{R} with 668 dates and 726 funds and
\code{F} with 668 dates and four factors. The helper retains rows 301--668,
corresponding to February 1987--September 2017, transposes returns, and keeps all
726 funds. It does not apply the separate original case study's less-than-20-percent
missingness screen. An intentionally narrow base-R reader supports only the
archive's little-endian MATLAB v5 ordinary real-double matrices (including
compressed blocks); no MATLAB installation or contributed R dependency is needed.
Unrelated archive scripts are never extracted or executed.

RF is taken from the official French July 2018 historical archive, not the
current factor download. In the checked current download, May 1994 RF is
0.32 percent rather than the article input's 0.31 percent. The pinned vintage
matches all 368 original RF observations, and the author's MAT arrays match
all original fund/factor observations. Canonical little-endian double-vector
checksums are checked before any input CSVs are written, independently of ZIP
timestamps, text formatting or newline conventions. CSV export uses 17 significant
digits. RF remains in percentage units for \code{read_fund_panel}'s default.

Downloading a public archive does not transfer redistribution rights. Users are
responsible for the terms of the original authors, CRSP-derived data and French
Data Library. The package hosts neither these raw returns nor a replacement
download mirror. URLs and contents can change; strict checks reject a silent
source substitution. Both archive arguments can be supplied for offline use.}
\value{A panel compatible with \code{read_fund_panel}, with an additional
\code{acquisition} component containing URLs, source-file hashes, numeric hashes,
transformation and vintage metadata. The three raw input CSVs and an
\code{acquisition_provenance.rds} file are written to \code{out_dir}.}
\examples{
\dontrun{
data <- fetch_empirical_data("raw-inputs")
fit <- reproduce_empirical(data, cores = 25)
calibration <- reproduce_calibration(data, replications = 1000, cores = 25)
}
}
\seealso{\code{read_fund_panel}, \code{reproduce_empirical}, \code{reproduce_calibration}}
