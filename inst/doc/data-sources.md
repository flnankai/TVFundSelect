# Raw empirical input acquisition

The package does not contain fund-level returns or the unpublished manuscript.
`fetch_empirical_data("data-private")` is an explicit, optional network action:
it downloads original-source raw data, transforms their orientation/units and
checks the article input numerically. It never downloads the present authors'
saved selection decisions, fitted models or result summaries, and never runs
any external MATLAB code.

## Original-source fund and factor archive

Gao, J., Peng, B. and Yan, Y. (2024), *Higher-Order Expansions and Inference for
Panel Data Models*, JASA 119, 2760-2771, DOI 10.1080/01621459.2023.2277411.
Bin Peng's author page links the public MATLAB replication archive:

- Author page: https://sites.google.com/view/binpengau/research-related
- Author archive: https://www.dropbox.com/scl/fo/v9qgrwn5xgx7ahwdyqhki/AFJ9_QuPz2BrHEXfl29BQlU?rlkey=sn2za70nnwojcw8jnogw34uhs&dl=1

Only `CaseStudy/mutual_fund.mat` is read. It contains the four factors F
(668 x 4) and net returns R (668 x 726). The final 368 months are retained and
R transposed to 726 x 368. Do NOT apply the source script's further missingness
screen: the present manuscript defines its own longest-history, balanced-panel
and training-window samples. The final fund and factor matrices agree exactly
with the article's archived CSV inputs. A deliberately narrow base-R MAT-v5
numeric reader avoids requiring MATLAB, Python or another R dependency; it is
not a general-purpose MAT-file import library and rejects unsupported objects.

The underlying fund returns originate in CRSP. Public download availability
does not itself establish unrestricted redistribution permission. Users are
responsible for their access/use permissions. We link the source; we do not
rehost those fund observations under the software license.

## Historical risk-free vintage

The one-month RF series is obtained from Kenneth French's historical archive:

https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/Data_Library/Historical_Archives/08%202018%20Update/ftp/F-F_Research_Data_Factors_CSV.zip

The selected monthly dates are 198702 through 201709. This July-2018 archive
has exactly the article's RF values. A current download differs in May 1994
(0.32 percent versus 0.31 percent in the article/historical archive), so the
helper deliberately pins the historical vintage and validates numeric values.
The factor series are retained from the authors' archive, not replaced by a
current French download.

## Offline inputs

- `fund.csv`: 726 rows x 368 columns, decimals, no header; zero means missing.
- `factor.csv`: 368 x 4, decimal Mkt-RF, SMB, HML and MOM; no intercept/header.
- `rf.csv`: 368 monthly RF values in percent, no header; converted once by /100.
- Dates: February 1987 through September 2017. Fund IDs are row labels only.

Reference original-file SHA256 (byte hashes, not a requirement on alternative
CSV formatting of the same numeric matrix):

| File | SHA256 |
|---|---|
| fund.csv | 507253A1020C97E52EE72C0A825DC3AFF607EF4C7B1B6B2DE3BC8702B68E24A5 |
| factor.csv | 219EC575B45040708BBEF62CE2A5802ECA73A8216211D3E9ABC40741FCAFDD60 |
| rf.csv | 8112E0F7731BEB6DFAE8E0DBCED78A247C2BCA6BCA864D5169D396609E6ED159 |

The fetch function additionally verifies canonical numeric checksums so that
harmless CSV formatting differences do not conceal changed observations. Its
output records URLs, transformation, numeric verification and input-file hashes.
If an upstream file changes, verification fails rather than silently using it.
