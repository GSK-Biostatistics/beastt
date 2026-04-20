# Bootstrap Covariate Data

Bootstrap Covariate Data

## Usage

``` r
bootstrap_cov(
  external_dat,
  n,
  imbal_var = NULL,
  imbal_prop = NULL,
  ref_val = 0
)
```

## Arguments

- external_dat:

  Data frame of the external data from which to bootstrap covariate
  vectors

- n:

  Number of rows in the output dataset

- imbal_var:

  Optional variable indicating which covariate's distribution should be
  altered to incorporate imbalance compared to the external data. If
  left `NULL`, the distributions of all covariates in the output dataset
  will match the distributions in the external dataset. The imbalance
  variable must be binary.

- imbal_prop:

  Optional imbalance proportion, required if an imbalance variable is
  specified. This defines the proportion of individuals with the
  reference value of the imbalance variable in the returned dataset.
  This can either be a single proportion or a vector of proportions, in
  which case a list of datasets is returned.

- ref_val:

  Optional value corresponding to the reference level of the binary
  imbalance variable, if specified

## Value

Data frame with the same number of columns as the external data frame
and n number of rows (if the length of `imbal_prop` is 0 or 1);
otherwise, a list of data frames with a length equal to that of
`imbal_prop`

## Details

Covariate data can be generated for `n` individuals enrolled in the
internal trial by bootstrap sampling entire covariate vectors from the
external data, thus preserving the correlation between the covariates.
If both `imbal_var` = `NULL` and `imbal_prop` = `NULL`, the function
returns a single data frame in which the distributions of each covariate
align with the covariate distributions from the external data (i.e.,
balanced covariate distributions across the two trials). Alternatively,
covariate imbalance can be incorporated into the generated sample with
respect to a binary covariate (`imbal_var`) such that a specified
proportion (`imbal_prop`) of individuals in the resulting sample will
have the reference level (`ref_val`) of this imbalance covariate. In
this case, stratified bootstrap sampling is employed with the imbalance
covariate as the stratification factor.

Multiple samples with varying degrees of imbalance can be generated
simultaneously by defining `imbal_prop` to be a vector of values. The
function then returns a list of data frames with a length equal to the
number of specified imbalance proportions.

## Examples

``` r
# Return one data frame with covariate distributions similar to external data
samp_balance <- bootstrap_cov(ex_binary_df, n = 1000)

# Return a list of two data frames that incorporate imbalance w.r.t. covariate 2
samp_imbalance <- bootstrap_cov(ex_binary_df, n = 1000, imbal_var = cov2,
                                imbal_prop = c(0.25, 0.5), ref_val = 0)
```
