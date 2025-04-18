% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps.R
\name{trim_ps}
\alias{trim_ps}
\title{Trim a \code{prop_scr} object}
\usage{
trim_ps(x, low = NULL, high = NULL, quantile = FALSE)
}
\arguments{
\item{x}{A \code{prop_scr} object}

\item{low}{Low cut-off such that all participants with propensity scores less
than this value (or quantile if \code{quantile = TRUE}) are removed.  If left
\code{NULL} no lower bound will be used}

\item{high}{High cut-off such that all participants with propensity scores
greater than this value (or quantile if \code{quantile = TRUE}) are removed. If
left \code{NULL} no upper bound will be used}

\item{quantile}{True/False value to determine if the cut-off values are based
directly on the propensity scores (false) or their quantiles (true). By default this is
false.}
}
\value{
a \code{prop_scr} object with a trimmed propensity score distribution
}
\description{
Trim a \code{prop_scr} object
}
\details{
This function uses R's default method of quantile calculation (type
7)
}
\examples{
library(dplyr)
ps_obj <- calc_prop_scr(internal_df = filter(int_binary_df, trt == 0),
                       external_df = ex_binary_df,
                       id_col = subjid,
                       model = ~ cov1 + cov2 + cov3 + cov4)
trim_ps(ps_obj, low = 0.3, high = 0.7)

}
