% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ps.R
\name{prop_scr_hist}
\alias{prop_scr_hist}
\title{Histogram of the Propensity Score Object}
\usage{
prop_scr_hist(
  x,
  variable = c("propensity score", "ps", "inverse probablity weight", "ipw"),
  ...
)
}
\arguments{
\item{x}{Propensity score object}

\item{variable}{Variable to plot, it must be either a propensity score or
inverse probability weight TODO MAKE THIS CLEARER}

\item{...}{Optional options for \code{geom_histogram}}
}
\value{
ggplot object
}
\description{
Histogram of the Propensity Score Object
}
