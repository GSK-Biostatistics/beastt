% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tte-sim.R
\name{calc_study_duration}
\alias{calc_study_duration}
\title{Calculate the analysis time based on a target number of events}
\usage{
calc_study_duration(
  study_time,
  observed_time,
  event_indicator,
  target_events = NULL,
  target_follow_up = NULL
)
}
\arguments{
\item{study_time}{Vector of study (accrual + observed) times}

\item{observed_time}{Vector of observed times}

\item{event_indicator}{Vector of boolean values (True/False or 1/0) indicating if the observed time value is an event or censoring}

\item{target_events}{Number of target events, if only using target follow-up time leave NULL}

\item{target_follow_up}{Target follow-up for each subject, if only using target events leave NULL}
}
\value{
Time of analysis
}
\description{
Calculate the analysis time based on a target number of events
}
\examples{
library(dplyr)
# Determining analysis time by reaching a number of events
ex_tte_df |> mutate(
  analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
                                     event_indicator = event, target_events = 30)
)
# Determining analysis time by minimum follow-up time
ex_tte_df |> mutate(
  analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
                                     event_indicator = event, target_follow_up = 12)
)
# Or use both and whichever happens first
ex_tte_df |> mutate(
  analysis_time = calc_study_duration(study_time = total_time, observed_time = y,
                                     event_indicator = event,
                                     target_events = 30, target_follow_up = 12)
)
}
