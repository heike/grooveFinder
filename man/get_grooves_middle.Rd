% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_grooves_middle.R
\name{get_grooves_middle}
\alias{get_grooves_middle}
\title{Use the center of a crosscut}
\usage{
get_grooves_middle(x, value, middle = 75, return_plot = F)
}
\arguments{
\item{x}{numeric vector of locations in microns}

\item{value}{numeric vector of surface measurements in microns}

\item{middle}{middle percent to use for the identification}

\item{return_plot}{return plot?}
}
\value{
list of groove vector and plot of crosscut, if return_plot is true
}
\description{
Use the center of a crosscut
}
\examples{
\dontrun{
# Set the data up to be read in, cleaned, etc.
library(bulletxtrctr)
library(x3ptools)

example_data <- bullet_pipeline(
  location = list(Bullet1 = c(hamby252demo$bullet1[3])),
  stop_at_step = "crosscut",
  x3p_clean = function(x) x \%>\%
      x3p_scale_unit(scale_by = 10^6) \%>\%
      rotate_x3p(angle = -90) \%>\%
      y_flip_x3p()
)

get_grooves_middle(example_data$ccdata[[1]]$x,
  example_data$ccdata[[1]]$value,
  return_plot = T
)
cc_locate_grooves(example_data$ccdata[[1]],
  method = "middle",
  return_plot = T
)
}
}
