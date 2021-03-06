% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pop.R
\name{sample_rtmvnorm_da}
\alias{sample_rtmvnorm_da}
\title{Sample disease activity levels and HAQ}
\usage{
sample_rtmvnorm_da(n, mean, sigma, lower, upper)
}
\arguments{
\item{n}{Number of samples.}

\item{mean}{Mean vector.}

\item{sigma}{Covariance matrix.}

\item{lower}{Vector of lower truncation points.}

\item{upper}{Vector of upper truncation points.}
}
\description{
Sample disease activity levels and HAQ from a truncated multivariate normal distribution.
}
\keyword{internal}
