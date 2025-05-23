% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/inspre.R
\name{est_noise_var}
\alias{est_noise_var}
\title{Estimates noise variances given graph and node variances.}
\usage{
est_noise_var(H, v, non_neg = FALSE)
}
\arguments{
\item{H}{D x D matrix.}

\item{v}{D-vector of floats. Total variance of each node.}

\item{non_neg}{Boolean. True to use non-negative least squares instead
of exact solution. Default FALSE.}
}
\description{
Can use H = (I-G)^-1 or an approximation H = I + G + G%^%2 + ...
}
