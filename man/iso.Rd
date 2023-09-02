\name{iso}
\alias{iso}
\title{Attributions of isotonic covariate effect}
\description{attributes the covariate with respect to the name, direction, anchor point, and class.}
\usage{
  iso(z, shape)
}
\arguments{
  \item{z}{a univariate covariate.}
  \item{shape}{a direction of \eqn{z} (\eqn{shape}="increasing" or "decreasing").}
}
\details{Internal function. The iso function attributes the covariate \eqn{z} for its name, shape , anchor point and class, where the anchor point is set to a median of \eqn{z}'s, and class is set to "iso covariate"}
\value{The value \eqn{z} with attribution of its name, shape and median anchor point.}
\author{Yunro Chung [cre]}
