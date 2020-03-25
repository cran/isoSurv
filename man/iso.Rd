\name{iso}
\alias{iso}
\title{Attributions of isotonic covariate effect}
\description{InternalIt attributes the covariate with respect to the name, direction and anchor constrain.}
\usage{
  iso(z, shape, K)
}
\arguments{
  \item{z}{a univariate covariate.}
  \item{shape}{a direction of z (shape="increasing" or "decreasing") (default is \eqn{"increasing"}).}
  \item{K}{an anchor point (default is \eqn{NA} for the median anchor \eqn{K=median(z)}).}
}
\details{Internal function. The iso function attributes the covariate z for its name, shape direction and anchor.}
\value{The value \eqn{z} with attribution of its name, shape and \eqn{K}.}
\author{Yunro Chung [cre]}
