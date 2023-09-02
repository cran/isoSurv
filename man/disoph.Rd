\name{disoph}
\alias{disoph}
\title{Fit Double Isotonic Proportional Hazards Model}
\description{Nonparametric estimation of monotone baseline hazard and monotone covariate effect functions in the proportional hazards model.}
\usage{
  disoph(formula, bshape, data, maxiter, eps)
}
\arguments{
  \item{formula}{formula object: response ~ iso(\eqn{z},shape="increasing")+\eqn{x_1+x_2+...+x_p}. The response must be right-censored survival outcome using the Surv function in the survival package. The iso function attributes the covariate \eqn{z}' name, shape and anchor point.}
  \item{bshape}{direnction of the baseline hazard function (\eqn{bshape}="increasing" or "decreasing").}
  \item{data}{data.frame includes variables named in the formula argument.}
  \item{maxiter}{maximum number of iteration (default is \eqn{10^4}).}
  \item{eps}{stopping convergence criteria (default is \eqn{10^-3}).}
}
\details{The disoph function computes (\eqn{\lambda0, \psi, \beta}) in the isotonic proportional hazards model, defined as
\deqn{\lambda(t|z,x)=\lambda0(t)exp(\psi(z)+\beta_1x_1+\beta_2x_2+...+\beta_px_p),}
based on the full likelihood, where \eqn{\lambda0} is a monotone increasing (or decreasing) baseline hazard function, \eqn{\psi} is a monotone increasing (or decreasing) covariate effect function, \eqn{z} is a univariate variable, (\eqn{x_1,x_2,...,x_p}) is a set of covariates, and \eqn{\beta=(\beta_1,\beta_2,...,\beta_p}) is a set of corresponding regression parameters. It allows to estimate \eqn{(\lambda0, \beta)} only if iso(\eqn{z},shape="increasing") is removed in the formula object. Likewise, It allows to estimate \eqn{(\lambda0, \psi)} only if \eqn{x} is removed in the formula object. Using the nonparametric maximum likelihood approaches, estimated \eqn{\lambda0} and \eqn{\psi} are right continuous increasing (or left continuos decreasing) step functions. Compared to the standard partial likelihood approach, the full likelihood approach in the disoph function additionally use shape-information on \eqn{\lambda0}, resulting in more efficient estimators especially for a finate sampe size.

For the anchor constraint, one point has to be fixed with \eqn{\psi(K)=0} to solve the identifiability problem, e.g. \eqn{\lambda0(t)exp(\psi(z))=(\lambda0(t)exp(-c))(exp(\psi(z)+c))} for any constant \eqn{c}. \eqn{K} is called an anchor point. By default, we set \eqn{K} as a median of values of \eqn{z}'s. The choice of anchor points are not important because, for example, different anchor points results in the same hazard ratios.
}
\value{A list of class fisoph:
\item{iso.bh}{data.frame with \eqn{t} and estimated \eqn{\lambda0}(t).}
\item{iso.cov}{data.frame with \eqn{z} and estimated \eqn{\psi}(z).}
\item{beta}{estimated \eqn{\beta_1,\beta_2,...,\beta_p}.}
\item{conv}{algorithm convergence status.}
\item{iter}{total number of iterations.}
\item{Zk}{anchor satisfying estimated \eqn{\psi(Zk)}=0.}
\item{shape.bh}{order restriction on \eqn{\lambda0}.}
\item{shape.cov}{order restriction on \eqn{\psi}.}
}

\references{Yunro Chung, Double Isotonic Proportional Hazards Models with Applications to Dose-Finding Studies. In preparation.}
\author{Yunro Chung [auth, cre]}

\examples{
#test1
test1=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4)
)
disoph.fit1=disoph(Surv(time, status)~iso(z,shape="inc"),bshape="inc",data=test1)
print(disoph.fit1)
plot(disoph.fit1)

#test2
test2=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4),
  trt=   c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0),
  x=     c(1, 2, 3, 4, 5, 6, 7, 8, 9, 1, 2, 3, 4, 5, 6)
)
disoph.fit2=disoph(Surv(time, status)~iso(z,shape="inc")+trt+x,bshape="inc",data=test2)
print(disoph.fit2)
plot(disoph.fit2)
}
