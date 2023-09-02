\name{isoph}
\alias{isoph}
\title{Fit Isotonic Proportional Hazards Model}
\description{Nonparametric estimation of a monotone covariate effect under the proportional hazards model.}
\usage{
  isoph(formula, data, maxiter, eps)
}
\arguments{
  \item{formula}{a formula object: response ~ iso(\eqn{z},shape="increasing")+\eqn{x_1+x_2+...+x_p}. The response must be right-censored survival outcome using the Surv function in the survival package. The iso function attributes the covariate \eqn{z}' name, shape and anchor point.}
  \item{data}{data.frame includes variables named in the formula argument.}
  \item{maxiter}{maximum number of iteration (default is \eqn{10^4}).}
  \item{eps}{stopping convergence criteria (default is \eqn{10^-3}).}
}
\details{The isoph function estimates (\eqn{\psi}, \eqn{\beta}) in the isotonic proportional hazards model, defined as
\deqn{\lambda(t|z,x)=\lambda0(t)exp(\psi(z)+\beta_1x_1+\beta_2x_2+...+\beta_px_p),}
based on the partial likelihood with unspecified baseline hazard function \eqn{\lambda0}, where \eqn{\psi} is a monotone increasing (or decreasing) covariate effect function, \eqn{z} is a univariate variable, \eqn{x=(x_1,x_2,...,x_p)} is a set of covariates, and \eqn{\beta=(\beta_1,\beta_2,...,\beta_p)} is a set of corresponding regression parameters. It allows to estimate \eqn{\psi} only if \eqn{x} is removed in the formula object. Using the nonparametric maximum likelihood approaches, estimated \eqn{\psi} is a right continuous increasing (or left continuos decreasing) step function.

For the anchor constraint, one point has to be fixed with \eqn{\psi(K)=0} to solve the identifiability problem, e.g. \eqn{\lambda0(t)exp(\psi(z))=(\lambda0(t)exp(-c))(exp(\psi(z)+c))} for any constant \eqn{c}. \eqn{K} is called an anchor point. By default, we set \eqn{K} as a median of values of \eqn{z}'s. The choice of anchor points are not important because, for example, different anchor points results in the same hazard ratios.}
\value{A list of class isoph:
\item{iso.cov}{data.frame with \eqn{z} and estimated \eqn{\psi}.}
\item{beta}{estimated \eqn{\beta_1,\beta_2,...,\beta_p}.}
\item{conv}{algorithm convergence status.}
\item{iter}{total number of iterations.}
\item{Zk}{anchor point satisfying \eqn{\psi(Zk)}=0.}
\item{shape}{Order-restriction imposed on \eqn{\psi}.}
}

\references{Yunro Chung, Anastasia Ivanova, Michael G. Hudgens, Jason P. Fine, Partial likelihood estimation of isotonic proportional hazards models, Biometrika. 2018, 105 (1), 133-148. doi:10.1093/biomet/asx064}
\author{Yunro Chung [aut, cre]}

\examples{
# test1
test1=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4)
)
isoph.fit1=isoph(Surv(time, status)~iso(z,shape="inc"),data=test1)
print(isoph.fit1)
plot(isoph.fit1)

# test2
test2=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4),
  trt=   c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0)
)
isoph.fit2=isoph(Surv(time, status)~iso(z,shape="inc")+trt, data=test2)
print(isoph.fit2)
plot(isoph.fit2)
}

