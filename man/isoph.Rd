\name{isoph}
\alias{isoph}
\title{Fit Isotonic Proportional Hazards Model}
\description{Nonparametric partial likelihood estimation of a monotone covariate effect in the proportional hazards model.}
\usage{
  isoph(formula, data, maxiter, eps)
}
\arguments{
  \item{formula}{a formula object: response ~ iso(covariate1)+covariate2+.... The response must be survival outcome using the Surv function in the survival package. The iso(covariate1,shape="increasing",K="median") is for isotonic estimation of covariate1 with shape="increasing" (or "decreasing") direction and K="median" anchor.}
  \item{data}{data.frame includes variables named in the formula argument.}
  \item{maxiter}{maximum number of iteration (default is \eqn{10^4}).}
  \item{eps}{stopping convergence criteria (default is \eqn{10^-3}).}
}
\details{The isoph function estimates \eqn{\psi} and \eqn{\beta} based on the partial likelihood under the isotonic proportional hazards model, defined as
\deqn{\lambda(t|z,w)=\lambda0(t)exp(\psi(z)+\beta w),}
where \eqn{\lambda0} is a baseline hazard function, \eqn{\psi} is an isotonic function, \eqn{z} is a univariate variable, \eqn{w} is a q by 1 dimensional covariate vector and \eqn{\beta} is a 1 by q dimensional regression parameter. One point has to be fixed with \eqn{\psi(K)=0}, where \eqn{K} is an anchor point. A direction of \eqn{\psi} is defined as monotone increasing or monotone decreasing in \eqn{z}. Pseudo iterative convex minorant algorithm is used to maximize the partial likelihood.}
\value{A list of class isoph:
\item{iso.cov}{data.frame with \eqn{z} and estimated \eqn{\psi}.}
\item{beta}{estimated \eqn{\beta}.}
\item{conv}{algorithm convergence status.}
\item{iter}{total number of iterations.}
\item{Zk}{anchor satisfying estimated \eqn{\psi(Zk)}=0. See the formular argument above.}
\item{shape}{order restriction on \eqn{\psi}. See the formular argument above.}
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

isoph.fit1 = isoph(Surv(time, status)~iso(z,shape="inc"),data=test1)
plot(isoph.fit1)

# test2
test2=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 1, 0, 1, 0, 1, 1, 1, 0, 1, 1, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4),
  trt=   c(1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 0, 1, 1, 0, 0)
)

isoph.fit2 = isoph(Surv(time, status)~iso(z,shape="inc")+trt, data=test2)
plot(isoph.fit2)
}

