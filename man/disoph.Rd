\name{disoph}
\alias{disoph}
\title{Fit Double Isotonic Proportional Hazards Model}
\description{Nonparametric full likleihood estimation of monotone baseline hazard and covariate effect functions in the proportional hazards model.}
\usage{
  disoph(formula, bshape, data, maxiter, eps)
}
\arguments{
  \item{formula}{a formula object: response ~ iso(covariate1). The response must be survival outcome using the Surv function in the survival package. The iso(covariate1,shape="increasing",K="median") is for isotonic estimation of covariate1 with shape="increasing" (or "decreasing") direction and K="median" anchor.}
  \item{bshape}{direnction of the baseline hazard function (bshape="increasing" or "decreasing").}
  \item{data}{data.frame includes variables named in the formula argument.}
  \item{maxiter}{maximum number of iteration (default is \eqn{10^4}).}
  \item{eps}{stopping convergence criteria (default is \eqn{10^-3}).}
}
\details{The disoph function estimates \eqn{\lambda0} and \eqn{\psi} based on the full likelihood under the isotonic proportional hazards model, defined as
\deqn{\lambda(t|z)=\lambda0(t)exp(\psi(z)),}
where \eqn{\lambda0} is an isotonic baseline hazard function and \eqn{\psi} is an isotonic function. One point has to be fixed with \eqn{ \psi(K)=0 }, where \eqn{K} is an anchor point. A direction of \eqn{\lambda0} is defined as monotone increasing or monotone decreasing in \eqn{t}. A direction of \eqn{\psi} is defined as monotone increasing or monotone decreasing in \eqn{z}. Back-and-forth iterative pool adjacent violators algorithm is used to maximize the full likelihood.}
\value{A list of class fisoph:
\item{iso.bh}{data.frame with \eqn{time} and estimated \eqn{\lambda0}.}
\item{iso.cov}{data.frame with \eqn{z} and estimated \eqn{\psi}.}
\item{conv}{algorithm convergence status.}
\item{iter}{total number of iterations.}
\item{Zk}{anchor satisfying estimated \eqn{\psi(Zk)}=0.}
\item{shape.bh}{order restriction on \eqn{\lambda0}.}
\item{shape.cov}{order restriction on \eqn{\psi}.}
}

\references{Yunro Chung, Double Isotonic Proportional Hazards Models with Monotone Baseline Hazard and Covariate Effect Functions. In preparation.}
\author{Yunro Chung [auth, cre]}

\examples{
test1=data.frame(
  time=  c(2, 5, 1, 7, 9, 5, 3, 6, 8, 9, 7, 4, 5, 2, 8),
  status=c(0, 1, 0, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0, 1),
  z=     c(2, 1, 1, 3, 5, 6, 7, 9, 3, 0, 2, 7, 3, 9, 4)
)

disoph.fit1=disoph(Surv(time, status)~iso(z,shape="inc"),bshape="inc",data=test1)
plot(disoph.fit1)
}
