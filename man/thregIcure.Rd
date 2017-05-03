\name{thregIcure}
\alias{thregIcure}
\title{fit threshold regression cure-rate model for interval-censored data}
\description{
Fit a threshold regression cure-rate model for interval-censored data and calculate the prediction of the initial health status value, the drift value of the health process and cure rate. Whether the accommodated cure rate exists by goodness-of-fit test compared with the 'thregI'.
 }
\usage{
thregIcure(formula, data)
}

\arguments{
  \item{formula}{a formula object, which has the similar format to that in the thregI(). There are two | operators, on the left-hand side | operator is used for the linear regression function of \eqn{ln{y_0}} and specifing explaination variables of \eqn{\mu}. On the right of right-hand side | operator is used for specifing explaination variables of \eqn{logit (p)}. A constant \eqn{ln{y_0}}, \eqn{\mu} or \eqn{logit (p)} is allowed.}

  \item{data}{a data frame has the same format as that in the thregI().}
}

\examples{
#load the data "hdsd"
data("hdsd", package="thregI")

#transform the Noadyn variable into factor variable f.noadyn
hdsd$f.noadyn=factor(hdsd$Noadyn)

#fit the threshold regression cure-rate model
#the covariates are TR360, Noadyn, Sex and Age
fit<-thregIcure(Surv(left, right, type='interval2')~f.noadyn+TR360|f.noadyn+TR360|f.noadyn+TR360+
     Sex+Age, data=hdsd)
fit
#note that p-value of the goodness of fit test return by thregIcure function
}

\keyword{cure rate}
\keyword{interval-censored data}
\keyword{survival analysis}
\keyword{threshold regression}
