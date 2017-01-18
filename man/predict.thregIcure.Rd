\name{predict.thregIcure}
\alias{predict.thregIcure}
\title{predict the initial health status value, the drift value of the health process and cure rate}
\description{
Display the predicted values of threshold regression coefficients, non-cure rate, density probability, survival probability and hazard rate at the desired time and specified scenario; or, display the predicted values for all subjects corresponding to subjects scenario values at the desired time.
 }
\usage{
\method{predict}{thregIcure}(object,timevalue,scenario,...)
      }
\arguments{
  \item{object}{a thregIcure object.}
  \item{timevalue}{specifies the desired time at which the predicted values of threshold regression coefficients, non-cure rate, density probability, survival probability, and hazard rate.}
  \item{scenario}{specifies the values of all predictors considered in the thregIcure() and the dummy of the categorical variable of var. The use of the scenario argument is the same as that in the predit.thregI().}
  \item{\dots}{for future methods}
}

\examples{
#load the data "hdsd"
data("hdsd", package="thregI")

# transform categorical variable Noadyn into factor variable f.noadyn
hdsd$f.noadyn=factor(hdsd$Noadyn)

# fit the threshold regression cure-rate model
# the covariates are TR360, Noadyn, Sex and Age
fit<-thregIcure(Surv(left, right, type='interval2')~f.noadyn|TR360|Sex+Age, data=hdsd)

# calculate the predicted values
# subject is ambulatory (Noadyn1=1), TR360 = 1.5, male and 30 years old
# note that f.noadyn1 is the name of non-reference group
predict.thregIcure(fit, timevalue = 1, scenario=f.noadyn1(1)+TR360(1.5)+Sex(1)+Age(30))
}
\keyword{cure rate}
\keyword{interval-censored data}
\keyword{prediction}
\keyword{survival analysis}
\keyword{threshold regression}

