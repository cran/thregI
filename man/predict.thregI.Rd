\name{predict.thregI}
\alias{predict.thregI}
\title{predict the initial health status value and the drift value of the health process}
\description{Display the predicted values of threshold regression coefficients, density probability, survival probability, and hazard rate at the desired time and specified scenario; or, display the predicted values for all subjects corresponding to subjects scenario values at the desired time.
 }
\usage{
\method{predict}{thregI}(object,timevalue,scenario,...)
}

\arguments{
  \item{object}{a thregI object.}
  \item{timevalue}{specifies the desired time at which the predicted values of threshold regression coefficients, density probability, survival probability, and hazard rate.}
  \item{scenario}{specifies the values of all predictors considered in the thregI() and the dummy of the categorical variable of var. The predicted values at a specified time value for all subjects are calculated by the covariate values for each subject are used as their corresponding scenario values if ignore this argument.}
   \item{\dots}{for future methods}
}

\examples{
#load the data "bcos"
data("bcos", package="thregI")

# transform the treatment variable into factor variable f.treatment
bcos$f.treatment=factor(bcos$treatment)

# fit the threshold regression model on the factor variable f.treatment
fit<-thregI(Surv(left, right, type='interval2')~f.treatment|f.treatment,data=bcos)

# calculate the predicted values for y0, mu, f, s and h
# for the specified scenario that the treatment classification is 1 at time fifth month.
predict.thregI(fit, timevalue = 5, scenario=f.treatmentRadiation(1))
}

\keyword{interval-censored data}
\keyword{prediction}
\keyword{survival analysis}
\keyword{threshold regression}
