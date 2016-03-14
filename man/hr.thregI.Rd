\name{hr.thregI}
\alias{hr}
\alias{hr.thregI}
\title{perform hazard ratio for threshold regression model}
\description{
Returns the hazard ratios at a selected time for specified scenarios of the threshold regression model.
 }
\usage{
%hr(object,var,timevalue,scenario)
\method{hr}{thregI}(object,var,timevalue,scenario)
}
\arguments{
  \item{object}{a thregI object.}
  \item{var}{specifies the categorical variable which is required for the hazard ratios. Transform this variable into factor variable that has been used in thregI(). The lower level of the factor variable of var is the reference level of hazard ratio.}
  \item{timevalue}{specifies the desired time at which the hazard ratios would like to be had. A vector is allowed for this argument.}
  \item{scenario}{specifies the values of all variables that considered in the thregI() except variable of var.}
}

\examples{
#load the data "bcos"
data("bcos", package="thregI")

#transform the "treatment" variable into factor variable f.treatment
bcos$f.treatment=factor(bcos$treatment)

#fit the threshold regression model on the factor variable f.treatment
fit<-thregI(Surv(left, right, type='interval2')~f.treatment|f.treatment,data=bcos)

#hazard ratio of the radiation group vs. the radiation with chemotherapy at fifth month
hr.thregI(fit, var = f.treatment, timevalue = 5)
}

\keyword{hazard ratio}
\keyword{interval-censored data}
\keyword{survival analysis}
\keyword{threshold regression}
