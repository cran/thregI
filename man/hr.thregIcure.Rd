\name{hr.thregIcure}
%\alias{hr}
\alias{hr.thregIcure}
\title{perform hazard ratio for threshold regression cure-rate model}

\description{
Returns the hazard ratios at a selected time for specified scenarios of the threshold regression cure-rate model.
 }
\usage{
%hr(object,var,timevalue,scenario)
\method{hr}{thregIcure}(object,var,timevalue,scenario)
}

\arguments{
  \item{object}{a thregIcure object.}
  \item{var}{specifies the name of the variable which is required to be categorical for which the hazard ratios are to be calculated. The use of the var argument is similar to that in the hr.thregI().}
  \item{timevalue}{specifies the desired time at which the hazard ratios would like to be had. A vector is allowed for this argument.}
  \item{scenario}{  specifies the values of all variables that considered in the thregIcure function except variable of var.}
}

\examples{
#load the data "hdsd"
data("hdsd", package="thregI")

#transform categorical variable Noadyn into factor variable f.noadyn
hdsd$f.noadyn=factor(hdsd$Noadyn)

#fit the threshold regression cure-rate model
#the covariates are TR360, Noadyn, Sex and Age
fit<-thregIcure(Surv(Left, Right, type='interval2')~f.noadyn|TR360|Sex+Age, data=hdsd)

#calculate the hazard ratios for categorical variable Noadyn
#specified scenario of TR360=1.5, Male, Age=30 at time one
hr.thregIcure(fit, var = f.noadyn, timevalue = 1, scenario=TR360(1.5)+Sex(1)+Age(30))
}

\keyword{cure rate}
\keyword{hazard ratio}
\keyword{interval-censored data}
\keyword{survival analysis}
\keyword{threshold regression}
