\name{plot.thregIcure}
\alias{plot.thregIcure}
\title{plot survival, hazard and density curves}

\description{
Display the graphs of the estimated survival, hazard or density function at different levels of a categorical variable which has been included in the threshold regression cure-rate model by thregIcure(). There are three options, "sv", "hz" and "ds" are for survival, hazard and density function, respectively.
 }

\usage{
\method{plot}{thregIcure}(x,var,scenario,graph,nolegend=0,nocolor=0,...)
      }

\arguments{
  \item{x}{a thregIcure object.}
  \item{var}{specifies the name of the variable which is required to be categorical. The use of the var argument is the same as that in the plot.thregI().}
  \item{scenario}{specifies a scenario for predicted plots.}
  \item{graph}{specifies the type of curves to be generated. The "hz" option is to plot hazard function accommodated a cure rate, the "sv" option is to plot survival function accommodated a cure rate and the "ds" option is to plot density function accommodated a cure rate.
}
  \item{nolegend}{ The use of the nolegend argument is the same as that in the plot.thregI().}
 \item{nocolor}{The use of the nolegend argument is the same as that in the plot.thregI().}
 \item{\dots}{for future methods}
}


\examples{
#load the data "hdsd"
data("hdsd", package="thregI")

#transform categorical variable Noadyn into factor variable f.noadyn
hdsd$f.noadyn=factor(hdsd$Noadyn)

#fit the threshold regression cure-rate model
#the covariates are TR360, Noadyn, Sex and Age
fit<-thregIcure(Surv(left, right, type='interval2')~f.noadyn|TR360|Sex+Age, data=hdsd)

#plot estimated hazard, survival and density function
#subject is ambulatory (Noadyn1=1), TR360 = 1.5, male and 30 years old
plot.thregIcure(fit, var=f.noadyn, scenario=TR360(1.5)+Sex(1)+Age(30), graph = "sv", nocolor = 1)
plot.thregIcure(fit, var=f.noadyn, scenario=TR360(1.5)+Sex(1)+Age(30), graph = "hz", nocolor = 1)
plot.thregIcure(fit, var=f.noadyn, scenario=TR360(1.5)+Sex(1)+Age(30), graph = "ds", nocolor = 1)
}

\keyword{cure rate}
\keyword{interval-censored data}
\keyword{plot}
\keyword{survival analysis}
\keyword{threshold regression}
