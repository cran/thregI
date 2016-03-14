\name{plot.thregI}
\alias{plot.thregI}
\title{plot survival, hazard and density curves}
\description{
Display the graphs of the estimated survival, hazard or density function at different levels of a categorical variable which has been included in the threshold regression model by thregI(). There are three options, "sv", "hz" and "ds" are for survival, hazard and density function, respectively.
 }
\usage{
\method{plot}{thregI}(x,var,scenario,graph,nolegend=0,nocolor=0,...)
}
\arguments{
  \item{x}{a thregI object.}
  \item{var}{specifies the name of the variable which is required to be categorical. For each level of which (survival, hazard, density) the plots would be generated at given scenario specified by the scenario argument.}
  \item{scenario}{specifies a scenario for predicted plots.}
  \item{graph}{specifies the type of curves to be generated. The "sv" option is to plot survival function, the "hz" option is to plot hazard function and the "ds" option is to plot density function.}
  \item{nolegend}{set nolegend to be 1 when no need for legend. Users can add legends by themselves after set nolegend=1.}
 \item{nocolor}{set nocolor to be 1 if users would like to have all curves in black.}
 \item{\dots}{for future methods}
}

\examples{
#load the data "bcos"
data("bcos", package="thregI")

# transform the "treatment" variable into factor variable f.treatment
bcos$f.treatment=factor(bcos$treatment)

# fit the threshold regression model on the factor variable f.treatment
fit<-thregI(Surv(left, right, type='interval2')~f.treatment|f.treatment,data=bcos)

# estimated survival function at all levels of a categorical variable
plot.thregI(fit, var = f.treatment, graph = "sv", nocolor = 1)

# estimated hazard function at all levels of a categorical variable
plot.thregI(fit, var = f.treatment, graph = "hz", nocolor = 1)

# estimated density function at all levels of a categorical variable
plot.thregI(fit, var = f.treatment, graph = "ds", nocolor = 1)
}

\keyword{interval-censored data}
\keyword{plot}
\keyword{survival analysis}
\keyword{threshold regression}
