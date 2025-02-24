% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/orsf_control.R
\name{orsf_control_fast}
\alias{orsf_control_fast}
\title{Accelerated ORSF control}
\usage{
orsf_control_fast(method = "efron", do_scale = TRUE, ...)
}
\arguments{
\item{method}{(\emph{character}) a character string specifying the method
for tie handling. If there are no ties, all the methods are
equivalent. Valid options are 'breslow' and 'efron'. The Efron
approximation is the default because it is more accurate when dealing
with tied event times and has similar computational efficiency compared
to the Breslow method.}

\item{do_scale}{(\emph{logical}) if \code{TRUE}, values of predictors will be
scaled prior to each instance of Newton Raphson scoring, using summary
values from the data in the current node of the decision tree.}

\item{...}{Further arguments passed to or from other methods (not currently used).}
}
\value{
an object of class \code{'orsf_control'}, which should be used as
an input for the \code{control} argument of \link{orsf}.
}
\description{
Fast methods to identify linear combinations of predictors while
fitting an \link{orsf} model.
}
\details{
code from the  \href{https://github.com/therneau/survival/blob/master/src/coxfit6.c}{survival package}
was modified to make this routine.

Adjust \code{do_scale} \emph{at your own risk}. Setting \code{do_scale = FALSE} will
reduce computation time but will also make the \code{orsf} model dependent
on the scale of your data, which is why the default value is \code{TRUE}.
}
\seealso{
linear combination control functions
\code{\link{orsf_control}()},
\code{\link{orsf_control_cph}()},
\code{\link{orsf_control_custom}()},
\code{\link{orsf_control_net}()}
}
\concept{orsf_control}
