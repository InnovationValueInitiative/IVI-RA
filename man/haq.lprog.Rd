% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{haq.lprog}
\alias{haq.lprog}
\title{Linear HAQ progression rates}
\format{A list with two elements.
\describe{
\item{tx}{A matrix with 5 columns (sname, est, se, lower, upper) for the name of the treatment,
the annual progression rate for that treatment, the standard error of the progression,
the .025\% quantile of the progression rate and the 97.5\% quantile of the progression rate,
respectively.}
\item{diff.age}{A matrix with 5 columns (age, est, se, lower, upper) for the 
age range (\eqn{age < 40}, \eqn{40 \le age < 65}, \eqn{ age \ge 65}),
the difference between progression rate for a given age range and the overall rate,
the standard error of the progression rate difference, the .025\% quantile of the progression
rate difference, and the 97.5\% quantile of the progression rate difference, respectively. }
}}
\source{
Michaud, Kaleb, Gene Wallenstein, and Frederick Wolfe. "Treatment and nontreatment
 predictors of health assessment questionnaire disability progression in 
 rheumatoid arthritis: a longitudinal study of 18,485 patients." 
 Arthritis care & research 63.3 (2011): 366-372.
 
Wolfe, Frederick, and Kaleb Michaud. "The loss of health status in rheumatoid 
arthritis and the effect of biologic therapy: a longitudinal observational study." 
Arthritis research & therapy 12.2 (2010): 1.
}
\usage{
haq.lprog
}
\description{
Constant annual rate of HAQ progression by age and treatment.
}
\keyword{datasets}
