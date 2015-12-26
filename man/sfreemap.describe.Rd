\name{sfreemapr.describe}
\alias{sfreemapr.describe}
\title{describe a \code{phylo} object modified by sfreemapr.map}
\usage{
sfreemapr.describe(tree)
}
\arguments{
    \item{tree}{A result of a sfreemapr.map execution (\code{phylo} or \code{multiPhylo})}

}
\description{
    This function returns the mean number of transitions and dwelling times of tree (which can be multiple trees as a multiPhylo object).
}
%\details{
%    Details
%}
\value{
    A named list with values \code{transitions} and \code{dwelling_times}.
}

%\references{
%    \code{\link{sfreemapr.map}}
%}

\author{Diego Pasqualin \email{dpasqualin@inf.ufpr.br}}
\seealso{
    \code{\link{sfreemapr.map}}
}
\keyword{utilities}
%\keyword{simulation}
