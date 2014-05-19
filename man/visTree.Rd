\name{visTree}
\alias{visTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
visTree
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
visTree(cond.tree, X, Y, range = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{cond.tree}{
%%     ~~Describe \code{cond.tree} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
}
  \item{range}{
%%     ~~Describe \code{range} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (cond.tree, X, Y, range = NULL) 
{
    splittree <- list_node(cond.tree)
    structure <- strsplit(list_node(cond.tree), split = ";")
    if (length(unlist(structure)) == 1) {
        stop("Tree has only a single node; nothing to visualize.")
    }
    n.terminals <- ifelse(is.null(range), length(unlist(structure)), 
        length(range))
    if (is.null(range)) {
        index <- 1:n.terminals
    }
    else {
        index <- min(range):min(max(range), length(unlist(structure)))
    }
    par(mfrow = c(2, ceiling(length(index)/2)), mar = c(2, 1, 
        3, 1))
    sapply(unlist(structure)[index], function(S) {
        plot.minmax(minmax.mat(S, colnames(X)), X, Y)
    })
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
