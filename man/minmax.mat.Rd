\name{minmax.mat}
\alias{minmax.mat}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
minmax.mat
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
minmax.mat(str, varnms)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{str}{
%%     ~~Describe \code{str} here~~
}
  \item{varnms}{
%%     ~~Describe \code{varnms} here~~
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
function (str, varnms) 
{
    comps <- strsplit(str, ",")
    MMM <- matrix(data = rep(c(-Inf, Inf), length(varnms)), nrow = length(varnms), 
        ncol = 2, byrow = TRUE)
    rownames(MMM) <- varnms
    for (i in 1:(length(comps[[1]]) - 1)) {
        nodestr <- strsplit(trim(comps[[1]][i]), " ")
        node.varnm <- trim(nodestr[[1]][1])
        node.dir <- trim(nodestr[[1]][2])
        node.split <- trim(nodestr[[1]][3])
        var.row <- which(varnms == node.varnm)
        if (node.dir == "<=") {
            MMM[var.row, 2] <- as.numeric(node.split)
        }
        else {
            MMM[var.row, 1] <- as.numeric(node.split)
        }
    }
    y <- comps[[1]][length(comps[[1]])]
    return(list(M = MMM, y = y))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
