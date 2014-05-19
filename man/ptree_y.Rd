\name{ptree_y}
\alias{ptree_y}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ptree_y
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptree_y(newtree, node_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newtree}{
%%     ~~Describe \code{newtree} here~~
}
  \item{node_id}{
%%     ~~Describe \code{node_id} here~~
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
function (newtree, node_id) 
{
    p <- nodes(newtree, node_id)[[1]]$prediction
    if (length(p) == 2) {
        return(p[2])
    }
    return(p)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
