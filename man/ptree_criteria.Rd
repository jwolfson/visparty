\name{ptree_criteria}
\alias{ptree_criteria}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
ptree_criteria
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ptree_criteria(newtree, node_id, left)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newtree}{
%%     ~~Describe \code{newtree} here~~
}
  \item{node_id}{
%%     ~~Describe \code{node_id} here~~
}
  \item{left}{
%%     ~~Describe \code{left} here~~
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
function (newtree, node_id, left) 
{
    if (nodes(newtree, node_id)[[1]]$terminal) {
        return("(error: terminal node)")
    }
    if (nodes(newtree, node_id)[[1]]$psplit$ordered) {
        sp <- nodes(newtree, node_id)[[1]]$psplit$splitpoint
        vn <- nodes(newtree, node_id)[[1]]$psplit$variableName
        if (left) {
            op <- "<="
        }
        else {
            op <- ">"
        }
        return(paste(vn, op, sp))
    }
    else {
        psplit <- nodes(newtree, node_id)[[1]]$psplit
        if (left) {
            l <- as.logical(psplit$splitpoint)
        }
        else {
            l <- as.logical(!psplit$splitpoint)
        }
        r <- paste(attr(psplit$splitpoint, "levels")[l], sep = "", 
            collapse = "','")
        return(paste(psplit$variableName, " in ('", r, "')", 
            sep = ""))
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
