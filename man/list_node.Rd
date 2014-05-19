\name{list_node}
\alias{list_node}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{ list_node
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
list_node(newtree, node_id = 1, start_criteria = character(0))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{newtree}{
%%     ~~Describe \code{newtree} here~~
}
  \item{node_id}{
%%     ~~Describe \code{node_id} here~~
}
  \item{start_criteria}{
%%     ~~Describe \code{start_criteria} here~~
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
function (newtree, node_id = 1, start_criteria = character(0)) 
{
    if (nodes(newtree, node_id)[[1]]$terminal) {
        prediction <- ptree_y(newtree, node_id)
        ypred <- paste(start_criteria, ",y =", prediction, ";")
        return(ypred)
    }
    left_node_id <- ptree_left(newtree, node_id)
    right_node_id <- ptree_right(newtree, node_id)
    if (is.null(left_node_id) != is.null(right_node_id)) {
        print("left node ID != right node id")
    }
    ypred <- character(0)
    if (!is.null(left_node_id)) {
        new_criteria <- paste(start_criteria, ptree_criteria(newtree, 
            node_id, T), sep = ",")
        if (1 == node_id) 
            new_criteria <- ptree_criteria(newtree, node_id, 
                T)
        ypred <- list_node(newtree, left_node_id, new_criteria)
    }
    if (!is.null(right_node_id)) {
        new_criteria <- paste(start_criteria, ptree_criteria(newtree, 
            node_id, F), sep = ",")
        if (1 == node_id) 
            new_criteria <- ptree_criteria(newtree, node_id, 
                F)
        ypred <- paste(ypred, list_node(newtree, right_node_id, 
            new_criteria))
    }
    return(ypred)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
