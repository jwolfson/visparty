\name{makeTransparent}
\alias{makeTransparent}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
makeTransparent
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
makeTransparent(..., alpha = 0.5)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{\dots}{
%%     ~~Describe \code{\dots} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
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
function (..., alpha = 0.5) 
{
    if (alpha < 0 | alpha > 1) 
        stop("alpha must be between 0 and 1")
    alpha = floor(255 * alpha)
    newColor = col2rgb(col = unlist(list(...)), alpha = FALSE)
    .makeTransparent = function(col, alpha) {
        rgb(red = col[1], green = col[2], blue = col[3], alpha = alpha, 
            maxColorValue = 255)
    }
    newColor = apply(newColor, 2, .makeTransparent, alpha = alpha)
    return(newColor)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
