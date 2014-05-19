\name{plot.minmax}
\alias{plot.minmax}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plot.minmax
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
plot.minmax(My, X, Y)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{My}{
%%     ~~Describe \code{My} here~~
}
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{Y}{
%%     ~~Describe \code{Y} here~~
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
function (My, X, Y) 
{
    mymat <- My$M
    my.y <- My$y
    my.y.val <- as.numeric(strsplit(trim(my.y), " ")[[1]][3])
    my.y.pct <- ecdf(Y)(my.y.val)
    var.nms <- rownames(mymat)
    act.vars <- apply(mymat, 1, function(x) {
        !all(abs(x) == Inf)
    })
    max.y <- sum(act.vars) + 1
    rbw <- rainbow(n = nrow(mymat))
    node.index <- 1:length(Y)
    for (i in 1:nrow(mymat)) {
        node.index <- intersect(node.index, which(X[, i] > mymat[i, 
            1] & X[, i] <= mymat[i, 2]))
    }
    H <- hist(ecdf(Y)(Y[node.index]), breaks = seq(0, 1, by = 0.1), 
        plot = FALSE)
    scale.factor <- max.y/max(H$density)
    plot(NA, xlim = c(0, 1), ylim = c(0, max.y), ylab = "", xlab = "Percentile", 
        main = paste("Mean =", signif(my.y.val, 2), ", n =", 
            length(node.index)), bty = "n", yaxt = "n")
    barplot(scale.factor * H$density, width = 0.1, col = rgb(0, 
        0, 0, 0.05), border = rgb(0, 0, 0, 0.1), add = TRUE, 
        space = 0, yaxt = "n")
    mu.Y <- mean(Y[node.index])
    print(mean(Y[node.index]))
    segments(ecdf(Y)(mu.Y), 0, ecdf(Y)(mu.Y), max.y, col = rgb(0, 
        0, 0, 0.5), lwd = 2)
    j <- 1
    for (i in which(act.vars)) {
        F.x <- ecdf(X[, var.nms[i]])
        lo <- ifelse(mymat[i, 1] == -Inf, 0, F.x(mymat[i, 1]))
        hi <- ifelse(mymat[i, 2] == Inf, 1, F.x(mymat[i, 2]))
        polygon(c(lo, lo, hi, hi), c(j - 0.5, j + 0.5, j + 0.5, 
            j - 0.5), col = makeTransparent(rbw[i], 0.5), border = NA)
        text(mean(c(lo, hi)), j, rownames(mymat)[i])
        j <- j + 1
    }
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
