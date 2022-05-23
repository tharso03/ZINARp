#' EXPLORATORY DATA ANALYSIS FOR ZINAR(p) PROCESSES
#'
#' This function generates a graph for exploring ZINAR(p) processes.
#'
#'@param x A vector containing a discrete non-negative time series data set.
#'
#'@return Plot time series graph, relative frequency bar plot, autocorrelation function graph and partial autocorrelation function graph on a common plot.
#'
#'@export

explore_zinarp <- function(x){

  if(sum(x%%1!=0)>0 | sum(x<0)>0){
    stop('x must be a ZINAR(p) process.')
  }

    df <- data.frame(
      values = x
    )

  oldpar <- graphics::par(no.readonly = TRUE)
  on.exit(graphics::par(oldpar))

  graphics::par(mfrow = c(2,2))


  stats::ts.plot(df$values, main = 'Time Series Plot', gpars = list(xlab = 'Time', ylab = 'Value'))
  graphics::barplot(table(x)/length(x),
                    xlab = 'Values',
                    ylab = 'Relative Frequency',
                    ylim = c(0,max(table(x)/length(x)) + 0.1), main = 'Relative Frequency')
  stats::acf(x, main = 'ACF')
  stats::pacf(x, main = 'PACF')

}






