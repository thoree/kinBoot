#'
#' Bootstrap confidence intervals
#'
#' Some alternative confidence intervals are implemented. PR TEST.
#'
#' @param x Double vector. Bootstrap sample.
#' @param method Character. See details.
#' @param conf.level Double i n(0,1). Confidence level.

#'
#' @return A vector of length 2
#' in which the first element is the lower bound and the second element is the upper bound.
#'
#' @details The implemented methods are `perc`, the standard using the `quantile` function R
#' to find the percentile; `bca`, see `coxed::bca`
#' (his method fails if the input is a vector of constant values), and 'student' based on
#' the normal approximation.
#'
#' @export
#'
#' @examples
#' library(coxed)
#' x = rnorm(100000)
#' interval(x, method = "perc")
#' interval(x, method = "bca")
#' interval(x, method = "student")
#' interval(x, method = "Not implemented")
#'
#' # Example Efron & Tibshirani p 154
#' x = c(52, 104, 146, 10, 51, 30, 40, 27, 46)
#' interval(x, method = "bca", conf.level = 0.9)


interval = function(x, method = "perc", conf.level = 0.95){
  low = (1 - conf.level)/2
  high = 1 - low
  if(method == "perc")
    interval = quantile(x, c(low,  high))
  else if(method == "bca")
    interval = coxed::bca(x, conf.level = conf.level)
  else if(method == "student"){
      est = mean(x)
      n = length(x)
      s = sd(x)
      tq = qt(low,df = n-1,lower.tail = FALSE)
      interval = c(est - tq*s, est + tq*s)
    }
  else  stop(method)
  interval
}





