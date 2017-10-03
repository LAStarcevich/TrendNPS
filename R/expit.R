#' @export
#' 
#' @title "Expit" transform, or the backwards logit transform.
#'   
#' @description Apply the standard backwards logit transform to a numeric vector
#'   of values.
#'   
#' @param x	A numeric vector.
#'   
#' @details The \code{expit} function simply applies the standard backwards
#'   logit transform to the provided numeric vector \code{x}: \eqn{1 / (1 +
#'   exp(-x))}.
#'   
#' @author Leigh Ann Starcevich and Jason Mitchell of Western EcoSystems 
#'   Technology, Inc.
#'   
#' @seealso \code{TrendNPS::logit}
#'   
#' @examples 
#' x <- c(-0.8472979,-1.3862944,2.1972246)
#' y <- expit(x)
#' 
#' y
#' # 0.3 0.2 0.9
expit <- function(x) {
  y <- 1 / (1 + exp(-x))
  return(y)
}