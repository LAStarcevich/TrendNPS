#' @export
#' 
#' @title Logit transform.
#'   
#' @description Apply the standard logit transform to a numeric vector of
#'   values.
#'   
#' @param x	A numeric vector.
#'   
#' @details The \code{logit} function simply applies the standard logit
#'   transform to the provided numeric vector \code{x}: \eqn{log(x / (1 - x))}.
#'   
#' @author Leigh Ann Starcevich and Jason Mitchell of Western EcoSystems 
#'   Technology, Inc.
#'   
#' @seealso \code{TrendNPS::expit}
#'   
#' @examples 
#' x <- c(0.3,0.2,0.9)
#' y <- logit(x)
#' 
#' y
#' # -0.8472979 -1.3862944  2.1972246
logit <- function(x){
  y <- log(x / (1 - x))
  return(y)
}