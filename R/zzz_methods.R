
# NOTE: If setMethod is included in the same file with the corresponting function 'predict.nblda', then aliases
# are not required to be defined.

#' @rdname predict
#' @aliases predict,nblda-method
#'
#' @export
setMethod(f = "predict", signature = signature(object = "nblda"), definition = predict.nblda)


###### Show Methods  ########

#' @title Show Method for the \code{nblda} Class
#'
#' @description  Show the results of an object returned from \code{\link{trainNBLDA}} on R console.
#'
#' @param object an object of class \code{nblda} to be printed.
#'
#' @author Dincer Goksuluk
#'
#' @examples
#' \dontrun{
#' set.seed(2128)
#' counts <- generateCountData(n = 20, p = 10, K = 2, param = 1, sdsignal = 0.5, DE = 0.8,
#'                             allZero.rm = FALSE, tag.samples = TRUE)
#' x <- t(counts$x + 1)
#' y <- counts$y
#' xte <- t(counts$xte + 1)
#' ctrl <- control(folds = 2, repeats = 2)
#'
#' fit <- trainNBLDA(x = x, y = y, type = "mle", tuneLength = 10,
#'                   metric = "accuracy", train.control = ctrl)
#'
#' show(fit)
#' }
#'
#' @name show
#' @rdname show
#'
#' @importFrom methods show
#' @method show nblda
show.nblda <- function(object){
  cat("DONE...")
}

#' @rdname show
#'
#' @export
setMethod(f = "show", signature = signature(object = "nblda"), definition = show.nblda)


