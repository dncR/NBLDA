# All classes..

###### 1. Class: nblda #####
#' \code{nblda} object
#'
#' @description This object is the main class for NBLDA package. It stores inputs, results and call info for the trained model.
#'
#' @details Objects can be created by calls of the form \code{new("nblda", ...)}. This type
#' of objects is created as a result of \code{trainNBLDA} function of \code{NBLDA} package.
#' It is then used in \code{predict} function for predicting class labels of new samples.
#'
#' @section Slots:
#'
#' \describe{
#'   \item{\code{input}:}{add description here}
#'   \item{\code{result}:}{add description here}
#'   \item{\code{call}:}{add description here}
#' }
#'
#' @note nblda.trained icerisinde yer alan slotlarin aciklamalari burada verilecek.
#'
#' @author Dincer Goksuluk
#'
#' @docType class
#' @name nblda-class
#' @rdname nblda-class
#'
#' @exportClass nblda
setClass("nblda",
         slots = c(input = "list",
                   result = "list",
                   call = "list"))

# setValidity("discrete.train", function(object){
#   TRUE    ##### SLOTLARIN VALIDITY'SI YAPILACAK.
# })


setOldClass(c("confusionMatrix", "train"))
# setClassUnion("MLSeq.train", c("train", "voom.train", "discrete.train"))
# setClassUnion("confMat", c("confusionMatrix"))


###### 2. SubClass: nblda.trained  #######


###### 3. SubClass: nblda.metaData  #######
