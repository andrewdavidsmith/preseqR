#' Shakespeare's word count frequencies
#'
#' The Shakespeare's word type frequencies data was from Efron, B., &
#' Thisted, R. (1976).
#'
#' @docType data
#'
#' @usage data(Shakespeare)
#'
#' @details
#' A two-column matrix.  The first column is the frequency j =
#' 1,2,...; and the second column is n_j, the number of unique words
#' appeared j times in Shakespeare's work.
#'
#' @format An object of class \code{data.frame}
#' \describe{
#'  \item{j:}{The count of appearances}
#'  \item{n_j:}{The number of words appearing this many times}
#' }
#'
#' @references
#' Efron, B., & Thisted, R. (1976). Estimating the number of unseen
#' species: How many words did Shakespeare know? Biometrika, 63(3),
#' 435-447.
#'
#' @keywords Shakespeare
#'
#' @examples
#'
#' # load library
#' library(preseqR)
#'
#' # load this data
#' data(Shakespeare)
#'
#' # check out the data
#' head(Shakespeare)
"Shakespeare"
