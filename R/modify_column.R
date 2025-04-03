#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
modify_column <- function(obj, column_to_update, original, update) {
  obj[[column_to_update]] <- ifelse(obj[[column_to_update]][,1] %in% original, update, obj[[column_to_update]][,1])
  return(obj)
}