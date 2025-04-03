#' .. content for \description{} (no empty lines) ..
#'
#' .. content for \details{} ..
#'
#' @title

#' @return
#' @author dylanmr
#' @export
add_new_annotation <- function(obj, new_name, name, string) {

  obj[[new_name]] <- ifelse(grepl(string, obj[[name]][,1]), "ko", "wt")
  return(obj)

}
