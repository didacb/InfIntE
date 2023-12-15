#' Compute the I stastitc
#' Obtain the I statistic from the compression values result of abduction
#' @param abduced data frame where each row contains a pair of OTU, the effect on the abundance and the compression supporting the effect. Result of abduce function.
#'
#' @return data frame where each row contains a pair of OTU, the effect on the abundance and the I statistic supporting the effect.
#' @export
#'
#' @examples
#' get_I_values(abduced)
get_I_values <- function(abduced) {

  #usethis::use_data_table()
  .datatable.aware = TRUE
  #Save col names
  nms<- colnames(abduced)
  
  #Transform to dt
  abduced<- setDT(abduced)
  
  #Sum same interactions
  abduced<- abduced[,.(sum(comp)), by=.(sp1,sp2,lnk)]
  colnames(abduced)<- nms
  
  #one minus opposite
  abduced<- abduced[,.(lnk[comp==max(comp)][1], minus(comp)), by=.(sp1,sp2)]
  
  #back to df
  abduced<- as.data.frame(abduced)
  colnames(abduced)<- nms
  
  return(abduced)
}
