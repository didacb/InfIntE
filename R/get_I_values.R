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
  # Obtain final value
  abduced <- plyr::ddply(.data = abduced, .variables = (sp1, sp2, lnk), summarise,
                         comp = max(comp))
  abduced <- plyr::ddply(.data = abduced, .variables = (sp1, sp2), summarise,
    lnk = lnk[comp == max(comp)][1], comp = if (length(comp) > 1) {
      max(comp) - min(comp)
    } else {
      comp
    }
  )
  return(abduced)
}
