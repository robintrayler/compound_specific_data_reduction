#' @title parse isodat data
#'
#' @export
#'


parse_isodat <- function(isodat_csv,
                  type = c('EA',
                           'TCEA',
                           'Gasbench')) {
  data <- switch(type,
                 EA = parse_EA_csv(isodat_csv),
                 TCEA = parse_TCEA_csv(isodat_csv),
                 Gasbench = print('Gasebench'))
  return(data)
}
