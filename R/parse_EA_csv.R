#' @title parse EA data
#'
#' @export
#'

parse_GC_csv <- function(data) {
  data <- data %>%
    select(identifier_1 = `Identifier 1`,
           sample = `Identifier 2`,
           amount = Comment,
           seq_nr = Row,
           retention_time = Rt,
           end = End,
           start = Start,
           area_44 = `Area 44`,
           d13C_measured = `d 13C/12C`) %>%
    return()
}
