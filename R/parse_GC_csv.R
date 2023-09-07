parse_GC_csv <- function(data) {
  data <- data %>%
    select(id1 = `Identifier 1`,
           sample = `Identifier 2`,
           amount = Comment,
           seq_nr = Row,
           retention_time = Rt,
           end = End,
           start = Start,
           amp_44 = `Ampl. 44`,
           d13C_measured = `d 13C/12C`) %>%
    return()
}