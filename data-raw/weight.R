## code to prepare `DATASET` dataset goes here
weight <-
  structure(
    list(
      traits = c("sypp", "dtf", "rpp", "ppr", "ppp", "spp", "pw"),
      ew = c(1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000, 1.0000),
      h2 = c(0.6947, 0.3244, 0.6861, 0.7097, 0.8336, 0.2201, 0.5226)
    ),
    class = "data.frame",
    row.names = c(NA,-7L)
  )
usethis::use_data(weight, overwrite = TRUE)
