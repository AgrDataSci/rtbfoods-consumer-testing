library("ClimMobTools")
library("readr")

key <- ""

dat <- getDataCM(key, "SPGSP22",
                 userowner = "nextgennacrritricot",
                 as.data.frame = TRUE,
                 pivot.wider = TRUE,
                 tidynames = TRUE)

write_excel_csv(dat, file = "data/tricot-consumer-sweetpotato-uganda.csv")
