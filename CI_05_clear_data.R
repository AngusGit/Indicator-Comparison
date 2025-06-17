source("../R/functions.R")


unlink("../data/modelled/*.*", recursive = TRUE)
unlink("../data/summarised/*.*", recursive = TRUE)

unlink("../output/figures/*.*", recursive = TRUE)

unlink("../docs/Reports/*.html")
unlink("../data/model_stage.RData")
unlink("../data/analysis_stage.RData")

## cli_alert_success("Workspace cleared in preparation")          
