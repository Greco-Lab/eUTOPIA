library(shiny)
library(shinyMethyl)

#Load temp shinyMethyl Rdata
tempRdata <- Sys.getenv("SHINY_METHYL_TEMP")
message("tempRdata:")
message(tempRdata)
tempRdata <- "../eUTOPIA/RGsetSummary.RData"
load(tempRdata)

# Create Shiny app ----
message("Running App...")

#appList <- list(ui=uiShinyMethyl, server=serverShinyMethyl)
#runApp(appList)

message("Shiny Instance Started!")
shinyApp(ui=uiShinyMethyl, server=serverShinyMethyl)
