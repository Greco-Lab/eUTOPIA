args=(commandArgs(TRUE))

message("args:")
message(args)

#tempRdata <- eval(parse(text=args[1]))
tempRdata <- args[1]
message("tempRdata:")
message(tempRdata)
load(tempRdata)

library(shiny)
library(shinyMethyl)

#shinyMethyl::runShinyMethyl(RGsetSummary)
message("Running App...")
appList <- list(ui=uiShinyMethyl, server=serverShinyMethyl)
runApp(appList, launch.browser=TRUE, host="0.0.0.0")
message("Shiny Instance Started!")
