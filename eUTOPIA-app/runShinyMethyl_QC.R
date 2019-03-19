args=(commandArgs(TRUE))

print("args:")
print(args)

#tempRdata <- eval(parse(text=args[1]))
tempRdata <- args[1]
print("tempRdata:")
print(tempRdata)
load(tempRdata)

library(shiny)
library(shinyMethyl)

#shinyMethyl::runShinyMethyl(RGsetSummary)
print("Running App...")
appList <- list(ui=uiShinyMethyl, server=serverShinyMethyl)
runApp(appList, launch.browser=TRUE)
print("Shiny Instance Started!")
