library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyFiles)
library(DT)

dashboardPage(
	dashboardHeader(title="AROMA (Agilent RObust Microarray Analysis) Pipeline", titleWidth="25%"),
	dashboardSidebar(disable=FALSE,
                sidebarMenu(
                            menuItem("Upload", tabName="upload",
				h5("Select Directory"),
				textOutput("dirText"),
				shinyDirButton("dirButton", label="Browse...", title="Select Directory", buttonType="default"),
				fileInput("fPheno", label="Phenotype Table File"),
				actionButton("upload_submit", "Load Phenotype Data")
			    ),
                            menuItem("Normalization", tabName="norm"),
                            menuItem("Batch Correction", tabName="batch"),
                            menuItem("Differential Expression", tabName="diff"),
                            menuItem("Results", tabName="res")
                )
         ),
	dashboardBody(
                #titlePanel("AROMA (Agilent RObust Microarray Analysis) Pipeline"),
                tabItems(
                        tabItem(tabName="upload",
                                h3("Upload data to the pipeline"),
                                fluidRow(
                                        column(6,
                                                wellPanel(
                                                        fluidRow(
                                                                column(6,
                                                                        h5("Select Directory")#,
                                                                        #textOutput("dirText"),
                                                                        #shinyDirButton("dirButton", label="Browse...", title="Select Directory")
                                                                ),column(6,
                                                                        h5("Select Directory")#,
                                                                        #fileInput("fPheno", label="Phenotype Table File")
                                                                #),column(4,
								#	uiOutput("selVarI")
                                                                )
                                                        )
                                                )
                                        )
                                ),fluidRow(
                                        column(12,
						h5("Select Directory")
                                                #actionButton("upload_submit", "Load Phenotype Data")
                                        )
                                ),fluidRow(
                                        column(6,
                                                hr(),
                                                h3("Display Area")
                                        )
                                ),fluidRow(
                                        tabBox(id="uDisplay", title="", width=12,
                                                tabPanel(value="pdTab", title="Phenotype Data",
                                                        fluidRow(
                                                                column(9,
                                                                        valueBoxOutput('totalSampleBox'),
                                                                        valueBoxOutput('filteredSampleBox'),
                                                                        valueBoxOutput('removedSampleBox')
                                                                ),column(3,
                                                                        p("Select rows from the table to remove unwanted samples."),
                                                                        actionButton("filterPh", "Remove Samples")
                                                                )
                                                        ),fluidRow(
                                                                column(12,
                                                                        DT::dataTableOutput("filtered")
                                                                )
                                                        )
                                                #), 
                                                #tabPanel("Monitor Technical Variations",
                                                #        #tabsetPanel(
                                                #        tabBox(id="mtv", title="", width=12,
                                                #                tabPanel(value="confPlotTab", title="Confounding Plot", fluidRow(column(8,plotOutput('confPlot')), column(4,helpText("Confounding Plot <HELP TEXT>")))),
                                                #                tabPanel(value="princePlotTab", title="Prince Plot",fluidRow(column(8,plotOutput('princePlot')), column(4,helpText("Prince Plot <HELP TEXT>")))),
                                                #                tabPanel(value="hcTab", title="Hierarchical Clustering",fluidRow(column(8,plotOutput('hcPlot')), column(4,helpText("Hierarchical Clustering Plot <HELP TEXT>"))))
                                                #        )
                                                )
                                        )
                                )
                        ),
                        tabItem(tabName="norm",
                                h3("Normalize and Filter"),
                                fluidRow(
                                	column(6,
                                                wellPanel(
                                                        fluidRow(
                                                                column(6,
									uiOutput("selQDist")
                                                                ),column(6,
									uiOutput("selPerProbes")
                                                                )
                                                        )
                                                )
                                        )
                                ),fluidRow(
                                	column(12,
                                                actionButton("norm_submit", "Run Normalization")
                                        )
                                ),fluidRow(
                                        column(6,
                                                hr(),
                                                h3("Display Area")
                                        )
                                ),fluidRow(
                                        tabBox(id="nDisplay", title="", width=12,
                                                tabPanel(value="mtvTab", title="Monitor Technical Variation",
							fluidRow(column(12,
								tabBox(id="mtv", title="", width=12,
									tabPanel(value="confPlotTab", title="Confounding Plot", 
										fluidRow(
											column(8,
												plotOutput('confPlot')
											)
										)
									),tabPanel(value="princePlotTab", title="Prince Plot",
										fluidRow(
											column(8,
												plotOutput('princePlot')
											)
										)
									),tabPanel(value="hcTab", title="Hierarchical Clustering",
										fluidRow(
											column(8,
												plotOutput('hcPlot')
											)
										)
									)
								)
							))
                                                )
                                        )
                                )
                        ),
                        tabItem(tabName="batch",
				h3("Perform Batch Correction"),
                                fluidRow(
					column(8,
                                                wellPanel(
                                                        fluidRow(
                                                                column(3,
									uiOutput("selVarI")
                                                                ),column(3,
									uiOutput("selCoVar")
                                                                ),column(3,
									uiOutput("selBatch")
                                                                ),column(3,
                                                                        fileInput("fAnno", label="Microarray Annotation"),
									uiOutput("selID")
                                                                )
                                                        )
                                                )
                                        )
                                ),fluidRow(
                                	column(12,
                                                actionButton("batch_submit", "Run Batch Correction")
                                        )
                                ),fluidRow(
                                        column(6,
                                                hr(),
                                                h3("Display Area")
                                        )
                                ),fluidRow(
                                        tabBox(id="bDisplay", title="", width=12,
                                                tabPanel(value="batchTab", title="Corrective Effect",
							fluidRow(column(12,
								tabBox(id="bTBox", title="", width=12,
									tabPanel(value="preCorTab", title="Before Correction", 
										fluidRow(
											column(12,
												plotOutput('preCorMDS')
											)
										)
									),tabPanel(value="postCorTab", title="After Correction", 
										fluidRow(
											column(12,
												plotOutput('postCorMDS')
											)
										)
									),tabPanel(value="postAggTab", title="After Aggregation", 
										fluidRow(
											column(12,
												plotOutput('postAggMDS')
											)
										)
									),tabPanel(value="postHcTab", title="Hierarchical Clustering",
										fluidRow(
											column(12,
												plotOutput('postHcPlot')
											)
										)
									),tabPanel(value="postPrincePlotTab", title="Prince Plot",
										fluidRow(
											column(12,
												plotOutput('postPrincePlot')
											)
										)
									)
								)
							))
                                                )
                                        )
                                )
                        ),
                        tabItem(tabName="diff",
				h3("Differential Expression Analysis"),
                                fluidRow(
					column(6,
                                                wellPanel(
                                                        fluidRow(
                                                                column(4,
									uiOutput("selTreatment")
                                                                ),column(4,
									uiOutput("selControl")
                                                                ),column(4,
									uiOutput("selComps")
                                                                )
                                                        ),fluidRow(
								column(4,
									actionButton("add_comp_submit", "Add Comparison")
								)
                                                        )
                                                )
                                        )
                                ),fluidRow(
                                	column(12,
                                                actionButton("de_submit", "Run Differential Expression Analysis")
                                        )
                                ),fluidRow(
                                        column(6,
                                                hr(),
                                                h3("Display Area")
                                        )
				)
                        ),
                        tabItem(tabName="res",
                                p("Pathway Analysis")
                        )
                )
        )
)
