library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyFiles)
library(DT)

dashboardPage(
	dashboardHeader(title="AROMA (Agilent RObust Microarray Analysis) Pipeline", titleWidth="25%"),
	dashboardSidebar(disable=FALSE,
		div(id="showPhenoOptions", h4("LOAD PHENOTYPE DATA")),
		hidden(div(id="phenoOptions",
			fluidRow(
				column(12, align="center",
					actionButton("import_pheno_submit", "Launch Import Wizard")
				)
			)
		)),
		hr(),
		div(id="showRawOptions", h4("LOAD RAW DATA")),
		hidden(div(id="rawOptions",
			fluidRow(
				column(12, align="center",
					textOutput("dirText")
				)
			),fluidRow(
				column(1, align="center",
                                        p(" ")
				),column(4, align="left",
					shinyDirButton("dirButton", label="Browse", title="Select Directory", buttonType="default")
				),column(1, align="center",
                                        p(" ")
				),column(5, align="left",
					actionButton("upload_raw_submit", "Upload")
				)
			)
		)),
		hr(),
		div(id="showFiltOptions", h4("PROBE FILTERING")),
		hidden(div(id="filtOptions",
			fluidRow(
				column(12,
					uiOutput("slideFiltDist"),
					uiOutput("slidePerSamples")
				)
			),fluidRow(
				column(12, align="center",
					textOutput("percProbesText")
				),column(12, align="center",
					textOutput("numProbesText")
				)
			)
		)),
		hr(),
		div(id="showNormOptions", h4("NORMALIZATION")),
		hidden(div(id="normOptions",
			fluidRow(
				column(12,
					uiOutput("selNormMethod")
				)
			),fluidRow(
				column(12, align="center",
					actionButton("norm_submit", "Run Normalization")
				)
			)
		)),
		hr(),
		div(id="showBatchOptions", h4("BATCH CORRECTION")),
		hidden(div(id="batchOptions",# style="height:600px;overflow-y:scroll;overflow-x:hidden",
			fluidRow(
				column(12,
					uiOutput("selVarI"),
					uiOutput("selCoVar"),
					uiOutput("selBatch"),
                                        checkboxInput("rmSva", "Remove Unknown Batches", value=TRUE)
				)
			),fluidRow(
				column(12, align="center",
					actionButton("import_ann_submit", "Annotation Import Wizard")
				)
			),fluidRow(
				column(12, align="center",
					actionButton("batch_submit", "Run Batch Correction")
				)
			)
		)),
		hr(),
		div(id="showDiffOptions", h4("DIFFERENTIAL ANALYSIS")),
		hidden(div(id="diffOptions",
			fluidRow(
				column(12,
					uiOutput("selTreatment"),
					uiOutput("selControl"),
					uiOutput("selComps")
				)
			),fluidRow(
				column(12, align="center",
					actionButton("add_comp_submit", "Add Comparison")
				)
			),fluidRow(
				column(12, align="center",
					uiOutput("selPvAdjMethod")
				)
			),fluidRow(
				column(12, align="center",
					hr(),
					actionButton("de_submit", "Run Differential Analysis")
				)
			)
		))
         ),
	dashboardBody(
		useShinyjs(),
		shinyBS::bsModal("importPhenoModal", "Import Phenotype Data", "import_pheno_submit", size="large",
			fluidRow(
				column(4,
					fileInput("fPheno", label="Phenotype File")
				),column(4,
					uiOutput("selSep")
				),column(4,
					textInput("sepT", "Other Seperator", value=":")
				)
			),fluidRow(
				column(6,
					actionButton("load_pheno_submit", "Load")
				)
			),hr(),
			fluidRow(
				column(6,
					uiOutput("selFileNameCol")
				),column(6,
					uiOutput("selDyeCol")
				)
			),hr(),
			fluidRow(
				column(12,
					textOutput("phRowsText"),
					textOutput("phColsText")
				)
			),fluidRow(
				column(12,
					DT::dataTableOutput("phenoDT")
				)
			),fluidRow(
				column(12,
					actionButton("upload_pheno_submit", "Upload")
				)
			)
		),
		shinyBS::bsModal("importAnnModal", "Import Annotation", "import_ann_submit", size="large",
			fluidRow(
				column(12, 
					selectInput("annType", "Select Annotation Type", 
						choices=c(
							"Ensembl Biomart"="mart", 
							"Annotation File"="file"
						), 
						selected="mart"
					) 
				),column(12, 
					hidden(div(id="annFile",
                                                fluidRow(
                                                        column(4,
                                                                fileInput("fAnno", label="Microarray Annotation")
                                                        ),column(4,
                                                                uiOutput("selSepAnno")
                                                        ),column(4,
                                                                textInput("sepTAnno", "Other Seperator", value=":")
                                                        )
                                                ),fluidRow(
                                                        column(6,
                                                                actionButton("load_ann_submit", "Load")
                                                        )
                                                ),hr(),
                                                fluidRow(
                                                        column(12,
                                                                textOutput("annRowsText"),
                                                                textOutput("annColsText")
                                                        )
                                                ),fluidRow(
                                                        column(12,
                                                                uiOutput("selID")
                                                        )
                                                )
					)),div(id="annQuery",
						selectInput("org", "Select Organism", 
							choices=c(
								"Human"="hsapiens_gene_ensembl", 
								"Mouse"="mmusculus_gene_ensembl"
							), 
							selected="hsapiens_gene_ensembl"
						),uiOutput("selProbeID"),
						selectInput("mapID", "Select Mapping ID", 
							choices=c(
								"Ensembl Gene ID"="ensembl_gene_id", 
								"Ensembl Transcript ID"="ensembl_transcript_id",
								"EntrezGene ID"="entrezgene",
								"Gene Symbol"="hgnc_symbol"
							), 
							selected="ensembl_gene_id"
						)
						#uiOutput("mapID")
					)
				)
			),fluidRow(
				column(12,
					DT::dataTableOutput("annDT")
				)
			)
		),
		fluidRow(
			column(6,
				#tags$head(tags$style(HTML(".small-box {height: 70px}"))),
				infoBoxOutput('totalSampleBox'),
				infoBoxOutput('filteredSampleBox'),
				infoBoxOutput('removedSampleBox')
			),column(6,
				valueBoxOutput('varIValueBox'),
				valueBoxOutput('coVarValueBox'),
				valueBoxOutput('batchValueBox')
			)
		),fluidRow(
			tabBox(id="display", title="", width=12,
				tabPanel(value="pdTab", title="Phenotype Data",
					fluidRow(
						column(6,
							p("Select rows from the table to remove unwanted samples."),
							actionButton("filterPh", "Remove Samples")
						)
					),fluidRow(
						column(12,
							DT::dataTableOutput("filtered")
						)
					)
				),tabPanel(value="normTab", title="Normalization",
					fluidRow(column(12,
						tabBox(id="mtv", title="", width=12,
							tabPanel(value="boxPlotTab", title="Box Plot", 
								fluidRow(
									column(6,
										plotOutput('preBoxPlot')
									),column(6,
										plotOutput('postBoxPlot')
									)
								)
							),tabPanel(value="densityPlotTab", title="Density Plot", 
								fluidRow(
									column(6,
										plotOutput('preDensityPlot')
									),column(6,
										plotOutput('postDensityPlot')
									)
								)
							),tabPanel(value="MDPlotTab", title="Mean-Difference Plot", 
								fluidRow(
									column(6,
										plotOutput('preMDPlot')
									),column(6,
										plotOutput('postMDPlot')
									)
								)
							),tabPanel(value="confPlotTab", title="Confounding Plot", 
								fluidRow(
                                                                        column(6,
										sliderInput("confWidth", "Plot Width", min=20, max=100, value=80, step=10)
                                                                        ),column(6,
										sliderInput("confHeight", "Plot Height", min=20, max=100, value=80, step=10)
                                                                        )
								),fluidRow(
									column(8,
										plotOutput('confPlot', height="auto", width="auto")
									)
								)
							),tabPanel(value="princePlotTab", title="Prince Plot",
								fluidRow(
									column(8,
										plotOutput('princePlot', height="auto")
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
				),tabPanel(value="batchTab", title="Batch Correction",
					fluidRow(column(12,
						tabBox(id="bTBox", title="", width=12,
							tabPanel(value="mdsTab", title="MDS Plots", 
                                                                fluidRow(
                                                                        column(6,
                                                                                sliderInput("mdsWidth", "Plot Width", min=20, max=100, value=80, step=10)
                                                                        ),column(6,
                                                                                sliderInput("mdsHeight", "Plot Height", min=20, max=100, value=80, step=10)
                                                                        )
                                                                ),fluidRow(column(12,
                                                                        tabBox(id="mdsBox", title="", width=12,
                                                                                tabPanel(value="preCorTab", title="Before Correction MDS", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('preCorMDS', height="auto", width="auto")
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postCorTab", title="After Known Correction MDS", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postCorMDS')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaCorTab", title="After Unknown Correction MDS", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postSvaCorMDS')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postAggTab", title="After Aggregation MDS", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postAggMDS')
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                ))
							),tabPanel(value="hcTab", title="Hierarchical Clustering", 
                                                                fluidRow(column(12,
                                                                        tabBox(id="hcBox", title="", width=12,
                                                                                tabPanel(value="postHcTab", title="Known Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postHcPlot')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaHcTab", title="Unknown Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postSvaHcPlot')
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                ))
							),tabPanel(value="princeTab", title="Prince Plot", 
                                                                fluidRow(column(12,
                                                                        tabBox(id="princeBox", title="", width=12,
                                                                                tabPanel(value="postPrincePlotTab", title="Known Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postPrincePlot')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaPrincePlotTab", title="Unknown Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postSvaPrincePlot')
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                ))
                                                        )
						)
					))
				),tabPanel(value="diffTab", title="Differential Expression",
					fluidRow(column(12,
						tabBox(id="diffTBox", title="", width=12,
							tabPanel(value="diffTableTab", title="Differential Expression Tables", 
								fluidRow(
									column(12,
										uiOutput('selCompDeTable')
									)
								),fluidRow(
									column(12,
										DT::dataTableOutput("deTable")
									)
								)
							),tabPanel(value="diffVenn", title="P.Value Adjustment Venn", 
								fluidRow(
									column(12,
										plotOutput('venn')
									)
								)
							),tabPanel(value="diffVolcano", title="Volcano Plot", 
								fluidRow(
									column(3,
										sliderInput("lfcThr", "LogFC Threshold", min=0, max=3, value=0.8, step=0.1)
									),column(3,
										sliderInput("pvThr", "P.Value Threshold", min=1, max=10, value=5, step=1)
									),column(3,
										sliderInput("xAxis", "X Axis Bounds", min=1, max=10, value=6, step=1)
									),column(3,
										sliderInput("yAxis", "Y Axis Bounds", min=1, max=10, value=8, step=1)
									)
								),fluidRow(
									column(12,
										checkboxInput("gName", "Display Gene Names", value=FALSE)
									)
								),fluidRow(
									column(12,
										plotOutput('volcanoPlot')
									)
								)
							)
						)
					))
				)
			)
		)
        )
)
