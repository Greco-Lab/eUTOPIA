library(shiny)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyFiles)
library(DT)

jsCode <- "
shinyjs.init = function() {
	function resizePlot() {
		var plotImg = this.getElementsByTagName('IMG')[0];
		nWidth=this.clientWidth-5;
		nHeight=this.clientHeight-5;
		nWidthStr=nWidth + 'px';
		nHeightStr=nHeight + 'px';
		plotImg.style.width=nWidthStr;
		plotImg.style.height=nHeightStr;
		//alert('here')
	}
	//var plotDiv = document.getElementsByClassName('shiny-plot-output')[0];
	//addResizeListener(plotDiv, resizePlot);
	//var plotDivs = document.getElementsByClassName('shiny-plot-output');
	var plotDivs = document.getElementsByClassName('sizeable');
	alert(plotDivs.length)
	var i;
	for (i=0; i<plotDivs.length; i++) {
	    addResizeListener(plotDivs[i], resizePlot);
	    //plotDivs[i].onresize = function() {resizePlot()};
	}
}"

dashboardPage(
	dashboardHeader(title="AROMA (Agilent RObust Microarray Analysis) Pipeline", titleWidth="25%"),
	dashboardSidebar(disable=FALSE,
		useShinyjs(),
		extendShinyjs(text=jsCode),
		tags$head(tags$script(src="resizing.js")),
		#tags$script(src="resizeCalling.js"),
                inlineCSS("
                        .main-sidebar { background-color: white !important; }
                        .sidebar { color: black; max-height: 900px; overflow-y: scroll;}
                        //.panel { background-color: #222d32; }
                        .panel-title a { font-weight: bold; color: white !important; }
                        .panel-warning .panel-heading { background-color: #00c0ef; }
                        .panel-warning { border-color: #8de9ff; }
                        .panel-danger .panel-heading { background-color: #dd4b39; }
                        .panel-success .panel-heading { background-color: #00a65a; }
                        //.info-box-icon { height: 40px; width: 40px; }
                        .sizeable { resize: both; overflow: auto; }
                        //.shiny-plot-output { resize: both; overflow: auto; }
                        .plot-container { position: relative; }
                        #loading-spinner {
                                position: absolute;
                                left: 50%;
                                top: 50%;
                                z-index: -1;
                                margin-top: -33px;  /* half of the spinner's height */
                                margin-left: -33px; /* half of the spinner's width */
                        }
                        #hcPlot.recalculating { z-index: -2; }
                        .shiny-plot-output .recalculating { z-index: -2; }

                "),
                bsCollapse(id="bsSidebar", open="LOAD PHENOTYPE DATA",
                        bsCollapsePanel("LOAD PHENOTYPE DATA", style="warning",
                                fluidRow(
                                        column(12, align="center",
                                                actionButton("import_pheno_submit", "Launch Import Wizard")
                                        )
                                )
                        ),bsCollapsePanel("LOAD RAW DATA", style="danger",
                                fluidRow(
                                        column(1, align="center",
                                                p("")
                                        ),column(11, align="left",
                                                shinyDirButton("dirButton", label="Browse", title="Select Directory", buttonType="default")
                                        )
                                ),fluidRow(
                                        column(1, align="center",
                                                p("")
                                        ),column(11, align="left",
                                                textOutput("dirText")
                                        )
                                ),fluidRow(
                                        column(1, align="center",
                                                p("")
                                        ),column(11, align="left",
                                                actionButton("upload_raw_submit", "Upload")
                                        )
                                )
                        ),bsCollapsePanel("PROBE FILTERING", style="danger",
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
                        ),bsCollapsePanel("NORMALIZATION", style="danger",
                                fluidRow(
                                        column(12,
                                                uiOutput("selNormMethod")
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                actionButton("norm_submit", "Run Normalization")
                                        )
                                )
                        ),bsCollapsePanel("BATCH CORRECTION", style="danger",
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
                        ),bsCollapsePanel("DIFFERENTIAL ANALYSIS", style="danger",
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
                        )
                )
         ),
	dashboardBody(
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
                                                                p(id="annMaps", "Mappability: 'Waiting for upload!'")
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
										div(id="test", class="sizeable", plotOutput('preBoxPlot'))
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
										div(id="test", class="sizeable", plotOutput('princePlot', height="auto"))
									)
								)
							),tabPanel(value="hcTab", title="Hierarchical Clustering",
								fluidRow(
									column(8,
										div(id="hcPlotDiv", class="plot-container",
                                                                                        tags$img(src="spinner.gif", id="loading-spinner"),
                                                                                        plotOutput('hcPlot')
                                                                                )
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
                                                                                tabPanel(value="preCorTab", title="Before Correction", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('preCorMDS', height="auto", width="auto")
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postCorTab", title="Known Corrected", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postCorMDS')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaCorTab", title="Unknown Corrected", 
                                                                                        fluidRow(
                                                                                                column(12,
                                                                                                        plotOutput('postSvaCorMDS')
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postAggTab", title="After Aggregation", 
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
					fluidRow(column(4,
                                                        uiOutput('selCompDeTable')
                                                )
					),fluidRow(column(2,
                                                        downloadButton("exportDE", "Export Differential Tables")
                                                ),column(2,
                                                        downloadButton("exportRpt", "Generate Report")
                                                )
                                        ),
					fluidRow(column(12,
						tabBox(id="diffTBox", title="", width=12,
							tabPanel(value="diffTableTab", title="Differential Expression Tables", 
								fluidRow(
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
