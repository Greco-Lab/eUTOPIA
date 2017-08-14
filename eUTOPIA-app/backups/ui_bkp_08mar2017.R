suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyFiles))
suppressMessages(library(DT))

jsCode <- "
shinyjs.init = function() {
	var doit;
	function resizePlot() {
		var id=this.id;
		//alert(id);
		sizeableDiv = document.getElementById(id)
		//alert(sizeableDiv.id);
		//var plotImg = this.getElementsByTagName('IMG')[0];
		var plotDiv = sizeableDiv.getElementsByClassName('shiny-plot-output')[0];
		var plotImg = plotDiv.getElementsByTagName('IMG')[0];
		//nWidth=this.clientWidth-5;
		//nHeight=this.clientHeight-5;

		var csSizeable = getComputedStyle(sizeableDiv);
                var csPlot = getComputedStyle(plotDiv);
		var paddingXSizeable = parseFloat(csSizeable.paddingLeft) + parseFloat(csSizeable.paddingRight);
                var paddingXPlot = parseFloat(csPlot.paddingLeft) + parseFloat(csPlot.paddingRight);
		var paddingYSizeable = parseFloat(csSizeable.paddingTop) + parseFloat(csSizeable.paddingBottom);
		var paddingYPlot = parseFloat(csPlot.paddingTop) + parseFloat(csPlot.paddingBottom);

		var borderXSizeable = parseFloat(csSizeable.borderLeftWidth) + parseFloat(csSizeable.borderRightWidth);
		var borderXPlot = parseFloat(csPlot.borderLeftWidth) + parseFloat(csPlot.borderRightWidth);
		var borderYSizeable = parseFloat(csSizeable.borderTopWidth) + parseFloat(csSizeable.borderBottomWidth);
		var borderYPlot = parseFloat(csPlot.borderTopWidth) + parseFloat(csPlot.borderBottomWidth);
		
		//nWidth = sizeableDiv.clientWidth;
		//nHeight = sizeableDiv.clientHeight;
		// Element width and height minus padding and border
		nWidthDiv = sizeableDiv.offsetWidth - paddingXSizeable - borderXSizeable;
		nHeightDiv = sizeableDiv.offsetHeight - paddingYSizeable - borderYSizeable;
		nWidthImg = nWidthDiv - paddingXPlot - borderXPlot;
		nHeightImg = nHeightDiv - paddingYPlot - borderYPlot;
		
		if(nHeightImg==0){
                	nHeightDiv = nWidthDiv;
                	nHeightImg = nWidthImg;
                }
		nWidthDivStr = nWidthDiv + 'px';
                nHeightDivStr = nHeightDiv + 'px';
                nWidthImgStr = nWidthImg + 'px';
                nHeightImgStr = nHeightImg + 'px';
                
                plotDiv.style.width = nWidthDivStr;
                plotDiv.style.height = nHeightDivStr;
                plotImg.style.width = nWidthImgStr;
                plotImg.style.height = nHeightImgStr;

		//clearTimeout(doit);
		//doit = setTimeout(function() {
		//	//alert(this.id);
		//	var id=sizeableDiv.id;
		//	var paramH=id + '_height';
		//	var paramW=id + '_width';
		//	Shiny.onInputChange(paramH, nHeightImg);
		//	Shiny.onInputChange(paramW, nWidthImg);
		//}, 100)
	}

	function replot(){
        	var id=this.id;
                sizeableDiv = document.getElementById(id)

                var plotDiv = sizeableDiv.getElementsByClassName('shiny-plot-output')[0];
                var plotImg = plotDiv.getElementsByTagName('IMG')[0];
                
                var csSizeable = getComputedStyle(sizeableDiv);
                var csPlot = getComputedStyle(plotDiv);
                var paddingXSizeable = parseFloat(csSizeable.paddingLeft) + parseFloat(csSizeable.paddingRight);
                var paddingXPlot = parseFloat(csPlot.paddingLeft) + parseFloat(csPlot.paddingRight);
		var paddingYSizeable = parseFloat(csSizeable.paddingTop) + parseFloat(csSizeable.paddingBottom);
		var paddingYPlot = parseFloat(csPlot.paddingTop) + parseFloat(csPlot.paddingBottom);

		var borderXSizeable = parseFloat(csSizeable.borderLeftWidth) + parseFloat(csSizeable.borderRightWidth);
		var borderXPlot = parseFloat(csPlot.borderLeftWidth) + parseFloat(csPlot.borderRightWidth);
		var borderYSizeable = parseFloat(csSizeable.borderTopWidth) + parseFloat(csSizeable.borderBottomWidth);
		var borderYPlot = parseFloat(csPlot.borderTopWidth) + parseFloat(csPlot.borderBottomWidth);

                // Element width and height minus padding and border
		nWidthDiv = sizeableDiv.offsetWidth - paddingXSizeable - borderXSizeable;
		nHeightDiv = sizeableDiv.offsetHeight - paddingYSizeable - borderYSizeable;
		nWidthImg = nWidthDiv - paddingXPlot - borderXPlot;
		nHeightImg = nHeightDiv - paddingYPlot - borderYPlot;

                var paramH=id + '_height';
                var paramW=id + '_width';
                Shiny.onInputChange(paramH, nHeightImg);
                Shiny.onInputChange(paramW, nWidthImg);
        }

	//var plotDiv = document.getElementsByClassName('shiny-plot-output')[0];
	//addResizeListener(plotDiv, resizePlot);
	//var plotDivs = document.getElementsByClassName('shiny-plot-output');
	var plotDivs = document.getElementsByClassName('sizeable');
	//alert(plotDivs.length)
	var i;
	for (i=0; i<plotDivs.length; i++) {
		addResizeListener(plotDivs[i], resizePlot);
		plotDivs[i].addEventListener('mouseup', replot);
		//plotDivs[i].onresize = function() {resizePlot()};
	}
}"

appCSS <- "
	.main-sidebar { background-color: white !important; }
	.sidebar { color: black; max-height: 900px; overflow-y: scroll;}
	//.panel { background-color: #222d32; }
	.panel-title a { font-weight: bold; color: white !important; }
	.panel-warning .panel-heading { background-color: #00c0ef; }
	.panel-warning { border-color: #8de9ff; }
	.panel-danger .panel-heading { background-color: #dd4b39; }
	.panel-success .panel-heading { background-color: #00a65a; }
	//.info-box-icon { height: 40px; width: 40px; }
	.sizeable { resize: both; overflow: auto; height: auto; width: 60%; padding: 5px 5px 5px 5px; box-shadow: 2px 2px 2px #888888; display: inline-block; }
        .sizeable:active { width: 0; height: 0; }
	//.shiny-plot-output { height: 100% !important; }

        //Add temp jacket div for plots
	.jacket { height: auto; width: 60%; }
        //Or use the plot div itself for CSS effect
	//.shiny-plot-output { width: 60%; }

	//For Loading Screen
	#loading-content {
		position: fixed;
		background: #000000 !important;
		z-index: 10000;
		left: 0;
		right: 0;
		top: 0;
		padding-top: 50px;
		bottom: 0;
		font-size: 50px;
		text-align: center;
		color: #FFFFFF;
	}

	#loading-gif { 
		position: absolute;
		opacity: 0.7; 
		display: block;
		margin-left: auto;
		margin-right: auto;
		z-index: 10000;
		left: 30%;
		top: 5%;
	}

	//For plot height
	//.shiny-plot-output { height:100vh !important; }

	//For plot loading
	//.plot-container { position: relative; }
	//#loading-spinner {
	//	position: absolute;
	//	left: 50%;
	//	top: 50%;
	//	z-index: -1;
	//	margin-top: -33px;  /* half of the spinner's height */
	//	margin-left: -33px; /* half of the spinner's width */
	//}
	//#hcPlot.recalculating { z-index: -2; }
	//#.shiny-plot-output .recalculating { z-index: -2; }
"

dashboardPage(
	dashboardHeader(title="AROMA (Agilent RObust Microarray Analysis) Pipeline", titleWidth="25%"),
	dashboardSidebar(disable=FALSE,
		useShinyjs(),
		extendShinyjs(text=jsCode),
		tags$head(tags$script(src="resizing.js")),
		inlineCSS(appCSS),
                bsCollapse(id="bsSidebar", open="LOAD PHENOTYPE DATA",
                        bsCollapsePanel("LOAD PHENOTYPE DATA", style="warning",
                                fluidRow(
                                        column(12, align="center",
                                                #actionButton("import_pheno_submit", "Phenotype Import Wizard")
                                                shinyBS::bsButton("import_pheno_submit", label="Import Phenotype Data", style="danger", icon=icon("exclamation-circle")),
						shinyBS::bsTooltip("import_pheno_submit", "Launch a graphical window, to configure import of phenotype data from a file!", placement="bottom")
                                        )
                                )
                        ),bsCollapsePanel("LOAD RAW DATA", style="danger",
                                fluidRow(
                                        column(1, align="center",
                                                p("")
                                        ),column(11, align="left",
                                                shinyDirButton("dirButton", label="Browse", title="Select Directory", buttonType="default"),
						shinyBS::bsTooltip("dirButton", "Browse local system directories and select the directory containing the RAW data files!", placement="bottom")
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
                                                uiOutput("selNormMethod"),
                                                uiOutput("selNormMethod2")
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
                                                checkboxInput("rmSva", "Remove Unknown Batches", value=TRUE),
						shinyBS::bsTooltip("rmSva", "Hidden batch variables are identified and removed. These are not present in the user provided phenotype data and are not associated with the 'Variables of Interest' or 'Co-variates'!", placement="top")
                                        )
                                ),fluidRow(
                                        column(12, align="left",
                                                #actionButton("import_ann_submit", "Annotation Import Wizard"),
                                                shinyBS::bsButton("import_ann_submit", label="Import Annotation", style="danger", icon=icon("exclamation-circle")),
						p(" ")
                                        )
                                ),fluidRow(
                                        column(12, align="left",
                                                actionButton("batch_submit", "Run Batch Correction")
                                        )
                                )
                        ),bsCollapsePanel("DIFFERENTIAL ANALYSIS", style="danger",
                                fluidRow(
                                        column(12,
                                                uiOutput("selTreatment")
                                        )
                                ),fluidRow(
                                        column(12,
                                                uiOutput("selControl")
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                actionButton("add_comp_submit", "Add Comparison")
                                        )
                                ),fluidRow(
                                        column(12,
                                                uiOutput("selComps")
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
				column(1,
					actionButton("load_pheno_submit", "Load")
				),column(2, align="left",
					textOutput("phRowsText"),
					textOutput("phColsText")
				)
			),hr(),
			fluidRow(
				column(6,
					uiOutput("selFileNameCol")
				),column(6,
					uiOutput("selDyeCol")
				)
			),fluidRow(
				column(12, align="right",
					#actionButton("upload_pheno_submit", "Import")
					shinyBS::bsButton("upload_pheno_submit", label="Import", style="info", icon=icon("hand-o-right"))
				)
			),hr(),
			fluidRow(
				column(12,
					DT::dataTableOutput("phenoDT")
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
                                                        column(12, align="right",
                                                                #actionButton("load_ann_submit", "Import")
								shinyBS::bsButton("load_ann_submit", label="Import", style="info", icon=icon("hand-o-right"))
                                                        )
                                                ),hr(),
                                                fluidRow(
                                                        column(12,
                                                                textOutput("annRowsText"),
                                                                textOutput("annColsText")
                                                        )
                                                ),fluidRow(
                                                        column(12,
                                                                #p(id="annMaps", "Mappability: 'Waiting for upload!'")
                                                                uiOutput("htmlAnnMaps")
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
		#shinyBS::bsModal("resizeImageModal", "Resize Image", "", size="small",
		#	fluidRow(
		#		column(12, 
		#		)
		#	)
		#),
		div(id="loading-content",
			img(id="loading-gif", src="screen-loading.gif")
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
									column(12,
										#div(id="confPlotDiv", class="sizeable", 
										div(id="confPlotDiv", class="jacket", 
											plotOutput('confPlot', height="auto", width="auto")
										)
									)
								)
							),tabPanel(value="princePlotTab", title="Prince Plot",
								fluidRow(
									column(12,
										#div(id="princePlotDiv", class="sizeable", 
										div(id="princePlotDiv", class="jacket", 
											plotOutput('princePlot', height="auto", width="auto")
										)
									)
								)
							),tabPanel(value="hcTab", title="Hierarchical Clustering",
								fluidRow(
									column(12,
										#div(id="hcPlotDiv", class="plot-container sizeable",
                                                                                #        tags$img(src="spinner.gif", id="loading-spinner"),
										#div(id="hcPlotDiv", class="sizeable",
										div(id="hcPlotDiv", class="jacket", 
                                                                                        plotOutput('hcPlot', height="auto", width="auto")
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
                                                                fluidRow(column(12,
                                                                        tabBox(id="mdsBox", title="", width=12,
                                                                                tabPanel(value="preCorTab", title="Before Correction", 
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="preCorMDSDiv", class="sizeable", 
													div(id="preCorMDSDiv", class="jacket", 
														plotOutput('preCorMDS', height="auto", width="auto")
													)	
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postCorTab", title="Known Corrected", 
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postCorMDSDiv", class="sizeable", 
													div(id="postCorMDSDiv", class="jacket", 
														plotOutput('postCorMDS', height="auto", width="auto")
													)	
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaCorTab", title="Unknown Corrected", 
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postSvaCorMDSDiv", class="sizeable", 
													div(id="postSvaCorMDSDiv", class="jacket", 
														plotOutput('postSvaCorMDS', height="auto", width="auto")
													)	
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postAggTab", title="After Aggregation", 
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postAggMDSDiv", class="sizeable", 
													div(id="postAggMDSDiv", class="jacket", 
														plotOutput('postAggMDS', height="auto", width="auto")
													)	
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
													#div(id="postHcPlotDiv", class="sizeable", 
													div(id="postHcPlotDiv", class="jacket", 
														plotOutput('postHcPlot', height="auto", width="auto")
													)	
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaHcTab", title="Unknown Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postSvaHcPlotDiv", class="sizeable", 
													div(id="postSvaHcPlotDiv", class="jacket", 
														plotOutput('postSvaHcPlot', height="auto", width="auto")
													)	
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
													#div(id="postPrincePlotDiv", class="sizeable", 
													div(id="postPrincePlotDiv", class="jacket", 
                                                                                                        plotOutput('postPrincePlot', height="auto", width="auto")
													)	
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postSvaPrincePlotTab", title="Unknown Corrected",
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postSvaPrincePlotDiv", class="sizeable", 
													div(id="postSvaPrincePlotDiv", class="jacket", 
                                                                                                        plotOutput('postSvaPrincePlot', height="auto", width="auto")
													)	
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
                                                ),column(2,
                                                        downloadButton("exportDE", "Export Differential Tables")
                                                ),column(2,
                                                        downloadButton("exportRpt", "Generate Report")
                                                )
					),fluidRow(column(3,
                                                        uiOutput("slideLfcThr")
                                                #),column(3,
                                                #        uiOutput("slidePvThr")
                                                ),column(3,
                                                        uiOutput("slideAdjPvThr")
                                                #),column(3,
                                                #        uiOutput("slideScoreThr")
                                                #)
					#),fluidRow(column(3,
					        ),column(2,
                                                        fluidRow(column(12,
                                                                #htmlOutput("textLfcThr"),
                                                                #htmlOutput("textPvThr"),
                                                                #htmlOutput("textAdjPvThr"),
                                                                #htmlOutput("textScoreThr")
                                                                htmlOutput("htmlValTable")
                                                        ))
					        ),column(4,
                                                        fluidRow(column(12,
                                                                htmlOutput("htmlDegInfo")
                                                        ))
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
							),tabPanel(value="diffVenn", title="Differential Gene Sets Intersection", 
								fluidRow(
									column(4,
                                                                               uiOutput("selIntersectComps")
									),column(8,
                                                                                #div(id="intersectPlotDiv", class="sizeable", 
                                                                                div(id="intersectPlotDiv", class="jacket", 
                                                                                        plotOutput('intersectPlot', height="auto", width="auto")
                                                                                )
									)
								)
							),tabPanel(value="diffVolcano", title="Volcano Plot", 
								fluidRow(
									column(3,
										sliderInput("xAxis", "X Axis Bounds", min=1, max=7, value=5, step=1)
									),column(3,
										sliderInput("yAxis", "Y Axis Bounds", min=1, max=30, value=20, step=1)
									)
								),fluidRow(
									column(12,
										checkboxInput("gName", "Display Gene Names", value=FALSE)
									)
								),fluidRow(
									column(12,
										plotOutput('volcanoPlot', height="auto", width="auto")
									)
								)
							)
                                                )
					))
				),tabPanel(value="visExp", title="Visualize Expression",
					fluidRow(column(12,
						tabBox(id="visExpTBox", title="", width=12,
                                                        tabPanel(value="expBoxplot", title="Boxplot", 
                                                                fluidRow(
                                                                        column(4,
                                                                                uiOutput("selExpGenes"),
                                                                                uiOutput("selConditions")
                                                                        ),column(8,
                                                                                plotOutput("expressionBoxPlot", height="auto", width="auto")
                                                                        )
                                                                )
                                                        ),tabPanel(value="expHeatmap", title="Heatmap", 
                                                                fluidRow(
                                                                        column(4,
                                                                                uiOutput("slidePercDE"),
                                                                                uiOutput("selHeatComps"),
                                                                                uiOutput("selConditionsHeat")
                                                                        ),column(8,
                                                                                plotOutput("expressionHeatmap", height="auto", width="auto")
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
