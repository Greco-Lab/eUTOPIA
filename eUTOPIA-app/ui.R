suppressMessages(library(shiny))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyFiles))
suppressMessages(library(DT))
suppressMessages(library(rhandsontable))
suppressMessages(library(shinycssloaders))

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

        //Function to hide unwanted bsCollapsePanel
        Shiny.addCustomMessageHandler('hideBSCollapsePanelMessage',
		function(message){
			//alert('The value is ' + message.name);
			//var html = $(\"div[value='\" + message.name + \"']\").html();
			//var obj = $(\"div[value='\" + message.name + \"']\");
			$(\"div[value='\" + message.name + \"']\").addClass(\"hidden\");
			//alert('The obj HTML is ' + html);
			//alert('The obj is ' + obj);
		}
	);

}"

appCSS <- "
	//.main-header { z-index: 100000; }
	.main-sidebar { background-color: white !important; width: 15%; }
	.sidebar { color: black; max-height: 900px; overflow-y: scroll; }
	.sidebar a { color: black !important; }
	.content-wrapper { margin-left: 15%; }
	//.panel { background-color: #222d32; }
	.panel-title a { font-weight: bold; color: white !important; }
	.panel-warning .panel-heading { background-color: #00c0ef; }
	.panel-warning { border-color: #8de9ff; }
	.panel-danger .panel-heading { background-color: #dd4b39; }
	.panel-success .panel-heading { background-color: #00a65a; }
	.panel-info .panel-heading { background-color: #7e46ff; }
	//.info-box-icon { height: 40px; width: 40px; }
	.sizeable { resize: both; overflow: auto; height: auto; width: 60%; padding: 5px 5px 5px 5px; box-shadow: 2px 2px 2px #888888; display: inline-block; }
        .sizeable:active { width: 0; height: 0; }
	//.shiny-plot-output { height: 100% !important; }

        //Add temp jacket div for plots
	.jacket { height: auto; width: 60%; }
        //Or use the plot div itself for CSS effect
	//.shiny-plot-output { width: 60%; }

	.shiny-input-container { padding-left: 1px !important; }

	.contentDiv {
		background: #ECF0F5;
		padding: 5px;
		border: solid;
		border-width: 1px;
		border-radius: 5px;
	}

        #buttonsDiv {
                overflow: hidden;
                width: auto;
        }

	#loading-content {
		//position: fixed;
		position: absolute;
		//background: #663399 !important;
		background: white !important;
		//background: #cfd7d7 !important;
		opacity: 0.8;
		z-index: 1000000;
		left: 0;
		right: 0;
		top: 0;
		//padding-top: 50%;
		bottom: 0;
		font-size: 50px;
		text-align: center;
		color: #black;
	}

	#loading-gif {
		//position: absolute;
		opacity: 0.8;
		display: block;
		margin-left: auto;
		margin-right: auto;
		vertical-align: middle;
		z-index: 1000000;
		//left: 30%;
		//top: 20%;
	}

	.container-fluid { padding-left: 0px; padding-right: 0px; }

	#sel-array-type {
		position: absolute;
		background: #000000;
		opacity: 1;
		z-index: 10000;
		left: 0;
		right: 0;
		height: 100%;
		text-align: center;
		color: #FFFFFF;
	}

	#widget-wrap {
		position: absolute;
		top: 40%;
		left: 40%;
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

        #htmlAnnMaps { color: green; font-weight: bold; }
"

fluidPage(
	useShinyjs(),
	extendShinyjs(text=jsCode),
	tags$head(tags$script(src="resizing.js")),
	inlineCSS(appCSS),
	div(id="loading-content",
		img(id="loading-gif", src="screen-loading.gif"),
                #p("LOADING...")
                #h1(textOutput("loading_text")),
                p(id="loadingText", "LOADING"),
                p("...")
	),
	div(id="sel-array-type",
		div(id="widget-wrap",
			fluidRow(column(12, align="center",
				selectInput("arrType", "Select Platform",
					choices=c(
						"Agilent Expression (2 color)"="ag_exp2",
						"Agilent Expression (1 color)"="ag_exp1",
						"Affymetrix Expression"="af_exp",
						"Illumina Methylation"="il_methyl",
						"RNA-Seq"="rna_seq"
					),
					selected="ag_exp2"
				),
				actionButton("submit_array_type", "Begin")
			))
		)
	),
	dashboardPage(
	dashboardHeader(title="eUTOPIA: solUTion for Omics data PreprocessIng and Analysis", titleWidth="35%"),
	dashboardSidebar(disable=FALSE,
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
                                                #shinyDirButton("dirButton", label="Browse", title="Select Directory", buttonType="default"),
						#shinyBS::bsTooltip("dirButton", "Browse local system directories and select the directory containing the RAW data files!", placement="bottom")
                                                div(id="outerDgxDiv", class="form-group shiny-input-container",
                                                        tags$label("Select Directory"),
                                                        div(id="innerDgxDiv", class="input-group",
                                                                div(id="dgxBtnDiv", class="input-group-btn",
                                                                        shinyDirButton("dirButton", label="Browse", title="Select Directory", buttonType="default"),
                                                                        shinyBS::bsTooltip("dirButton", "Browse local system directories and select the directory containing the RAW data files!", placement="bottom")
                                                                ),
                                                                tags$input(id="dirTextDisp", class="form-control", placeholder="No Directory", readonly="readonly", type="text")
                                                        )
                                                )
                                        )
                                #),fluidRow(
                                #        column(1, align="center",
                                #                p("")
                                #        ),column(11, align="left",
                                #                textOutput("dirText")
                                #        )
                                ),fluidRow(
					column(12,
						hr(),
						div(id="agSourceDiv",
							fluidRow(
                                                                column(12,
                                                                        selectInput("selAgSource", "Source Image Analysis Program",
                                                                            choices=c(
                                                                                "agilent",
                                                                                "agilent.median",
                                                                                "agilent.mean",
                                                                                "genepix",
                                                                                "genepix.median"
                                                                            ),
                                                                            selected="agilent.median"
                                                                        )
                                                                )
                                                        )
						),
						hidden(div(id="affAnnDiv", class="contentDiv",
							h4("Select Annotation"),
							fluidRow(
								column(1, align="center",
									p("")
								),column(11,
									actionButton("launch_cdf_modal", label="Install Missing CDF Annotation")
								)
							),fluidRow(column(12,
									uiOutput("selAffCDF")
								)
							)
						)),
						hidden(div(id="seqAnnDiv",
							fluidRow(column(12,
									checkboxInput(inputId="chkPaired", label="Paired End", value=FALSE),
									uiOutput("selCoresSeq"),
									h4("Load Annotation"),
									fileInput("fGTF", label="GTF Annotation File")
								)
							)
						))
					)
                                ),fluidRow(
                                        #column(1, align="center",
                                        #        p("")
                                        #),column(11, align="left",
                                        column(12, align="center",
                                                actionButton("upload_raw_submit", "Upload")
                                        )
                                )
                        ),bsCollapsePanel("QUALITY CONTROL", style="danger",
                                fluidRow(
                                        column(12, align="center",
                                        #column(12,
                                                div(id="download_QC_div",
							hidden(div(id="affy_QC_div",
                                                                fluidRow(
                                                                        column(12, align="center",
                                                                                uiOutput("selPoolSize")
                                                                        ),column(12, align="center",
                                                                                checkboxInput(inputId="chkParallel", label="Multithreaded", value=FALSE),
                                                                                #shinyjs::disabled(uiOutput("selCores"))
                                                                                uiOutput("selCores")
                                                                        )
                                                                )
                                                        )),
                                                        downloadButton("exportQC", "QC Report")
                                                ),
                                                hidden(div(id="methyl_QC_div",
                                                        actionButton("qc_methyl_submit", "QC Report")
                                                ))
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                actionButton("qc_skip_submit", "Skip QC")
                                        )
                                )
                        ),bsCollapsePanel("FEATURE FILTERING", style="danger",
                                fluidRow(
                                        column(12,
                                                sliderInput("filtDist", "Quantile Based Cutoff", min=0.10, max=1.0, value=0.75, step=0.05, round=2),
                                                selectInput("detectPV", "P.value detection threshold", choices=c(0.01, 0.05), selected=0.01), #For Illumina Methylation
                                                sliderInput("perSamples", "Percentage of Samples", min=1, max=100, value=75, step=1),
						uiOutput("selVarIFilt"), #For RNA-Seq
                                                uiOutput("selFiltMethod"), #For RNA-Seq
						numericInput("filtCPM", "Counts Per Million", value=1, min=1, max=NA), #For RNA-Seq
                                                selectInput("filtPvAdj", "P.value Adjustment",
							choices=c(
								"Holm"="holm",
								"Hochberg"="hochberg",
								"Hommel"="hommel",
								"Bonferroni"="bonferroni",
								"Benjamini & Hochberg"="BH",
								"Benjamini & Yekutieli"="BY",
								"False Detection Rate"="fdr",
								"None"="none"
							),
							selected="none"
						), #For RNA-Seq
                                                hr()
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                actionButton("filt_submit", "Filter"),
                                                hr()
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
                                                uiOutput("selNormMethod2"),
                                                hidden(div(id="ilChkDiv",
                                                        checkboxInput(inputId="chkRmSNP", label="Remove SNPs", value=FALSE),
                                                        checkboxInput(inputId="chkRmSexChr", label="Remove Sex Chromosomes", value=FALSE)
                                                ))
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                actionButton("norm_submit", "Run Normalization")
                                        )
                                )
                        ),bsCollapsePanel("BATCH CORRECTION", style="danger",
                                fluidRow(
                                        column(12,
                                                ##uiOutput("selVarI"),
                                                ##uiOutput("selCoVar"),
                                                ##uiOutput("selBatch"),
                                                ##checkboxInput("rmSva", "Remove Unknown Batches", value=TRUE),
						##shinyBS::bsTooltip("rmSva", "Hidden batch variables are identified and removed. These are not present in the user provided phenotype data and are not associated with the 'Variables of Interest' or 'Co-variates'!", placement="top")
                                                #selectInput("corrType", "Correction With", choices=c("SVA+ComBat"="sc", "ComBat"="c", "Nothing"="n"), selected="sc")
                                                selectInput("corrType", "Correction With", choices=c("SVA+ComBat"="sc", "ComBat"="c", "SVA"="s", "Nothing"="n"), selected="n")
                                        )
                                ),fluidRow(
                                        column(12, align="left",
                                                div(id="svaDiv",
                                                        shinyBS::bsButton("launch_sva_modal", label="Launch SVA Module", style="danger", icon=icon("exclamation-circle")),
                                                        #shinyBS::bsTooltip("launch_sva_modal", "Launch a graphical window, to configure and run Surrogate Variable Analysis (SVA)!", placement="bottom"),
                                                        p(" ")
                                                )
                                        )
                                ),fluidRow(
                                        column(12, align="left",
                                                div(id="combatDiv",
                                                        shinyBS::bsButton("launch_combat_modal", label="Launch ComBat Module", style="danger", icon=icon("exclamation-circle")),
                                                        #shinyBS::bsTooltip("launch_combat_modal", "Launch a graphical window, to configure and run batch correction with ComBat!", placement="bottom"),
                                                        p(" ")
                                                )
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                div(id="skipDiv",
                                                        actionButton("corr_skip_submit", "SKIP")
                                                )
                                        )
                                )
                        ),bsCollapsePanel("ANNOTATION", style="danger",
                                fluidRow(
                                        column(12, align="left",
                                                selectInput("annType", "Select Annotation Type",
                                                        choices=c(
                                                                "From Raw Data"="raw",
                                                                "Annotation File"="file"
                                                        ),
                                                        selected="file"
                                                )
                                        )
                                ),fluidRow(
                                        column(12,
                                                div(id="launchAnnModalDiv",
                                                        fluidRow(column(12, align="left",
                                                                shinyBS::bsButton("launch_ann_modal", label="Annotate", style="danger", icon=icon("exclamation-circle")),
                                                                p(" ")
                                                        ))
                                                )
                                        )
                                ),fluidRow(
                                        column(12,
                                                hidden(div(id="submitAnnDiv",
                                                        fluidRow(column(12, align="left",
                                                                shinyBS::bsButton("submit_ann", label="Annotate", style="success", icon=icon("angle-right")),
                                                                p(" ")
                                                        ))
                                                ))
                                        )
                                )
                        ),bsCollapsePanel("DIFFERENTIAL ANALYSIS", style="danger",
                                fluidRow(
																			hidden(div(id="diffMethodDiv",
																							fluidRow(
																									column(12, align="left",
																											selectInput("repType", "Replicates",
																													choices=c(
																															"Technical"="technical",
																															"Biological"="biological"
																													),
																													selected=c("noiseq", "edger"),
																													multiple=TRUE
																											)
																									)
																							),fluidRow(
																									column(12, align="left",
																											selectInput("diffMethod", "Method",
																													choices=c(
																															"NOISeq"="noiseq",
																															"edgeR"="edger"
																													),
																													selected="noiseq"
																											) #For RNA Seq
																									)
																							),fluidRow(
																									column(12, align="left",
																											selectInput("joinType", "Combine Methods",
																													choices=c(
																															"Union"="full",
																															"Intersect"="inner"
																													),
																													selected="full"
																											) #For RNA Seq
																									)
																							)
																			))
																	),fluidRow(
                                        column(12,
                                                div(id="limmaModelDiv", class="contentDiv",
                                                        h4("Design Model"),
                                                        fluidRow(
                                                                column(12,
                                                                        uiOutput("selVarILimma"),
                                                                        uiOutput("selCoVarLimma"),
                                                                        checkboxInput(inputId="chkKeepNA", label="Keep NA Variables", value=FALSE)
                                                                        #uiOutput("selBatchLimma")
                                                                )
                                                        )
                                                )
                                        )
                                ),fluidRow(
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
                        ),bsCollapsePanel("REPORTING", style="info",
                                fluidRow(
                                        column(12,
                                                div(id="exportExprDiv", class="contentDiv",
                                                        h4("Export Expression Matrix"),
                                                        fluidRow(column(12, align="center",
                                                                downloadButton("exportFiltMat", "Filtered"),
                                                                p(" ")
                                                        )),fluidRow(column(12, align="center",
                                                                downloadButton("exportNormMat", "Normalized"),
                                                                p(" ")
                                                        )),fluidRow(column(12, align="center",
                                                                downloadButton("exportCorrMat", "Corrected"),
                                                                p(" ")
                                                        )),fluidRow(column(12, align="center",
                                                                downloadButton("exportAggMat", "Aggregated")
                                                        ))
                                                ),
                                                p(" ")
                                        )
                                ),fluidRow(
                                        column(12, align="center",
                                                downloadButton("exportRpt", "Analysis Report")
                                        )
                                )
                        )
                )
         ),
	dashboardBody(
		shinyBS::bsModal("importPhenoModal", "Import Phenotype Data", "import_pheno_submit", size="large",
			fluidRow(
				column(3,
					fileInput("fPheno", label="Phenotype File")
				),column(3,
					uiOutput("selSep")
				),column(3,
					textInput("sepT", "Other Seperator", value=":")
				),column(3,
					uiOutput("selQuote")
				)
			),fluidRow(
				column(1,
					actionButton("load_pheno_submit", "Preview")
				),column(2, align="left",
					textOutput("phRowsText"),
					textOutput("phColsText")
				)
			),hr(),
			fluidRow(
				column(12,
					hidden(div(id="phenoPreviewDiv",
                                                fluidRow(
                                                        column(12,
                                                                rhandsontable::rHandsontableOutput("phenoTypesRH")
                                                        )
                                                ),hr(),
                                                fluidRow(
                                                #        column(12,
                                                #                DT::dataTableOutput("phenoDT")
                                                #        )
                                                #),fluidRow(
                                                        column(4,
                                                                uiOutput("selFileNameCol")
                                                        ),column(4,
                                                                uiOutput("selSampleIDCol")
                                                        ),column(4,
                                                                uiOutput("selDyeCol")
                                                        )
                                                ),fluidRow(
                                                        column(12, align="right",
                                                                shinyBS::bsButton("upload_pheno_submit", label="Import", style="info", icon=icon("hand-o-right"))
                                                        )
                                                )
                                        ))
                                )
                        )
		),
		shinyBS::bsModal("importAnnModal", "Import Annotation", "launch_ann_modal", size="large",
			fluidRow(
				column(12,
				#	selectInput("annType", "Select Annotation Type",
				#		choices=c(
				#			"Ensembl Biomart"="mart",
				#			"Annotation File"="file"
				#		),
				#		selected="file"
				#	)
				#),column(12,
					#hidden(div(id="annFile",
					div(id="annFile",
                                                fluidRow(
                                                        column(4,
                                                                fileInput("fAnno", label="Microarray Annotation")
                                                        ),column(4,
                                                                uiOutput("selSepAnno")
                                                        ),column(4,
                                                                textInput("sepTAnno", "Other Seperator", value=":")
                                                        )
                                                ),fluidRow(
                                                        column(12,
                                                                #actionButton("load_ann_submit", "Import")
								#shinyBS::bsButton("load_ann_submit", label="Import", style="info", icon=icon("hand-o-right"))
                                                                actionButton("load_ann_submit", "Preview")
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
					#)),div(id="annQuery",
					#	selectInput("org", "Select Organism",
					#		choices=c(
					#			"Human"="hsapiens_gene_ensembl",
					#			"Mouse"="mmusculus_gene_ensembl"
					#		),
					#		selected="hsapiens_gene_ensembl"
					#	),uiOutput("selProbeID"),
					#	selectInput("mapID", "Select Mapping ID",
					#		choices=c(
					#			"Ensembl Gene ID"="ensembl_gene_id",
					#			"Ensembl Transcript ID"="ensembl_transcript_id",
					#			"EntrezGene ID"="entrezgene",
					#			"Gene Symbol"="hgnc_symbol"
					#		),
					#		selected="ensembl_gene_id"
					#	)
					#	#uiOutput("mapID")
					)
				)
                        ),fluidRow(
				column(12, align="right",
					shinyBS::bsButton("upload_ann_submit", label="Import", style="info", icon=icon("hand-o-right"))
				)
			),hr(),
			fluidRow(
				column(12,
					DT::dataTableOutput("annDT")
				)
			)
		),
		shinyBS::bsModal("installCDFModal", "Import Annotation", "launch_cdf_modal", size="large",
			fluidRow(
				column(12,
					textInput("cdfURL", "URL of the CDF source R package", value=""),
					actionButton("install_cdf_submit", label="Install")
				)
			)
		),
                shinyBS::bsModal("svaModal", "Configure and Run SVA", "launch_sva_modal", size="large",
			fluidRow(
				column(12,
                                        uiOutput("selVarISva"),
                                        uiOutput("selCoVarSva"),
                                        #uiOutput("selBatchSva"),
                                        actionButton("sva_submit", "Run SVA")
				)
			)
		),
                shinyBS::bsModal("combatModal", "Configure and Run ComBat", "launch_combat_modal", size="large",
			fluidRow(
				column(12,
                                        uiOutput("selVarICombat"),
                                        uiOutput("selCoVarCombat"),
                                        uiOutput("selBatchCombat"),
                                        actionButton("combat_submit", "Run ComBat")
				)
			)
		),
		#div(id="loading-content",
		#	img(id="loading-gif", src="screen-loading.gif")
		#),
		fluidRow(
			column(6,
				#tags$head(tags$style(HTML(".small-box {height: 70px}"))),
				infoBoxOutput('totalSampleBox'),
				infoBoxOutput('filteredSampleBox'),
				infoBoxOutput('removedSampleBox')
			#),column(6,
			##	valueBoxOutput('varIValueBox'),
			##	valueBoxOutput('coVarValueBox'),
			##	valueBoxOutput('batchValueBox')
                        #        div(id="buttonsDiv", class="contentDiv",
                        #                fluidRow(
                        #                        column(4,
                        #                                actionButton("load_expr_submit", "Load Expression Data")
                        #                        ),column(4,
                        #                                downloadButton("exportQC", "QC Report")
                        #                        ),column(4,
                        #                                #actionButton("analysis_report_submit", "Analysis Report")
                        #                                downloadButton("exportRpt", "Analysis Report")
                        #                        )
                        #                ),fluidRow(
                        #                        column(4,
                        #                                downloadButton("exportMat", "Export Expression Matrix")
                        #                        )
                        #                )
                        #        )
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
						tabBox(id="normTBox", title="", width=12,
							tabPanel(value="boxPlotTab", title="Box Plot",
								fluidRow(
									column(6,
										shinycssloaders::withSpinner(plotOutput('preBoxPlot'), type=6)
									),column(6,
										shinycssloaders::withSpinner(plotOutput('postBoxPlot'), type=6)
									)
								)
							),tabPanel(value="densityPlotTab", title="Density Plot",
								fluidRow(
									column(6,
										shinycssloaders::withSpinner(plotOutput('preDensityPlot'), type=6)
									),column(6,
										shinycssloaders::withSpinner(plotOutput('postDensityPlot'), type=6)
									)
								)
							),tabPanel(value="MDPlotTab", title="Mean-Difference Plot",
								fluidRow(
									column(6,
										shinycssloaders::withSpinner(plotOutput('preMDPlot'), type=6)
									),column(6,
										shinycssloaders::withSpinner(plotOutput('postMDPlot'), type=6)
									)
								)
							)
						)
					))
				),tabPanel(value="tvTab", title="Technical Variation",
					fluidRow(column(12,
						tabBox(id="mtv", title="", width=12,
							tabPanel(value="confPlotTab", title="Confounding Plot",
								fluidRow(
									column(12,
										#div(id="confPlotDiv", class="sizeable",
										div(id="confPlotDiv", class="jacket",
											shinycssloaders::withSpinner(plotOutput('confPlot', height="auto", width="auto"), type=6)
										)
									)
								)
							),tabPanel(value="princePlotTab", title="Prince Plot",
                                                                fluidRow(
                                                                        column(12,
                                                                                tabBox(id="princeBox", title="", width=12,
                                                                                        tabPanel(value="prePrincePlotTab", title="Before Correction",
                                                                                                fluidRow(
                                                                                                        column(12,
                                                                                                                #div(id="princePlotDiv", class="sizeable",
                                                                                                                div(id="princePlotDiv", class="jacket",
                                                                                                                        shinycssloaders::withSpinner(plotOutput('princePlot', height="auto", width="auto"), type=6)
                                                                                                                )
                                                                                                        )
                                                                                                )
                                                                                        ),tabPanel(value="postPrincePlotTab", title="After Correction",
                                                                                                fluidRow(
                                                                                                        column(12,
                                                                                                                #div(id="postPrincePlotDiv", class="sizeable",
                                                                                                                div(id="postPrincePlotDiv", class="jacket",
                                                                                                                        shinycssloaders::withSpinner(plotOutput('postPrincePlot', height="auto", width="auto"), type=6)
                                                                                                                )
                                                                                                        )
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                )
							),tabPanel(value="hcTab", title="Hierarchical Clustering",
								fluidRow(
									column(12,
										#div(id="hcPlotDiv", class="plot-container sizeable",
                                                                                #        tags$img(src="spinner.gif", id="loading-spinner"),
										#div(id="hcPlotDiv", class="sizeable",
                                                                                tabBox(id="hcBox", title="", width=12,
                                                                                        tabPanel(value="preHcTab", title="Before Correction",
                                                                                                fluidRow(
                                                                                                        column(12,
                                                                                                                div(id="hcPlotDiv", class="jacket",
                                                                                                                        shinycssloaders::withSpinner(plotOutput('hcPlot', height="auto", width="auto"), type=6)
                                                                                                                )
													)
                                                                                                )
                                                                                        ),tabPanel(value="postHcTab", title="After Correction",
                                                                                                fluidRow(
                                                                                                        column(12,
                                                                                                                #div(id="postHcPlotDiv", class="sizeable",
                                                                                                                div(id="postHcPlotDiv", class="jacket",
                                                                                                                        shinycssloaders::withSpinner(plotOutput('postHcPlot', height="auto", width="auto"), type=6)
                                                                                                                )
													)
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                )
                                                        ),tabPanel(value="mdsTab", title="MDS Plots",
                                                                fluidRow(
                                                                        column(2,
                                                                                uiOutput("selMdsLabel")
                                                                        ),column(2, align="left",
                                                                                uiOutput("selMdsColor")
                                                                        )
                                                                ),fluidRow(column(12,
                                                                        tabBox(id="mdsBox", title="", width=12,
                                                                                tabPanel(value="preCorTab", title="Before Correction",
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="preCorMDSDiv", class="sizeable",
													div(id="preCorMDSDiv", class="jacket",
														shinycssloaders::withSpinner(plotOutput('preCorMDS', height="auto", width="auto"), type=6)
													)
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postCorTab", title="After Correction",
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postCorMDSDiv", class="sizeable",
													div(id="postCorMDSDiv", class="jacket",
														shinycssloaders::withSpinner(plotOutput('postCorMDS', height="auto", width="auto"), type=6)
													)
                                                                                                )
                                                                                        )
                                                                                ),tabPanel(value="postAggTab", title="After Aggregation",
                                                                                        fluidRow(
                                                                                                column(12,
													#div(id="postAggMDSDiv", class="sizeable",
													div(id="postAggMDSDiv", class="jacket",
														shinycssloaders::withSpinner(plotOutput('postAggMDS', height="auto", width="auto"), type=6)
													)
                                                                                                )
                                                                                        )
                                                                                )
                                                                        )
                                                                ))
							)
						)
					))
				#),tabPanel(value="diffTab", title="Differential Expression",
				),tabPanel(value="diffTab", title="Differential Results",
					fluidRow(column(3,
                                                        numericInput(inputId="lfcThr", "LogFC Threshold", value=0, min=0, max=20, step=0.1)
                                                ),column(3,
                                                        #numericInput(inputId="adjPvThr", paste0("Adj. P.Value Threshold"), value=0.05, min=0, max=1, step=0.001)
                                                        numericInput(inputId="pvThr", paste0("P.Value Threshold"), value=0.05, min=0, max=1, step=0.001)
                                                ),column(3,
                                                        selectInput("pvType", "P.Value Type", choices=c("Raw"="P.Value", "Adjusted"="adj.P.Val"), selected="adj.P.Val")
                                                ),column(3,
                                                        uiOutput('selCompDeTable')
                                                #),column(2,
                                                #        downloadButton("exportRpt", "Generate Report")
                                                #),column(2,
                                                #        downloadButton("exportMat", "Export Expression Matrix")
                                                )
					),fluidRow(column(9,
                                                        downloadButton("exportDE", "Export Differential Tables"),
                                                        checkboxInput(inputId="chkExportDE", label="Export Unfiltered Tables", value=FALSE)
                                                ),column(3,
                                                        actionButton("filterDE_submit", "Filter Differential Result")
                                                )
                                        ),
					fluidRow(column(12,
						tabBox(id="diffTBox", title="", width=12,
							#tabPanel(value="diffSummTableTab", title="Differential Expression Summary",
							tabPanel(value="diffSummTableTab", title="Differential Summary",
								fluidRow(
									column(12,
										DT::dataTableOutput("deSummTable")
                                                                        )
								)
							#),tabPanel(value="diffTableTab", title="Differential Expression Tables",
							),tabPanel(value="diffTableTab", title="Differential Tables",
								fluidRow(
									column(12,
										DT::dataTableOutput("deTable")
                                                                        )
								)
							#),tabPanel(value="diffVenn", title="Differential Gene Sets Intersection",
							),tabPanel(value="diffVenn", title="Differential Sets Intersection",
								fluidRow(
									column(4,
                                                                                uiOutput("selIntersectComps"),
                                                                                actionButton("venn_submit", "Replot")
									),column(8,
                                                                                ##div(id="intersectPlotDiv", class="sizeable",
                                                                                #div(id="intersectPlotDiv", class="jacket",
                                                                                #        shinycssloaders::withSpinner(plotOutput('intersectPlot', height="auto", width="auto"), type=6)
                                                                                #)
                                                                                div(id="vennContentDiv", class="contentDiv", align="center",
                                                                                        h4("Venn Representation"),
                                                                                        div(id="vennPlotDiv", class="jacket",
                                                                                                shinycssloaders::withSpinner(plotOutput('vennPlot', height="auto", width="auto"), type=6)
                                                                                        )
                                                                                ),
                                                                                h4(" "),
                                                                                div(id="upsetContentDiv", class="contentDiv", align="center",
                                                                                        h4("UpSet Representation"),
                                                                                        div(id="upsetPlotDiv", class="jacket",
                                                                                                shinycssloaders::withSpinner(plotOutput('upsetPlot', height="auto", width="auto"), type=6)
                                                                                        )
                                                                                )
									)
								)
							),tabPanel(value="diffVolcano", title="Volcano Plot",
								fluidRow(
									column(3,
										sliderInput("xAxis", "X Axis Bounds", min=1, max=7, value=5, step=1)
									),column(3,
										sliderInput("yAxis", "Y Axis Bounds", min=1, max=30, value=20, step=1)
									),column(3,
                                                                                actionButton("vol_submit", "Replot")
									)
								),fluidRow(
									column(12,
										#checkboxInput("gName", "Display Gene Names", value=FALSE)
										checkboxInput("gName", "Display Names", value=FALSE)
									)
								),fluidRow(
									column(12,
										shinycssloaders::withSpinner(plotOutput('volcanoPlot', height="auto", width="auto"), type=6)
									)
								)
							)
                                                )
					))
				),tabPanel(value="visExp", title="Visualize Expression/Methylation",
					fluidRow(column(12,
						tabBox(id="visExpTBox", title="", width=12,
                                                        tabPanel(value="expBoxplot", title="Box Plot",
                                                                fluidRow(
                                                                        column(4,
                                                                                #uiOutput("selExpGenes"),
                                                                                selectizeInput("expGenes", "Select Features", choices=c("NA"), multiple=TRUE),
                                                                                uiOutput("selVarIBoxPlot"),
                                                                                uiOutput("selConditionsBoxplot"),
                                                                                actionButton("expBoxplot_submit", "Plot")
                                                                        ),column(8,
                                                                                shinycssloaders::withSpinner(plotOutput("expressionBoxPlot", height="auto", width="auto"), type=6)
                                                                        )
                                                                )
                                                        ),tabPanel(value="expHeatmap", title="Heatmap",
                                                                fluidRow(
                                                                        column(4,
                                                                                #uiOutput("slidePercDE"),
                                                                                sliderInput("percDE", "Percentage of Highly Differential Features", min=0.1, max=5, value=0.1, step=0.1),
                                                                                uiOutput("selHeatComps"),
                                                                                uiOutput("selConditionsHeat"),
                                                                                actionButton("expHeatmap_submit", "Replot")
                                                                        ),column(8,
                                                                                shinycssloaders::withSpinner(plotOutput("expressionHeatmap", height="auto", width="auto"), type=6)
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
)
