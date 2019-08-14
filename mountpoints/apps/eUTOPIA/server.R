suppressMessages(library(shiny))
suppressMessages(library(DT))
suppressMessages(library(rhandsontable))
suppressMessages(library(shinyjs))
suppressMessages(library(shinyBS))
suppressMessages(library(shinydashboard))
suppressMessages(library(shinyFiles))
suppressMessages(library(ggplot2))
suppressMessages(library(ggrepel))
suppressMessages(library(rmarkdown))
suppressMessages(library(grid))
suppressMessages(library(UpSetR))
suppressMessages(library(reshape2))
suppressMessages(library(gplots))
suppressMessages(library(htmlTable))
suppressMessages(library(devtools))
suppressMessages(library(httr))
suppressMessages(library(randomcoloR))
suppressMessages(library(Biobase))

#Setup options for local library path
.libPaths("/home/shiny/Rlibs")

##Get source directory
#srcDir <- dirname(getSrcDirectory(function(x){x}))
#print("Print Source Directory")
##print(dirname(getSrcDirectory(function(x){x})))
#print(srcDir)
#array_pipeline_R <- file.path(srcDir, "agilent_twocolors_array_analysis_pipe_fixed.R")

srcDir <- "/srv/shiny-server/eUTOPIA/"
print("Print Source Directory")
print(srcDir)
array_pipeline_R <- file.path(srcDir, "agilent_twocolors_array_analysis_pipe_fixed.R")
source(array_pipeline_R)

options(shiny.maxRequestSize=500*1024^2)
set.seed(1) #Block randomness. Set seed

celTable <- NULL
phTable <- NULL
extractedList <- NULL

#Function to hide bsCollapsePanel
hideBSCollapsePanel <- function(session, panel.name)
{
	#print("In hideBSCollapsePanel")
	session$sendCustomMessage(type='hideBSCollapsePanelMessage', message=list(name=panel.name))
	#print("Exiting hideBSCollapsePanel")
}

#Function to get names of columns with datatype character
factorize_cols <- function(phTable, idx){
        for(i in idx){
                phTable[,i] <- factor(phTable[,i])
        }
        return(phTable)
}

shinyServer(
	function(input, output, session){
                gVars <- shiny::reactiveValues(
                        phTable=NULL,
                        rgList=NULL,
                        celDir=NULL,
                        totalSamples=NULL,
                        filteredSamples=NULL,
                        removedSamples=NULL,
                        norm.data=NULL,
                        pcChoices=NULL,
                        comb.data=NULL,
                        agg.data=NULL,
                        comps=list(),
                        loadedRaw=FALSE,
                        QC_passed=FALSE,
                        filtered=FALSE,
                        normalized=FALSE,
                        corrected=FALSE,
			filePaths=NULL, #New loaded raw features
			rawDir=NULL
                )

		gVars$sepChoices <- c("TAB", ",", ";", "SPACE", "OTHER")
		gVars$quoteChoices <- c(NA, "SINGLE", "DOUBLE")
		gVars$pvAdjChoices <- c("Holm"="holm", "Hochberg"="hochberg", "Hommel"="hommel", "Bonferroni"="bonferroni", "Benjamini & Hochberg"="BH", "Benjamini & Yekutieli"="BY", "False Detection Rate"="fdr", "None"="none")
		gVars$normChoices <- c("Between Arrays"="BA", "Quantile"="quantile", "Variance Stabilizing"="vsn", "Cyclic Loess"="cl")
		gVars$baChoices <- c("None"="none", "Scale"="scale", "Quantile"="quantile", "Cyclic Loess"="cyclicloess")



                ## Get available Affymetrix CDF annotations
		allPkgNames <- rownames(installed.packages())
		gVars$cdfPkgNames <- allPkgNames[grep(".*cdf", allPkgNames)]
		cdfQC <- function(x){if(length(x)>0){return(x)}else{return("NA")}}
		gVars$cdfChoices <- reactive({
			cdfPkgNames <- gVars$cdfPkgNames
			cdfQC(cdfPkgNames)
		})

		observeEvent(input$install_cdf_submit, {
			if(is.null(input$cdfURL) || input$cdfURL==""){
				return(NULL)
			}
			cdfURL <- input$cdfURL
			if(httr::http_error(cdfURL)){
				shinyjs::info(
					paste0("Invalid OR Expired URL!\n\nNo Installation Performed!\n\nURL: ", cdfURL)
				)
				return(NULL)
			}
			progress <- shiny::Progress$new()

			#start loading screen
			shinyjs::html(id="loadingText", "Installing CDF annotation")
			shinyjs::show(id="loading-content")

                        on.exit({
				progress$close()
				shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
			})
                        progress$set(message="CDF Package...", value=0)
                        progress$set(detail="Installing...", value=1/2)
			#devtools::install_url(cdfURL)
			withr::with_libpaths("/home/shiny/Rlibs", devtools::install_url(cdfURL))
			Sys.sleep(2)

                        progress$set(detail="Updating...", value=2/2)
			allPkgNames <- rownames(installed.packages())
			gVars$cdfPkgNames <- allPkgNames[grep(".*cdf", allPkgNames)]
		})

		#New load raw data features
		observeEvent(input$load_raw, {
			print("In loaded raw...")
			if(is.null(input$raw)){
				return(NULL)
			}

			if(is.null(gVars$rawDir)){
				rawDir <- tempdir()
				gVars$rawDir <- rawDir
			}else{
				rawDir <- gVars$rawDir
			}

			fileNameColName <- gVars$fileNameColName
			rawDF <- input$raw
			rawFileNames <- rawDF$name
			rawFilePaths <- rawDF$datapath
			print("rawFileNames:")
			print(rawFileNames)
			print("rawFilePaths:")
			print(rawFilePaths)
			
			rawFilePathsCopy <- file.path(rawDir, rawFileNames)
			file.copy(rawFilePaths, rawFilePathsCopy)
			if(any(file.exists(rawFilePaths))){
				file.remove(rawFilePaths[file.exists(rawFilePaths)])
			}
			gVars$filePaths <- setNames(object=rawFilePathsCopy, nm=rawFileNames)
			print("str(gVars$filePaths)")
			print(str(gVars$filePaths))
			print(paste0("length: ", length(gVars$filePaths)))

			phTable <- gVars$phTable

			fileNames <- as.vector(phTable[,fileNameColName])
                        if(input$arrType=="il_methyl"){
				idxLoaded <- NULL
				for(i in seq(length(fileNames))){
					fileName <- fileNames[i]
					fileNameIdats <- paste0(fileName, c("_Grn.idat", "_Red.idat"))
					print("fileNameIdats")
					print(fileNameIdats)
					idxTmp <- which(fileNameIdats %in% names(gVars$filePaths))
					print("idxTmp")
					print(idxTmp)
					if(length(idxTmp)==2){
						idxLoaded <- append(idxLoaded, i)
					}
				}
			}else{
				idxLoaded <- which(as.vector(phTable[,fileNameColName]) %in% names(gVars$filePaths))
			}

			if(length(idxLoaded)>0){
				print("length(idxLoaded)")
				print(length(idxLoaded))
				gVars$idxLoaded <- union(gVars$idxLoaded, idxLoaded)
				print("length(gVars$idxLoaded)")
				print(length(gVars$idxLoaded))
			}

			shinyBS::toggleModal(session, "loadRawModal", toggle="open")
		})

		output$samplesLoadedText <- renderText({
                        if(is.null(gVars$phTable)){
				str <- "(NA/NA)"
			}else{
				phTable <- gVars$phTable
				fileNameColName <- gVars$fileNameColName
				samplesTotal <- length(as.vector(phTable[,fileNameColName]))
				samplesLoaded <- ifelse(is.null(gVars$idxLoaded), 0, length(gVars$idxLoaded))
				str <- paste0("(", samplesLoaded, "/", samplesTotal, ")")
			}
			
			return(paste0("Samples - ", str))
		})

		output$sampleLoadedDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$phTable), "No phenotype file!")
			)
			
			phTable <- gVars$phTable
			fileNameColName <- gVars$fileNameColName
			samplesTotal <- length(as.vector(phTable[,fileNameColName]))
			samplesLoaded <- ifelse(is.null(gVars$idxLoaded), 0, length(gVars$idxLoaded))
			samplesMissing <-  samplesTotal - samplesLoaded
		
			samplesLoadedDF <- data.frame("Total"=samplesTotal, "Loaded"=samplesLoaded, "Missing"=samplesMissing, row.names="Samples")

			DT::datatable(samplesLoadedDF, filter="none",
				options = list(
					paging=FALSE,
					searching=FALSE,
					info=FALSE
				)
			)
		})

		output$rawLoadedDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$phTable), "No phenotype file!")
			)
			fileNameColName <- gVars$fileNameColName
			sampleColName <- gVars$sampleColName
			phTable <- gVars$phTable[,c(sampleColName,fileNameColName)]
			phTable <- data.frame(Loaded="NO", phTable, stringsAsFactors=FALSE)
			if(isFALSE(is.null(gVars$idxLoaded)) && length(gVars$idxLoaded)>0){
				phTable$Loaded[gVars$idxLoaded] <- "YES"
			}
			DT::datatable(phTable, filter="none",
				options = list(
					search=list(regex=TRUE, caseInsensitive=FALSE),
					scrollX=TRUE,
					ordering=F
				)
			) %>% DT::formatStyle(
				'Loaded',
				backgroundColor = DT::styleEqual(c("YES", "NO"), c('#4CA079', '#DD4B39'))
			)
		})

                observeEvent(input$load_ann_submit, {
			if(is.null(input$fAnno))
			return(NULL)

                        correctionLvl <- gVars$correctionLvl
                        if(is.null(correctionLvl)){
                                expr.data <- gVars$norm.data
                        }else if(correctionLvl==1){
                                expr.data <- gVars$comb.data
                        }else{
                                shinyjs::info("Could not locate expression data!")
                                return(NULL)
                        }
                        if(is.null(expr.data)){
                                shinyjs::info("Expression data is NULL!")
                                return(NULL)
                        }else if(nrow(expr.data)==0){
                                shinyjs::info("Expression data is empty!")
                                return(NULL)
                        }

                        gVars$annDF <- NULL
			annFile <- input$fAnno
			sepS <- input$sepSAnno
			sepT <- input$sepTAnno
			sepChar <- NULL
			if(sepS=="OTHER"){
				sepChar <- sepT
			}else{
				if(sepS=="TAB"){
					sepChar <- "\t"
				}else if(sepS=="SPACE"){
					sepChar <- " "
				}else{
					sepChar <- sepS
				}
			}

                        con <- file(annFile$datapath, "r", blocking=FALSE)
                        fileLines <- readLines(con)
			fileLines <- gsub("\t\t", "\tNA\t", fileLines)
			fileLines <- gsub("\t$", "\tNA", fileLines)
                        close(con)
                        colNumByRowDist <- table(sapply(fileLines, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=F))
                        if(length(colNumByRowDist) > 1){
                                shinyjs::info(paste0("Separating character '", sepChar, "', results in inconsistent number of columns!\n\nPlease check the input file format and select the correct separating character!"))
                                return(NULL)
                        }

                        annDF <- read.csv(annFile$datapath, header=TRUE, sep=sepChar, stringsAsFactors=FALSE, quote="")

                        chkIdx <- which(rownames(expr.data) %in% annDF[,1])
                        if(length(chkIdx) == 0){
                                shinyjs::info(paste0("No ProbeName from the raw data could be mapped to annotation file '", input$fAnno, "'\n\n, Please check and provide appropriate annotation file!"))
                                return(NULL)
                        }

                        chkPerc <- (length(chkIdx)/nrow(expr.data))*100
                        if(chkPerc < 85){
                                shinyjs::info(paste0("Mapped only ", chkPerc, "% of ProbeNamesfrom the raw data to annotation file '", input$fAnno, "'\n\n, Please check and provide appropriate annotation file!"))
                                return(NULL)
                        }

			gVars$annLoaded <- 1
			chkIdx <- which(annDF[,1] %in% rownames(expr.data))
                        mis <- apply(annDF[chkIdx,-1], 2, function(x){sum(sum(is.na(x)), is.null(x),sum(x==""))})
                        percMis <- round(100-(mis/nrow(annDF)*100))
                        percMisStr <- paste0(names(percMis), "[", percMis, "%", "]")
                        dispStr <- paste0("<p>", paste(percMisStr, collapse=", "), "</p>")

                        gVars$annDF <- annDF
                        gVars$annMapDispStr <- dispStr
		})

		output$annDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$annDF), "No Annotation file!")
			)

                        annDF <- gVars$annDF
			DT::datatable(annDF, filter="none",
				options = list(
					ch=list(regex=TRUE, caseInsensitive=FALSE),
					scrollX=TRUE,
					pageLength=2,
					lengthMenu=c(1,2,3),
					ordering=F
				)
			)
		},server=TRUE)

		gVars$idChoices <- reactive({
                        annDF <- gVars$annDF
			if (is.null(annDF))
                        return("Waiting for annotation file!")

			return(colnames(annDF))
		})

		gVars$inputPh <- eventReactive(input$load_pheno_submit, {
			if(is.null(input$fPheno))
			return(NULL)

                        phFile <- input$fPheno
			sepS <- input$sepS
			sepT <- input$sepT
			sepChar <- NULL
			if(sepS=="OTHER"){
				sepChar <- sepT
			}else{
				if(sepS=="TAB"){
					sepChar <- "\t"
				}else if(sepS=="SPACE"){
					sepChar <- " "
				}else{
					sepChar <- sepS
				}
			}

			quote <- input$quote
			if(is.na(quote) || quote=="NA"){
				quote <- ""
			}else if(quote=="SINGLE"){
				quote <- "'"
			}else if(quote=="DOUBLE"){
				quote <- '"'
			}

			rowNames <- NULL
                        con <- file(phFile$datapath, "r", blocking=FALSE)
                        fileLines <- readLines(con)
			fileLines <- gsub("\t\t", "\tNA\t", fileLines)
			fileLines <- gsub("\t$", "\tNA", fileLines)
                        close(con)
                        colNumByRowDist <- table(sapply(fileLines, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=F))
                        if(length(colNumByRowDist) > 1){
				fileLines2 <- fileLines[-1]
				colNumByRowDist2 <- table(sapply(fileLines2, function(x) {length(strsplit(x, sepChar)[[1]])}, USE.NAMES=F))
				if(length(colNumByRowDist2) > 1){
					shinyjs::info(paste0("Separating character '", sepChar, "', results in inconsistent number of columns!\n\nPlease check the input file format and select the correct separating character!"))
					gVars$phLoaded <- NULL
					return(NULL)
				}
				rowNames <- 1 #First column is taken as rownames if the first row has one less value. REDUNDANT
                        }

                        #Define strings to identify NAs in the phenotype file while reading
                        naStr <- paste0(c("N", "n"), c("A", "a", "a", "A"))
                        naStr1 <- paste0("-", naStr ,"-")
                        naStr2 <- paste0("<", naStr ,">")
                        naStrAll <- c(naStr, naStr1, naStr2)
                        nullStrAll <- c("NULL", "null", "Null")

                        phTable <- read.csv(phFile$datapath, header=TRUE, sep=sepChar, stringsAsFactors=FALSE, quote=quote, as.is=TRUE, strip.white=TRUE, na.strings=c(naStrAll, "", "-", nullStrAll))
			gVars$phLoaded <- 1

                        #Get column types
                        coltypes <- unlist(lapply(phTable, class))
                        print("coltypes")
                        print(coltypes)
                        print(table(coltypes))

                        #Fix logical columns bullshit
                        coltypes.logical.idx <- which(coltypes=="logical")
                        if(length(coltypes.logical.idx)>0){
                                print("Phenotype contains logical columns!")
                                print(coltypes.logical.idx)
                                for(idx in as.vector(coltypes.logical.idx)){
                                        print("Updating logical to character!")
                                        print(colnames(phTable)[idx])
                                        phTable[,idx] <- as.character(phTable[,idx])
                                }
                                coltypes <- unlist(lapply(phTable, class))
                        }

                        coltypes.charOnly.idx <- which(coltypes=="character")
                        coltypes.nonChar.idx <- which(!coltypes=="character")
                        coltypes.charOnly.len <- length(coltypes.charOnly.idx)
                        coltypes.nonChar.len <- length(coltypes.nonChar.idx)
			remInfo <- 0
			remStr <- ""
                        if(coltypes.charOnly.len>0){
                                phTable.charOnly <- phTable[, coltypes.charOnly.idx, drop=F]
                                print("dim(phTable.charOnly)")
                                print(dim(phTable.charOnly))
                                numCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.numeric(col))))){1}else{0}}))
                                intCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.integer(col))))){1}else{0}}))
                                doubleCheck <- unlist(lapply(phTable.charOnly, function(col){if(suppressWarnings(all(is.na(as.double(col))))){1}else{0}}))
                                allCheck <- numCheck+intCheck+doubleCheck
                                #remInfo <- 0
                                #remStr <- ""
                                if(all(allCheck==0)){
                                        checkFailed <- names(allCheck[allCheck==0])
                                        remInfo <- 1
                                        remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
                                        if(coltypes.nonChar.len>0){
                                                phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
                                                phTable.comb <- phTable.nonChar
                                        }else{
                                                remStr <- paste0(remStr, "No column survived filtering!!! Please define phenotype data columns with singular data type.")
                                                shinyjs::info(remStr)
                                                return(NULL)
                                        }
                                }else{
                                        if(any(allCheck==0)){
                                                checkFailed <- names(allCheck[allCheck==0])
                                                phTable.charOnly <- phTable.charOnly[,-which(colnames(phTable.charOnly) %in% checkFailed)]
                                                remInfo <- 1
                                                remStr <- paste0("Following columns are removed because they contain mixed character and numeric data types:\n[",paste0(checkFailed, collapse=", "), "]\n\n")
                                        }
                                        print("str(phTable) -- before:")
                                        print(str(phTable))
                                        phTable.charOnly <- as.data.frame(apply(phTable.charOnly, 2, function(x){res<-trimws(x); res<-gsub(" +", " ", res, perl=T); res<-gsub("[ -]", "_", res); return(res)}), stringsAsFactors=F)
                                        if(coltypes.nonChar.len>0){
                                                phTable.nonChar <- phTable[, coltypes.nonChar.idx, drop=F]
                                                print("dim(phTable.nonChar)")
                                                print(dim(phTable.nonChar))
                                                phTable.comb <- data.frame(phTable.charOnly, phTable.nonChar, stringsAsFactors=FALSE)
                                        }else{
                                                phTable.comb <- phTable.charOnly
                                        }
                                }
                                colOrgIdx <- sapply(colnames(phTable.comb), function(x){which(colnames(phTable) %in% x)})
                                phTable <- phTable.comb[,names(colOrgIdx[order(colOrgIdx)]), drop=F]
                        }

                        print("str(phTable) -- check:")
                        print(str(phTable))
                        #Remove columns with single level data
                        nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
                        nrlevels.singular <- which(nrlevels<2)

                        if(length(nrlevels.singular)>0){
                                remInfo <- 1
                                remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
                                if(length(nrlevels.singular)==ncol(phTable)){
                                        remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
                                        shinyjs::info(remStr)
                                        return(NULL)
                                }
                                col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
                                phTable <- phTable[,-col2rem, drop=F]
                        }

                        #Inform user with the columns removed from the data frame
                        if(remInfo==1){
                                shinyjs::info(remStr)
                        }
                        print("str(phTable) -- after:")
                        print(str(phTable))

			return(phTable)
		})

                gVars$phColTypes <- reactive({
                        if(is.null(gVars$inputPh()))
                        return(NULL)

                        phTable <- gVars$inputPh()
                        colNames <- list("c"=NULL,"n"=NULL)
                        coltypes <- unlist(lapply(phTable, class))
                        coltypes.charOnly.idx <- which(coltypes=="character")
                        coltypes.nonChar.idx <- which(!coltypes=="character")
                        coltypes.charOnly.len <- length(coltypes.charOnly.idx)
                        coltypes.nonChar.len <- length(coltypes.nonChar.idx)
                        if(coltypes.charOnly.len>0){
                                colNames[["c"]] <- colnames(phTable)[coltypes.charOnly.idx]
                        }
                        if(coltypes.nonChar.len>0){
                                colNames[["n"]] <- colnames(phTable)[coltypes.nonChar.idx]
                        }
                        print(str(colNames))
                        return(colNames)
                })

		output$phRowsText <- renderText({
                        if(is.null(gVars$inputPh())){
				nRow <- "NA"
			}else{
				nRow <- nrow(gVars$inputPh())
			}

			return(paste0("Samples: ", nRow))
		})

		output$phColsText <- renderText({
                        if(is.null(gVars$inputPh())){
				nCol <- "NA"
			}else{
				nCol <- ncol(gVars$inputPh())
			}

			return(paste0("Variables: ", nCol))
		})

		gVars$phColChoices <- reactive({
			if(is.null(gVars$phLoaded))
			return(c("NA"))

			choicesVec <- seq(1,ncol(gVars$inputPh()))
			choicesNames <- paste0("Variable ", choicesVec)
			names(choicesVec) <- choicesNames
			return(choicesVec)
		})

		output$phenoDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$inputPh()), "No phenotype file!")
			)

                        phTable <- gVars$inputPh()
                        colnames(phTable) <- paste0(colnames(phTable), " [", c(1:ncol(phTable)), "]")
			DT::datatable(phTable, filter="none",
				options = list(
					ch = list(regex=TRUE, caseInsensitive=FALSE),
					scrollX=TRUE,
					pageLength=2,
					lengthMenu=c(1,2,3),
					ordering=F
				)
			)
		},server=TRUE)

                output$phenoTypesRH <- rhandsontable::renderRHandsontable({
			shiny::validate(
				need(!is.null(gVars$inputPh()), "No phenotype file!")
			)

                        phTable <- gVars$inputPh()
                        colNames <- paste0(colnames(phTable), " [", c(1:ncol(phTable)), "]")
                        colTypesList <- gVars$phColTypes()
                        print(str(colTypesList))
                        colClass <- unlist(lapply(phTable, class))
                        print(str(colClass))
                        colTypesDF <- data.frame(Variable=colNames, Type="-", Class=colClass, t(head(phTable)), stringsAsFactors=FALSE)
                        colnames(colTypesDF)[c(4:ncol(colTypesDF))] <- paste0("Sample", c(1:(ncol(colTypesDF)-3)))
                        if(!is.null(colTypesList[["c"]])){
                                print("In colTypes C...")
                                cIdx <- which(colnames(phTable) %in% colTypesList[["c"]])
                                print(cIdx)
                                if(length(cIdx>0)){
                                       colTypesDF[cIdx,"Type"] <- "factor"
                                }
                        }
                        if(!is.null(colTypesList[["n"]])){
                                print("In colTypes N...")
                                nIdx <- which(colnames(phTable) %in% colTypesList[["n"]])
                                print(nIdx)
                                if(length(nIdx>0)){
                                       colTypesDF[nIdx,"Type"] <- "vector"
                                }
                        }
                        print("Making rhandsontable...")
			rhandsontable::rhandsontable(colTypesDF, rowHeaders=NULL, readOnly=TRUE, contextMenu=FALSE) %>%
                                hot_col("Type", readOnly=FALSE, type="dropdown", source=c("factor","vector"),
                                        renderer = "function (instance, td, row, col, prop, value, cellProperties) {
                                                Handsontable.renderers.NumericRenderer.apply(this, arguments);
                                                if (value == 'factor') {
                                                            td.style.background = 'lightblue';
                                                } else if (value == 'vector') {
                                                            td.style.background = 'lightgreen';
                                                }
                                        }"
                                )
		})

                observeEvent(input$upload_pheno_submit, {
			shiny::validate(
				need(!is.null(gVars$inputPh()), "No Phenotype File Provided!")
			)

                        phTable <- gVars$inputPh()
			arrType <- input$arrType
			dyeColID <- NULL
			sampleColID <- as.integer(input$sampleIDCol)
			fileNameColID <- as.integer(input$fileNameCol)
			dyeColName <- NULL
			sampleColName <- colnames(phTable)[as.integer(input$sampleIDCol)]
			fileNameColName <- colnames(phTable)[as.integer(input$fileNameCol)]

                        #sampleIDs <- phTable[,sampleColID]
                        sampleIDs <- phTable[,sampleColName]
                        if(any(duplicated(sampleIDs))){
                                shinyjs::info(paste0("Sample IDs are not unique!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has unique values!"))
                                return(NULL)
                        }else if(any(sampleIDs=="") || any(sampleIDs==" ") || any(is.na(sampleIDs))){
                                shinyjs::info(paste0("Sample IDs contain BLANK and/or NA values!\n\nPlease check the phenotype data and ensure that the selected Sample ID column has complete information!"))
                                return(NULL)
                        }

			if(arrType=="ag_exp2"){
				dyeColID <- as.integer(input$dyeCol)
				dyeColName <- colnames(phTable)[as.integer(input$dyeCol)]
				phTable <- phTable[order(phTable[,dyeColID], phTable[,fileNameColID]),]
			}

                        #Make sample ID column character
                        if(is.integer(phTable[,sampleColName])){
                                print("Updating sample ID column from integer to character!")
                                phTable[,sampleColName] <- as.character(phTable[,sampleColName])
                        }

                        #Create factorized phTable
                        print("Create factorized phTable")
                        print(is.null(input$phenoTypesRH))
                        print(str(input$phenoTypesRH))
                        phRH <- rhandsontable::hot_to_r(input$phenoTypesRH)
                        print(str(phRH))
                        factorIdx <- which(phRH$Type=="factor")
                        print(factorIdx)
                        factorCols <- NULL
                        if(length(factorIdx)>0){
                                factorCols <- colnames(phTable)[factorIdx]
                                phFactor <- factorize_cols(phTable, factorIdx)
                        }else{
                                phFactor <- phTable
                        }
                        print(str(phFactor))
                        print(str(phTable))

                        gVars$phTable <- phTable
                        gVars$phFactor <- phFactor
                        gVars$factorCols <- factorCols
			gVars$dyeColID <- dyeColID
			gVars$fileNameColID <- fileNameColID
			gVars$sampleColID <- sampleColID
			gVars$dyeColName <- dyeColName
			gVars$fileNameColName <- fileNameColName
			gVars$sampleColName <- sampleColName
                        gVars$totalSamples <- nrow(phTable)
                        gVars$filteredSamples <- nrow(phTable)
                        gVars$removedSamples <- 0
			gVars$removedSamplesInfo <- NULL

			shinyBS::toggleModal(session, "importPhenoModal", toggle="close")
                        shinyBS::updateButton(session, "import_pheno_submit", style="success", icon=icon("check-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", open="LOAD RAW DATA", style=list("LOAD PHENOTYPE DATA"="success", "LOAD RAW DATA"="warning"))
                })

		gVars$pcChoices <- reactive({
                        phTable <- gVars$phTable
			if (is.null(phTable))
                        return(c("NA"))

                        colNames <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
                        if(length(colNames)==0){
                                shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
                                return(c("NA"))
                        }
			return(colNames)
		})

                gVars$pcChoicesSV <- reactive({
                        phTable <- gVars$phTable
			if (is.null(phTable))
                        return(c("NA"))

                        if(!is.null(gVars$svaSV)){
                                svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                svaSVc <- gVars$svaSVc
                                phTable <- cbind(phTable, svaSV, svaSVc)
                        }

                        colNames <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
                        if(length(colNames)==0){
                                shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
                                return(c("NA"))
                        }
			return(colNames)
		})

		gVars$condChoices <- reactive({
			if(is.null(gVars$phTable))
			return(c("NA"))

			varI <- input$varILimma
			if (is.null(varI) || is.na(varI))
			return(c("NA"))

                        phTable <- gVars$phTable
                        conds <- levels(factor(phTable[,varI]))
			return(conds)
		})

                observeEvent(input$filterPh, {
                        if(is.null(gVars$phTable) || is.null(input$filtered_rows_selected))
                        return(NULL)

			dyeColID <- gVars$dyeColID
			fileNameColID <- gVars$fileNameColID
                        sampleColID <- gVars$sampleColID
			dyeColName <- gVars$dyeColName
			fileNameColName <- gVars$fileNameColName
                        sampleColName <- gVars$sampleColName
                        phTable <- gVars$phTable
                        phRowsSel <- as.integer(input$filtered_rows_selected)
			print("Sel Rows:")
			print(phRowsSel)
                        phRows2Remove <- NULL
                        removedSampleIDs <- NULL
			if (length(phRowsSel)>0){
				phRows2Remove <- phRowsSel
				print("Rows to Remove:")
				print(phRows2Remove)
                                removedSampleIDs <- phTable[phRows2Remove,sampleColName]
				removedSamplesInfo <- phTable[phRows2Remove,c(fileNameColName, dyeColName)]
				phTable <- phTable[-phRows2Remove,]
                                cat("Removed Sample IDs :", paste0(removedSampleIDs, collapse=","), "\n")
			}

                        #Update expression objects
                        if(!is.null(gVars$norm.data)){
                                print("Subsetting gVars$norm.data")
                                rmIdx <- which(colnames(gVars$norm.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$norm.data matrix!\n")
                                        gVars$norm.data <- gVars$norm.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$rgList)){
                                print("Subsetting gVars$rgList")
                                removedArrays <- phTable[phRows2Remove,fileNameColName]
                                print("removedArrays: ")
                                print(removedArrays)
                                removedArrays <- gsub("\\.[a-zA-Z]+$", "", removedArrays)
                                print("removedArrays: ")
                                print(removedArrays)
                                print("table(removedArrays): ")
                                print(table(removedArrays))
                                array2Remove <- names(which(table(removedArrays)>1))
                                if(length(array2Remove)>0){
                                        gVars$rgList <- gVars$rgList[,-array2Remove]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$rgList.norm)){
                                print("Subsetting gVars$rgList.norm")
                                rIdx <- which(colnames(gVars$rgList.norm$R) %in% removedSampleIDs)
                                rmCheck <- 0
                                if(length(rIdx)>0){
                                        cat("Removing index ", rIdx, " from R matrix!\n")
                                        gVars$rgList.norm$R <- gVars$rgList.norm$R[,-rIdx]
                                        rmCheck <- 1
                                }
                                gIdx <- which(colnames(gVars$rgList.norm$G) %in% removedSampleIDs)
                                if(length(gIdx)>0){
                                        cat("Removing index ", gIdx, " from G matrix!\n")
                                        gVars$rgList.norm$G <- gVars$rgList.norm$G[,-gIdx]
                                        rmCheck <- 1
                                }
                                if(rmCheck==0){
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$RGset)){
                                print("Subsetting gVars$RGset")
                                removedArrays <- phTable[phRows2Remove,fileNameColName]
                                rmIdx <- which(sampleNames(gVars$RGset) %in% removedArrays)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$RGset matrix!\n")
                                        gVars$RGset <- gVars$RGset[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$Mset)){
                                print("Subsetting gVars$Mset")
                                removedArrays <- phTable[phRows2Remove,fileNameColName]
                                rmIdx <- which(sampleNames(gVars$Mset) %in% removedArrays)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$Mset matrix!\n")
                                        gVars$Mset <- gVars$Mset[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$expr.data)){
                                print("Subsetting gVars$expr.data")
                                rmIdx <- which(colnames(gVars$expr.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$expr.data matrix!\n")
                                        gVars$expr.data <- gVars$expr.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$nc.data)){
                                print("Subsetting gVars$nc.data")
                                rmIdx <- which(colnames(gVars$nc.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$nc.data matrix!\n")
                                        gVars$nc.data <- gVars$nc.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$c.data)){
                                print("Subsetting gVars$c.data")
                                rmIdx <- which(colnames(gVars$c.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$c.data matrix!\n")
                                        gVars$c.data <- gVars$c.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$comb.data)){
                                print("Subsetting gVars$comb.data")
                                rmIdx <- which(colnames(gVars$comb.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$comb.data matrix!\n")
                                        gVars$comb.data <- gVars$comb.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$agg.data)){
                                print("Subsetting gVars$agg.data")
                                rmIdx <- which(colnames(gVars$agg.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$agg.data matrix!\n")
                                        gVars$agg.data <- gVars$agg.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }
                        if(!is.null(gVars$comb.sva.data)){
                                print("Subsetting gVars$comb.sva.data")
                                rmIdx <- which(colnames(gVars$comb.sva.data) %in% removedSampleIDs)
                                if(length(rmIdx)>0){
                                        cat("Removing index ", rmIdx, " from gVars$comb.sva.data matrix!\n")
                                        gVars$comb.sva.data <- gVars$comb.sva.data[,-rmIdx]
                                }else{
                                        print("Nothing was removed!")
                                }
                        }

                        #Remove columns with single level data
                        nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
                        nrlevels.singular <- which(nrlevels==1)

                        remInfo <- 0
                        remStr <- ""
                        if(length(nrlevels.singular)>0){
                                remInfo <- 1
                                remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
                                if(length(nrlevels.singular)==ncol(phTable)){
                                        remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
                                        shinyjs::info(remStr)
                                        return(NULL)
                                }
                                col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
                                phTable <- phTable[,-col2rem, drop=F]
                        }

                        #Inform user with the columns removed from the data frame
                        if(remInfo==1){
                                shinyjs::info(remStr)
                        }

                        factorCols <- gVars$factorCols
                        if(is.null(factorCols)){
                                phFactor <- phTable
                        }else{
                                tmpIdx <- which(factorCols %in% colnames(phTable))
                                if(length(tmpIdx)>0){
                                        factorCols <- factorCols[tmpIdx]
                                        factorIdx <- which(colnames(phTable) %in% factorCols)
                                        phFactor <- factorize_cols(phTable, factorIdx)
                                }else{
                                        phFactor <- phTable
                                }
                        }
                        gVars$phTable <- phTable
                        gVars$phFactor <- phFactor
                        gVars$filteredSamples <- nrow(phTable)
                        gVars$removedSamples <- gVars$totalSamples - nrow(phTable)
			gVars$removedSamplesInfo <- rbind(gVars$removedSamplesInfo, removedSamplesInfo)
			print("Removed Samples:")
			print(gVars$removedSamplesInfo)
                })

		output$filtered <- DT::renderDataTable({
                        if(is.null(gVars$phTable))
                        return(NULL)

                        phTable <- gVars$phTable
			DT::datatable(phTable, filter="none",
				options = list(
					search = list(regex=TRUE, caseInsensitive=FALSE),
					scrollX=TRUE,
					ordering=F
				)
			)
		},server=TRUE)

                observeEvent(input$upload_raw_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#Start loading screen
			shinyjs::html(id="loadingText", "UPLOADING RAW DATA")
			shinyjs::show(id="loading-content")

                        on.exit({
				progress$close()
				shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
			})

                        if(is.null(gVars$phTable) || is.null(gVars$rawDir)){
                                return(NULL)
                        }

			#Clear filter result reactive variables
                        gVars$filt.data <- NULL
                        gVars$keepMeth <- NULL

			#Clear normalization result reactive variables
                        gVars$expr.data <- NULL
                        gVars$norm.data <- NULL
                        gVars$rgList.norm <- NULL
                        gVars$Mset <- NULL

                        #Clear SVA result reative variables
                        gVars$varISva <- NULL
                        gVars$svaSV <- NULL
                        gVars$svaSVc <- NULL
                        gVars$svaStep <- NULL

                        #Clear batch correction result reactive variables
                        gVars$comb.data <- NULL
			gVars$agg.data <- NULL
                        gVars$conditions <- NULL
                        gVars$conditionsBoxplot <- NULL
                        gVars$correctionLvl <- NULL

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$QC_passed <- FALSE
                        gVars$filtered <- FALSE
                        gVars$normalized <- FALSE
                        gVars$corrected <- FALSE
                        gVars$loadedRaw <- FALSE

                        progress$set(message="Raw Data:", value=0)

                        phTable <- gVars$phTable
                        phFactor <- gVars$phFactor
			removedSamplesInfo <- gVars$removedSamplesInfo
			fileNameColID <- gVars$fileNameColID
                        sampleColID <- gVars$sampleColID
			dyeColID <- gVars$dyeColID
			fileNameColName <- gVars$fileNameColName
                        sampleColName <- gVars$sampleColName
			dyeColName <- gVars$dyeColName
			arrType <- input$arrType
			affCDF <- input$affCDF

			####New upload raw data feature

			rawDir <- gVars$rawDir
			filePaths <- gVars$filePaths
			print("str(filePaths)")
			print(str(filePaths))
			print(paste0("length: ", length(filePaths)))

			idxLoaded <- gVars$idxLoaded
			totalSamples <-  length(as.vector(phTable[,fileNameColName]))
			if(is.null(idxLoaded) || length(idxLoaded)==0){
				errStr <- "No raw data filei were uploaded!!!"
				shinyjs::info(errStr)
				return(NULL)
			}else if(totalSamples>length(idxLoaded)){
				msgStr <- paste0((totalSamples-length(idxLoaded)), " out of", totalSamples, " samples were not uploaded!!!\n\nMissing samples will be removed!!!")
				shinyjs::info(msgStr)

				removedSampleIDs <- phTable[-idxLoaded,sampleColName]
				cat("Removed Sample IDs :", paste0(removedSampleIDs, collapse=","), "\n")
				phTable <-  phTable[idxLoaded,,drop=FALSE]

				#Remove columns with single level data
				nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
				nrlevels.singular <- which(nrlevels==1)

				remInfo <- 0
				remStr <- ""
				if(length(nrlevels.singular)>0){
					remInfo <- 1
					remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
					if(length(nrlevels.singular)==ncol(phTable)){
						remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
						shinyjs::info(remStr)
						return(NULL)
					}
					col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
					phTable <- phTable[,-col2rem, drop=F]
				}

				#Inform user with the columns removed from the data frame
				if(remInfo==1){
					shinyjs::info(remStr)
				}

				factorCols <- gVars$factorCols
				if(is.null(factorCols)){
					phFactor <- phTable
				}else{
					tmpIdx <- which(factorCols %in% colnames(phTable))
					if(length(tmpIdx)>0){
						factorCols <- factorCols[tmpIdx]
						factorIdx <- which(colnames(phTable) %in% factorCols)
						phFactor <- factorize_cols(phTable, factorIdx)
					}else{
						phFactor <- phTable
					}
				}
				gVars$phTable <- phTable
				gVars$phFactor <- phFactor

				gVars$filteredSamples <- nrow(phTable)
				gVars$removedSamples <- gVars$totalSamples - nrow(phTable)
			}

			####New upload raw data feature

                        updateProgress(detail="Reading Raw Data...", value=1/3)
			fileNames <- unique(phTable[,fileNameColName])
                        fileNameCount <- length(fileNames)
                        if(arrType=="il_methyl"){
                                fileNamesTmp <- as.vector(sapply(fileNames, function(x){g<-paste0(x, "_Grn.idat"); r<-paste0(x, "_Red.idat"); return(c(g,r))}))
                                fileNameCount <- length(fileNamesTmp)

				filePaths <- filePaths[fileNamesTmp]
                        }else{
				filePaths <- filePaths[fileNames]
                        }

			if(arrType=="ag_exp2"){
                                ### FIXING TO KEEP SAMPLES IN PD [START] ###
                                smp2keepG <- NULL
                                smp2keepR <- NULL
                                sampleNamesGreen <- NULL
                                sampleNamesRed <- NULL

				#samplesInfo <- phTable[, c(fileNameColID, dyeColID, sampleColID)]
				samplesInfo <- phTable[, c(fileNameColName, dyeColName, sampleColName)]
                                print("head(samplesInfo):")
                                print(head(samplesInfo))
                                samplesInfo[,1] <- gsub("\\.[a-zA-Z]*$", "", samplesInfo[,1])
                                smp2keepID <- grep("Cy3", samplesInfo[,2], ignore.case=T)
                                if(length(smp2keepID)>0){
                                        smp2keepG <- samplesInfo[smp2keepID,1]
                                        sampleNamesGreen <- setNames(samplesInfo[smp2keepID,3], samplesInfo[smp2keepID,1])
                                        print("Green Samples to Keep")
                                        print(smp2keepG)
                                }
                                smp2keepID <- grep("Cy5", samplesInfo[,2], ignore.case=T)
                                if(length(smp2keepID)>0){
                                        smp2keepR <- samplesInfo[smp2keepID,1]
                                        sampleNamesRed <- setNames(samplesInfo[smp2keepID,3], samplesInfo[smp2keepID,1])
                                        print("Red Samples to Keep")
                                        print(smp2keepR)
                                }
                                #print("Samples Info:")
                                #print(head(samplesInfo))
                                #print(smp2keepG)
                                #print(smp2keepR)
                                ### FIXING TO KEEP SAMPLES IN PD [STOP] ###

                                sourceOpt <- input$selAgSource
                                cat("Agilent Source: ", sourceOpt, "\n")

				#Check file source match
				columns <- switch(sourceOpt,
					agilent.mean = list(
						G = "gMeanSignal",
						Gb = "gBGMedianSignal",
						R = "rMeanSignal",
						Rb = "rBGMedianSignal"
					),
					agilent = ,
					agilent.median = list(
						G = "gMedianSignal",
						Gb = "gBGMedianSignal",
						R = "rMedianSignal",
						Rb = "rBGMedianSignal"
					),
					genepix = ,
					genepix.mean = list(
						R = "F635 Mean",
						G = "F532 Mean",
						Rb = "B635 Median",
						Gb = "B532 Median"
					),
					genepix.median = list(
						R = "F635 Median",
						G = "F532 Median",
						Rb = "B635 Median",
						Gb = "B532 Median"
					),
					NULL
				)

				file <- file.path(rawDir, fileNames[1])
				con <- file(file, "r")
				Found <- FALSE
				i <- 0
				repeat {
					i <- i + 1
					txt <- readLines(con, n = 1)
					if (!length(txt)){
						Found <- FALSE
						break
					}
					Found <- TRUE
					for (a in columns) Found <- Found && length(grep(a, txt))
					if (Found)
						break
				}
				close(con)
				if(!Found){
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                        shinyjs::info(
                                                paste0("Source type mismatch with raw data files!\nPlease check and select the correct source!")
                                        )
                                        return(NULL)
				}

				#Read raw data
				rgList <- limma::read.maimages(files=fileNames, source=sourceOpt, path=rawDir, verbose=TRUE, skipNul=TRUE)
				updateProgress(detail="Preprocessing & Filtering...", value=2/3)
				gCh <- rgList$G
				rCh <- rgList$R

                                if(is.null(smp2keepG)){
                                        gCh <- NULL
                                }else{
                                        print("Samples Read Green Before Filter")
					colID <- which(colnames(gCh) %in% smp2keepG)
					gCh <- gCh[,c(colID), drop=FALSE]
                                        print("Samples Read Green After Filter")
                                        print(colnames(gCh))
                                        print("Samples Green Colnames After Update")
                                        colnames(gCh) <- sampleNamesGreen[colnames(gCh)]
                                        print(colnames(gCh))
					print(str(gCh))
                                }
                                if(is.null(smp2keepR)){
                                        rCh <- NULL
                                }else{
                                        print("Samples Read Red Before Filter")
                                        print(colnames(rCh))
					colID <- which(colnames(rCh) %in% smp2keepR)
					rCh <- rCh[,c(colID), drop=FALSE]
                                        print("Samples Read Red After Filter")
                                        print(colnames(rCh))
                                        print("Samples Red Colnames After Update")
                                        colnames(rCh) <- sampleNamesRed[colnames(rCh)]
                                        print(colnames(rCh))
					print(str(rCh))
                                }

				#data <- cbind(rgList$G, rgList$R)
				data <- cbind(gCh, rCh)
				rownames(data) <- rgList$genes$ProbeName
                                print("samplesInfo[,3] : ")
                                print(samplesInfo[,3])
                                data <- data[,samplesInfo[,3]] ## Reorder expression data after filtering
                                print("EXPR Colnames")
                                print(colnames(data))
                                ncIdx <- which(rgList$genes$ControlType==0)
                                cIdx <- which(rgList$genes$ControlType==-1)
				nc.data <- data[ncIdx,]
				c.data  <- data[cIdx,]
				updateProgress(detail="Completed!", value=2/3)

				gVars$rgList <- rgList
				gVars$nc.data <- nc.data
				gVars$c.data <- c.data
			}else if(arrType=="ag_exp1"){
				eListRaw <- limma::read.maimages(files=fileNames, source="agilent.median", path=rawDir, green.only=T, other.columns=c("gIsGeneDetected"),verbose=TRUE, skipNul=TRUE)
                                data <- eListRaw$E

				rownames(data) <- eListRaw$genes$ProbeName

				samplesInfo <- phTable[, c(fileNameColName, sampleColName)]
                                samplesInfo[,1] <- gsub("\\.[a-zA-Z]*$", "", samplesInfo[,1])
                                sampleNames <- setNames(samplesInfo[,2], samplesInfo[,1])
                                colnames(data) <- sampleNames[colnames(data)]
                                ncIdx <- which(eListRaw$genes$ControlType==0)
                                cIdx <- which(eListRaw$genes$ControlType==-1)
                                nc.data <- data[ncIdx,]
                                c.data <- data[cIdx,]

                                map <- data.frame(eListRaw$genes$ProbeName[ncIdx], eListRaw$genes$SystematicName[ncIdx], stringsAsFactors=FALSE)
                                dups <- which(duplicated(map[,1]))
                                if(length(dups)>0){
                                        map <- map[-dups,]
                                }

				updateProgress(detail="Preprocessing & Filtering...", value=2/3)
				gVars$nc.data <- nc.data
				gVars$c.data <- c.data
                                gVars$map <- map
                                gVars$eListRaw <- eListRaw
			}else if(arrType=="af_exp"){
                                tmpFile <- file.path(rawDir, fileNames[1])
                                celHeader <- affyio::read.celfile.header(tmpFile)
                                print("str(celHeader)")
                                print(str(celHeader))
                                cdfName <- gsub("-|_", "", celHeader$cdfName)
                                print("cdfName:")
                                print(cdfName)
                                allInstalledPkgs <- rownames(installed.packages())
                                cdfNameBioc <- paste0(tolower(cdfName), "cdf")
                                print("Bioconductor cdfName:")
                                print(cdfNameBioc)
                                if(!any(grepl(cdfName, affCDF, ignore.case=T))){
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                        shinyjs::info(
                                                paste0("CDF annotation mismatch with the CEL files!\n\nNeed CDF file corresponding to '", cdfName, "'\n\nUser provided CDF '", affCDF, "'")
                                        )
                                        return(NULL)
                                }

                                rownames(phTable) <- phTable[,sampleColName]
                                rownames(phFactor) <- phTable[,sampleColName]
                                pheno <- new("AnnotatedDataFrame", data=phFactor)

				eset <- affy::justRMA(filenames=fileNames, celfile.path=rawDir, phenoData=pheno, sampleNames=phTable[,sampleColName], normalize=TRUE, background=TRUE, cdfname=affCDF)
				exprs <- affy::exprs(eset)
				gVars$norm.data <- exprs
				gVars$expr.data <- exprs
				gVars$rma.affy.data <- exprs
                                gVars$normalized <- TRUE

                                tmpChoices <- rownames(gVars$expr.data)
                                updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                                gVars$expGenes <- NULL
                                shiny::updateTabsetPanel(session, "display", selected="normTab")
			}else if(arrType=="il_methyl"){
                                pdFactor <- phFactor
                                pdFactor$Basename <- factor(file.path(rawDir, fileNames))
                                RGset <- minfi::read.metharray.exp(targets=pdFactor)

				RGarray <- NULL
                                if(minfi:::.isEPIC(RGset)){
                                        RGarray <- "EPIC"
                                }else if(minfi:::.is450k(RGset)){
                                        RGarray <- "450k"
                                }

                                detP <- minfi::detectionP(RGset)

                                gVars$RGarray <- RGarray
                                gVars$RGset <- RGset
                                gVars$detP <- detP
			}
                        shinyBS::updateCollapse(session, "bsSidebar", open="QUALITY CONTROL", style=list("LOAD RAW DATA"="success", "QUALITY CONTROL"="warning"))
                        gVars$loadedRaw <- TRUE

                        #Update buttons and sidebar style for following steps
                        shinyBS::updateButton(session, "launch_sva_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_combat_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", style = list("PROBE FILTERING"="danger", "NORMALIZATION"="danger", "BATCH CORRECTION"="danger", "DIFFERENTIAL ANALYSIS"="danger", "ANNOTATION"="danger"))
		})

                observeEvent(input$qc_skip_submit, {
                        if(!gVars$QC_passed){
                                if(input$arrType=="af_exp"){
                                        shinyBS::updateCollapse(session, "bsSidebar", open="BATCH CORRECTION", style=list("QUALITY CONTROL"="success", "BATCH CORRECTION"="warning"))
                                }else if(input$arrType=="ag_exp2"){
                                        shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
                                }else if(input$arrType=="ag_exp1"){
                                        shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
                                }else if(input$arrType=="il_methyl"){
                                        shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
                                }
                                gVars$QC_passed <- TRUE
                        }
                })

                observeEvent(input$qc_methyl_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "CREATING QC REPORT")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        RGset <- gVars$RGset
                        if(is.null(RGset)){
                                shinyjs::info("Raw data not available for QC!")
                                return(NULL)
                        }

                        progress$set(message="QC:", value=0)

                        updateProgress(detail="summarizing", value=1/2)
                        RGsetSummary <- shinyMethyl::shinySummarize(RGset)
                        message("Creating UI...")
                        uiShinyMethyl <- shinyMethyl:::ui.shinyMethyl(RGsetSummary)
                        message("Creating Server...")
                        serverShinyMethyl <- shinyMethyl:::server.shinyMethyl(RGsetSummary)

                        updateProgress(detail="loading", value=2/2)

                        #tempRdata <- file.path(tempdir(), "RGsetSummary.RData")
			#tempRdata <- tempfile(pattern="RGsetSummary_", fileext=".RData")
                        tempRdata <- file.path(".", "RGsetSummary.RData")
                        save(list=c("uiShinyMethyl", "serverShinyMethyl"), file=tempRdata)
			Sys.setenv(SHINY_METHYL_TEMP=tempRdata)
			message('Sys.getenv("SHINY_METHYL_TEMP")')
			message(Sys.getenv("SHINY_METHYL_TEMP"))

                        message("Created shinyMethyl report!")
			shinyBS::toggleModal(session, "shinyMethylModal", toggle="open")

                        if(!gVars$QC_passed){
                                shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
                        }

                        gVars$QC_passed <- TRUE
                })

		observeEvent(input$filt_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "FILTERING DATA")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        #Clear normalization result reactive variables
                        gVars$expr.data <- gVars$rma.affy.data
                        gVars$norm.data <- NULL
                        gVars$rgList.norm <- NULL
                        gVars$Mset <- NULL

                        #Clear SVA result reative variables
                        gVars$varISva <- NULL
                        gVars$svaSV <- NULL
                        gVars$svaSVc <- NULL
                        gVars$svaStep <- NULL

                        #Clear batch correction result reactive variables
                        gVars$comb.data <- NULL
			gVars$agg.data <- NULL
                        gVars$conditions <- NULL
                        gVars$conditionsBoxplot <- NULL
                        gVars$correctionLvl <- NULL

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$filtered <- FALSE
                        gVars$normalized <- FALSE
                        gVars$corrected <- FALSE

                        arrType <- input$arrType
                        RGset <- gVars$RGset
                        detP <- gVars$detP
			if(arrType=="il_methyl"){
                                if(is.null(RGset))
                                return(NULL)
                        }else{
                                if(is.null(gVars$nc.data) || is.null(gVars$c.data))
                                return(NULL)
                        }

			nc.data <- gVars$nc.data
			c.data <- gVars$c.data
			qfilt <- as.numeric(input$filtDist)
			perc <- as.numeric(input$perSamples)

                        progress$set(message="Filtering:", value=0)

                        updateProgress(detail="...", value=1/2)
                        if(arrType=="il_methyl"){
                                percSample <- perc
                                pvThr <- as.numeric(input$detectPV)
                                sampleCount <- round(percSample*ncol(RGset)/100)
                                keep <- rowSums(detP<pvThr)>=sampleCount
                                print("~~~~~~~~~~~~~~~~~~~~RGset DIM:")
                                print(dim(RGset))
                                print("~~~~~~~~~~~~~~~~~~~~detP DIM:")
                                print(dim(detP))
                                print("~~~~~~~~~~~~~~~~~~~~detP HEAD:")
                                print(head(detP))
                                print("Sample Count:")
                                print(sampleCount)
                                print("PV Thr:")
                                print(pvThr)
                                print("Keep:")
                                print(table(keep))
                                gVars$keepMeth <- keep

                                filt.data <- RGset[keep,]

				MsetRaw <- minfi::preprocessRaw(RGset)
                                print("~~~~~~~~~~~~~~~~~~~~MsetRaw DIM:")
                                print(dim(MsetRaw))
                                print("~~~~~~~~~~~~~~~~~~~~MsetRawSummary:")
                                print(summary(minfi::getM(MsetRaw)[,1]))

                                MsetRawFilt <- MsetRaw[keep,]
                                print("~~~~~~~~~~~~~~~~~~~~MsetRaw Filt DIM:")
                                print(dim(MsetRawFilt))

                                #Map to genomic locations
                                GMsetFilt <- minfi::mapToGenome(MsetRawFilt)

				#Convert to M or Beta
                                valType <- "M"
                                methFilt.data <- minfi::getMethSignal(GMsetFilt, valType)
				gVars$methFilt.data <- methFilt.data
                        }else{
                                filt.data <- filt.by.neg.probes(nc.data, c.data, qdist=qfilt, perc=perc, verbose=T)
                        }
                        updateProgress(detail="Completed!", value=2/2)
                        shinyBS::updateCollapse(session, "bsSidebar", open="NORMALIZATION", style=list("PROBE FILTERING"="success", "NORMALIZATION"="warning"))
                        gVars$filt.data <- filt.data

                        gVars$filtered <- TRUE

                        #Update buttons and sidebar style for following steps
                        shinyBS::updateButton(session, "launch_sva_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_combat_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", style = list("BATCH CORRECTION"="danger", "DIFFERENTIAL ANALYSIS"="danger", "ANNOTATION"="danger"))
		})

		output$percProbesText <- renderText({
			if(is.null(gVars$filt.data))
			return("Probes Removed: NA")

			#filt.data <- gVars$filt.data()
			filt.data <- gVars$filt.data
			nc.data <- gVars$nc.data
                        arrType <- input$arrType
                        RGset <- gVars$RGset
                        detP <- gVars$detP

                        if(arrType=="il_methyl"){
                                keep <- gVars$keepMeth
                                countProbes <- sum(!keep)
                                percProbes <- round((countProbes/length(keep))*100, digits=2)
                        }else{
                                countProbes <- nrow(nc.data)-nrow(filt.data)
                                percProbes <- round((countProbes/nrow(nc.data))*100, digits=2)
                        }

			percProbesText <- paste0("Probes Removed: ", countProbes, " (", percProbes, ")")
			return(percProbesText)
		})

		output$numProbesText <- renderText({
			if(is.null(gVars$filt.data))
			return("Probes Remaining: NA")

                        arrType <- input$arrType
                        if(arrType=="il_methyl"){
                                keep <- gVars$keepMeth
                                countProbes <- sum(keep)
                        }else{
                            filt.data <- gVars$filt.data
                            countProbes <- nrow(filt.data)
                        }

			numProbesText <- paste0("Probes Remaining: ", countProbes)
			return(numProbesText)
		})

                observeEvent(input$norm_submit,{
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "NORMALIZING DATA")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

			if(is.null(gVars$filt.data)){
                                return(NULL)
                        }

                        #Clear SVA result reative variables
                        gVars$varISva <- NULL
                        gVars$svaSV <- NULL
                        gVars$svaSVc <- NULL
                        gVars$svaStep <- NULL

                        #Clear batch correction result reactive variables
                        gVars$comb.data <- NULL
			gVars$agg.data <- NULL
                        gVars$conditions <- NULL
                        gVars$conditionsBoxplot <- NULL
                        gVars$correctionLvl <- NULL

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$normalized <- FALSE
                        gVars$corrected <- FALSE

                        progress$set(message="Normalization:", value=0)

			#filt.data <- gVars$filt.data()
			filt.data <- gVars$filt.data
                        method <- input$normMethod
                        method2 <- input$normMethod2
			phTable <- gVars$phTable
			dyeColID <- gVars$dyeColID
                        sampleColID <- gVars$sampleColID
			dyeColName <- gVars$dyeColName
                        sampleColName <- gVars$sampleColName
                        arrType <- input$arrType

                        updateProgress(detail="Processing...", value=1/2)
                        Mset <- NULL
			if(arrType=="il_methyl"){
                                print("In if il_methyl NORM!")
                                #set.seed(1) #Block randomness. Set seed

                                RGset <- gVars$RGset
                                detP <- gVars$detP
                                print(apply(detP, 2, summary))
                                percSample <- as.numeric(input$perSamples)
                                pvThr <- as.numeric(input$detectPV)
                                sampleCount <- round(percSample*ncol(RGset)/100)
                                keep <- rowSums(detP<pvThr)>=sampleCount
                                print("~~~~~~~~~~~~~~~~~~~~RGset DIM:")
                                print(dim(RGset))
                                print("Sample Count:")
                                print(sampleCount)
                                print("PV Thr:")
                                print(pvThr)
                                print("Keep:")
                                print(table(keep))

                                MsetRaw <- minfi::preprocessRaw(RGset)
                                print("~~~~~~~~~~~~~~~~~~~~MsetRaw DIM:")
                                print(dim(MsetRaw))
                                print("~~~~~~~~~~~~~~~~~~~~MsetRawSummary:")
                                print(summary(minfi::getM(MsetRaw)[,1]))

				#Perform normalization
                                if(method=="SWAN"){
                                        print("Running preprocessSWAN...")
                                        Mset <- minfi::preprocessSWAN(RGset)
                                }else if(method=="Raw"){
                                        print("Running preprocessRaw...")
                                        Mset <- minfi::preprocessRaw(RGset)
                                }else if(method=="Quantile"){
                                        print("Running preprocessQuantile...")
                                        Mset <- minfi::preprocessQuantile(RGset)
                                }else if(method=="Illumina"){
                                        print("Running preprocessIllumina...")
                                        Mset <- minfi::preprocessIllumina(RGset)
                                }else if(method=="FunNorm"){
                                        print("Running preprocessFunnorm...")
                                        Mset <- minfi::preprocessFunnorm(RGset)
                                }else if(method=="Noob"){
                                        print("Running preprocessNoob...")
                                        Mset <- minfi::preprocessNoob(RGset)
                                }
                                print("~~~~~~~~~~~~~~~~~~~~Mset Norm DIM:")
                                print(dim(Mset))
                                Mset <- Mset[keep,]
                                print("~~~~~~~~~~~~~~~~~~~~Mset Filt DIM:")
                                print(dim(Mset))
                                print("~~~~~~~~~~~~~~~~~~~~MsetSummary:")
                                print(summary(minfi::getM(Mset)[,1]))

                                #Map to genomic locations
                                GMset <- minfi::mapToGenome(Mset)

                                #Remove SNP probes
                                if(input$chkRmSNP){
                                        GMset <- minfi::dropLociWithSnps(GMset)
                                }

                                #Remove non specific probes [OPTIONAL]
                                remStr <- NULL
                                cgcrossFile <- "cgcross.txt"
                                nonspec <- as.vector(read.table(cgcrossFile, header=T)[,1])
                                nonSpecIdx <- which(is.element(minfi::featureNames(GMset), nonspec))
                                if(length(nonSpecIdx)){
                                        GMset <- GMset[-nonSpecIdx,]
                                        remStr <- paste0("Removed ", length(nonSpecIdx), " nonspecific features!\n")
                                }

                                #Make autosomal, remove sex chromosome specific probes
                                if(input$chkRmSexChr){
                                        annType <- paste0(minfi::annotation(RGset)[["array"]], "anno.", minfi::annotation(RGset)[["annotation"]])
                                        ann <- minfi::getAnnotation(annType)
                                        sexChrIdx <- which(ann$chr=="chrX" | ann$chr=="chrY")
                                        remCount <- 0
                                        if(length(sexChrIdx)>0){
                                                print("Checking sex chromosome probes...")
                                                print(length(sexChrIdx))
                                                autosomalProbes <- rownames(ann)[-sexChrIdx]
                                                autsomalIdx <- which(is.element(featureNames(GMset), autosomalProbes))
                                                print("Checking autosomal probes...")
                                                print(length(autsomalIdx))
                                                remCount <- length(minfi::featureNames(GMset)) - length(autsomalIdx)
                                                remStr <- paste0(remStr, "\nRemoved ", remCount, " sex chromsome specific features!\n")
                                                if(length(autsomalIdx)==0){
                                                        remStr <- paste0(remStr, "\nNo features remaining after removing sex chromosome features!\nABORTED!!!")
                                                        shinyjs::info(remStr)
                                                        return(NULL)
                                                }
                                                GMset <- GMset[autsomalIdx,]
                                        }
                                }

                                if(!is.null(remStr)){
                                        shinyjs::info(remStr)
                                }

                                #Convert to M or Beta
                                valType <- "M"
                                norm.data <- minfi::getMethSignal(GMset, valType)
				colnames(norm.data) <- phTable[,sampleColName] ## Update colnames to match sampleIDs
                                print("Exiting if il_methyl NORM!")
                        }else{
			        norm.data <- run.norm(filt.data=filt.data, method=method, method2=method2, plot=FALSE)
                        }

                        rgList.norm <- NULL
			if(arrType=="ag_exp2"){
                                dyeColChk <- which(dyeColName %in% colnames(phTable))
                                if(length(dyeColChk)>0){
                                        gSampleCount <- table(phTable[,dyeColName])[[1]]
                                        rSampleCount <- table(phTable[,dyeColName])[[2]]
                                        ttlSampleCount <- gSampleCount+rSampleCount
                                        rgList.norm <- new("RGList", list(G=norm.data[,c(1:gSampleCount)], R=norm.data[,c((gSampleCount+1):ttlSampleCount)]))
                                }
                        }
                        updateProgress(detail="Completed!", value=2/2)

                        gVars$expr.data <- norm.data
                        gVars$norm.data <- norm.data
                        gVars$rgList.norm <- rgList.norm
                        gVars$Mset <- Mset
                        print("Completed NORM!")

                        tmpChoices <- rownames(gVars$expr.data)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

                        shiny::updateTabsetPanel(session, "display", selected="normTab")
                        shinyBS::updateCollapse(session, "bsSidebar", open="BATCH CORRECTION", style=list("NORMALIZATION"="success", "BATCH CORRECTION"="warning"))
                        gVars$normalized <- TRUE

                        #Update buttons and sidebar style for following steps
                        shinyBS::updateButton(session, "launch_sva_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_combat_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", style=list("DIFFERENTIAL ANALYSIS"="danger", "ANNOTATION"="danger"))
                })

                observeEvent(input$normMethod, {
                        method <- input$normMethod
                        baChoices <- gVars$baChoices
                        if(method=="BA"){
                                tmpChoices <- gVars$baChoices
                        }else{
                                tmpChoices <- c("NA")
                        }
                        updateSelectInput(session, "normMethod2", choices=tmpChoices)
                })

		observeEvent(input$sva_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "SURROGATE VARIABLE ANALYSIS")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        if(is.null(gVars$rawDir) || is.null(gVars$norm.data))
                        return(NULL)

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$corrected <- FALSE

                        progress$set(message="Batch Correction:", value=0)

                        varI <- input$varISva
                        if(is.null(input$coVarSva)){
				coVar <- NULL
                        }else{
				coVar <- as.list(input$coVarSva)
                        }

			batch <- as.list(input$batchSva)
                        data <- gVars$norm.data
			phTable <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(phTable) <- phTable[,gVars$sampleColName]
                        rownames(phFactor) <- phTable[,gVars$sampleColName]

			arrType <- input$arrType
                        correctionLvl <- NULL

			batchCorVar <- list(var.int=varI, covariates=coVar, batches=batch)
			print(str(batchCorVar))

                        updateProgress(detail="Getting Surrogate Variables...", value=1/3)
                        batches.sva <- get.sva.batch.effects(comb.data=data, pd=phFactor, vars=batchCorVar, cmd.ui=F)
			print("str(batches.sva):")
			print(str(batches.sva))

                        updateProgress(detail="Filtering Confounded Variables...", value=2/3)
                        assoc.cutoff <- 0.05
                        sv.filt.logic <- apply(batches.sva$pd, 2, function(x) {val<-x[batchCorVar$var.int]; if(is.finite(val)){val<assoc.cutoff}else{FALSE}})
                        sv.filt.names <- names(sv.filt.logic[sv.filt.logic==F])
                        print("sv.filt.names:")
                        print(sv.filt.names)
                        if(length(sv.filt.names)>0){
				#Remove columns with single level data
				filt <- batches.sva$sv[,sv.filt.names, FALSE]
				nrlevels <- apply(filt, 2, function(x){length(levels(factor(x)))})
				nrlevels.singular <- which(nrlevels==1)

				remInfo <- 0
				remStr <- ""
				if(length(nrlevels.singular)>0){
					remInfo <- 1
					remStr <- paste0(remStr, "Following surrogate variables are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
					if(length(nrlevels.singular)==ncol(filt)){
						remStr <- paste0(remStr, "\n\nNo surrogate variables survived filtering!!! Please proceed with known variables.")
						shinyjs::info(remStr)
						return(NULL)
					}
					sv2rem <- which(sv.filt.names %in% names(nrlevels.singular))
					sv.filt.names <- sv.filt.names[,-sv2rem, drop=F]
				}

                                if(length(sv.filt.names)>1){
                                        #Check for confounded sva variables
                                        sva.assoc.mat <- assoc.var.int(batches.sva$sv[sv.filt.names], batches.sva$sv[sv.filt.names])
                                        colnames(sva.assoc.mat) <- rownames(sva.assoc.mat)
                                        sva.assoc.mat[lower.tri(sva.assoc.mat, diag=T)] <- NA
                                        sva.assoc.DF <- as.data.frame(as.table(sva.assoc.mat), stringsAsFactors=F)
                                        sva.assoc.DF <- sva.assoc.DF[-which(is.na(sva.assoc.DF$Freq)),]
                                        print("dim(sva.assoc.DF):")
                                        print(dim(sva.assoc.DF))
                                        print("head(sva.assoc.DF):")
                                        print(head(sva.assoc.DF))
                                        rowSel <- which(sva.assoc.DF$Freq<0.01)
                                        if(length(rowSel)>0){
                                                sva.assoc.DF.conf <- sva.assoc.DF[rowSel,]
                                                sva.rm.names <- unique(sva.assoc.DF.conf[,2])
                                                sv.filt.names <- sv.filt.names[-which(sv.filt.names %in% sva.rm.names)]
                                        }
                                }
                                #Add selected sva sv and perform combat
                                sva.sv.filt <- batches.sva$sv[,sv.filt.names]
                                sva.svc.filt <- batches.sva$svc[,sv.filt.names]
                                print(class(sva.sv.filt))
                                if(!is.data.frame(sva.sv.filt) && !is.matrix(sva.sv.filt)){
                                        sva.sv.filt <- data.frame(sva.sv.filt)
                                        sva.svc.filt <- data.frame(sva.svc.filt)
                                }
                                colnames(sva.sv.filt) <- paste("svaD",c(1:ncol(sva.sv.filt)),sep=".")
                                colnames(sva.svc.filt) <- paste("svaC",c(1:ncol(sva.svc.filt)),sep=".")
                                print(class(sva.sv.filt))
                                print(dim(sva.sv.filt))
                                print(head(sva.sv.filt))

                                svaSV <- sva.sv.filt
                                svaSVc <- sva.svc.filt
                                svaStep <- 1
                        }else{
                                shinyjs::info(paste0("All surrogate variables are confounded with variable of interest! Proceed to ComBat with known variables OR rerun SVA with different model!"))
                                svaStep <- 0
				svaSV <- NULL
				svaSVc <- NULL
                        }
                        updateProgress(detail="Completed...", value=3/3)

                        gVars$varISva <- varI
                        gVars$svaSV <- svaSV
                        gVars$svaSVc <- svaSVc
                        gVars$svaStep <- svaStep

                        #Update reactive variables associated to batch correction
                        gVars$expr.data <- gVars$norm.data
                        gVars$comb.data <- NULL
			gVars$agg.data <- NULL
                        gVars$conditions <- NULL
                        gVars$conditionsBoxplot <- NULL
                        gVars$correctionLvl <- NULL

                        tmpChoices <- rownames(gVars$expr.data)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

			shinyBS::toggleModal(session, "svaModal", toggle="close")

			shinyBS::updateButton(session, "launch_sva_modal", style="success", icon=icon("check-circle"))
                        if(input$corrType=="s"){
                                #if(arrType=="ag_exp2"){
				if(grepl("ag_exp.*", input$arrType)){
                                        shinyBS::updateCollapse(session, "bsSidebar", open="ANNOTATION", style=list("BATCH CORRECTION"="success", "ANNOTATION"="warning"))
                                        shinyBS::updateCollapse(session, "bsSidebar", style = list("DIFFERENTIAL ANALYSIS"="danger"))
                                }else{
                                        shinyBS::updateCollapse(session, "bsSidebar", open="DIFFERENTIAL ANALYSIS", style=list("BATCH CORRECTION"="success", "DIFFERENTIAL ANALYSIS"="warning"))
                                }
                        }

                        gVars$corrected <- TRUE

                        shiny::updateTabsetPanel(session, "display", selected="tvTab")
                        shiny::updateTabsetPanel(session, "mtv", selected="confPlotTab")

			shinyBS::updateButton(session, "launch_combat_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", style=list("DIFFERENTIAL ANALYSIS"="danger", "ANNOTATION"="danger"))
		})

		observeEvent(input$combat_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "BATCH CORRECTION")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

			if(is.null(gVars$norm.data) || is.null(input$batchCombat))
                        return(NULL)

                        #Clear SVA result reative variables
                        if(input$corrType=="c"){
                                gVars$varISva <- NULL
                                gVars$svaSV <- NULL
                                gVars$svaSVc <- NULL
                                gVars$svaStep <- NULL
                        }

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$corrected <- FALSE

                        progress$set(message="ComBat:", value=0)

			varI <- input$varICombat
                        coVar <- input$coVarCombat
                        if(is.null(coVar)){
                                coVar <- NULL
                        }else{
                                coVar <- as.list(coVar)
                        }
			batch <- as.list(input$batchCombat)
			batchCorVar <- list(var.int=varI, covariates=coVar, batches=batch)
			print("str(batchCorVar)")
			print(str(batchCorVar))

                        data <- gVars$norm.data
			phTable <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(phTable) <- phTable[,gVars$sampleColName]
                        rownames(phFactor) <- phTable[,gVars$sampleColName]

                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phTable <- cbind(phTable, svaSV, svaSVc)
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }
			arrType <- input$arrType
                        correctionLvl <- NULL

                        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~In ComBat Submit!!!!!!!!!")
                        updateProgress(detail="Correcting Batch Effects...", value=1/3)
                        npc <- 10
                        if(npc>ncol(data)){
                                npc <- ncol(data)
                        }
			comb.data <- remove.batch.effects(data, phFactor, npc, batchCorVar, method="Combat", plot=FALSE, verbose=TRUE)
                        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~After Correcting!!!!!!!!!")
                        if(is.character(comb.data)){
                                shinyjs::info(comb.data)
                                return(NULL)
                        }
                        expr.data <- comb.data
                        correctionLvl <- 1

                        agg.data <- NULL
                        map <- NULL
                        
                        updateProgress(detail="Completed!", value=3/3)
                        print("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Batch Corrected!!!!!!!!!")

			gVars$comb.data <- comb.data
			gVars$agg.data <- agg.data
			gVars$expr.data <- expr.data
                        gVars$varICombat <- varI
                        gVars$conditions <- levels(phFactor[,varI])
                        gVars$conditionsBoxplot <- levels(phFactor[,varI])
                        gVars$correctionLvl <- correctionLvl

                        tmpChoices <- rownames(gVars$expr.data)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

			shinyBS::updateButton(session, "launch_combat_modal", style="success", icon=icon("check-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))

                        #if(arrType=="ag_exp2"){
                        if(grepl("ag_exp.*", input$arrType)){
                                shinyBS::updateCollapse(session, "bsSidebar", open="ANNOTATION", style=list("BATCH CORRECTION"="success", "ANNOTATION"="warning"))
                                shinyBS::updateCollapse(session, "bsSidebar", style=list("DIFFERENTIAL ANALYSIS"="danger"))
                        }else{
                                shinyBS::updateCollapse(session, "bsSidebar", open="DIFFERENTIAL ANALYSIS", style=list("BATCH CORRECTION"="success", "DIFFERENTIAL ANALYSIS"="warning"))
                        }

                        gVars$corrected <- TRUE

			shinyBS::toggleModal(session, "combatModal", toggle="close")
                        shiny::updateTabsetPanel(session, "display", selected="tvTab")
                        shiny::updateTabsetPanel(session, "mtv", selected="mdsTab")
                        shiny::updateTabsetPanel(session, "mdsBox", selected="postCorTab")
		})

                observeEvent(input$corr_skip_submit, {
                        progress <- shiny::Progress$new()
                        updateProgress <- function(value=NULL, detail=NULL){
                                if (is.null(value)) {
                                        value <- progress$getValue()
                                        value <- value + (progress$getMax() - value) / 5
                                }
                                progress$set(value = value, detail = detail)
                        }
                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        print("In Skip Submit!!!")
                        if(is.null(gVars$expr.data)){
                                print("expr.data is null!!!")
                                return(NULL)
                        }

                        #Clear SVA result reative variables
                        gVars$varISva <- NULL
                        gVars$svaSV <- NULL
                        gVars$svaSVc <- NULL
                        gVars$svaStep <- NULL

                        #Clear differential analysis result reactive variables
                        gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        #Reset step status
                        gVars$corrected <- FALSE

			arrType <- input$arrType
                        print(paste0("arrType: ", arrType))

                        if(!grepl("ag_exp.*", input$arrType)){
                                shinyBS::updateCollapse(session, "bsSidebar", open="DIFFERENTIAL ANALYSIS", style=list("BATCH CORRECTION"="success", "DIFFERENTIAL ANALYSIS"="warning"))
                        }else{
                                shinyBS::updateCollapse(session, "bsSidebar", open="ANNOTATION", style=list("BATCH CORRECTION"="success", "ANNOTATION"="warning"))
                                shinyBS::updateCollapse(session, "bsSidebar", style=list("DIFFERENTIAL ANALYSIS"="danger"))
                                print("IN if context for ag2, updating collapse panel!")

                                #Update reactive variables associated to batch correction
                                gVars$agg.data <- NULL
                        }

                        #Update reactive variables associated to batch correction
                        gVars$expr.data <- gVars$norm.data
                        gVars$comb.data <- NULL
                        gVars$conditions <- NULL
                        gVars$conditionsBoxplot <- NULL
                        gVars$correctionLvl <- NULL

                        tmpChoices <- rownames(gVars$expr.data)
                        print("Exp heatmap choices :")
                        print(tmpChoices)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

                        gVars$corrected <- TRUE

			shinyBS::updateButton(session, "launch_sva_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_combat_modal", style="danger", icon=icon("exclamation-circle"))
			shinyBS::updateButton(session, "launch_ann_modal", style="danger", icon=icon("exclamation-circle"))
                })

                observeEvent(input$upload_ann_submit, {
                        updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "ANNOTATION & AGGREGATION")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        if(is.null(gVars$expr.data))
                        return(NULL)

                        #Clear differential analysis result reactive variables
			gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        shinyjs::html(id="loadingText", "ANNOTATION & AGGREGATION")

			progress <- shiny::Progress$new()

                        progress$set(message="Annotation:", value=0)

                        correctionLvl <- gVars$correctionLvl
			arrType <- input$arrType
                        if(is.null(correctionLvl)){
                                expr.data <- gVars$norm.data
                        }else if(correctionLvl==1){
                                expr.data <- gVars$comb.data
                        }else{
                                shinyjs::info("Could not locate expression data!")
                                return(NULL)
                        }
                        if(is.null(expr.data)){
                                shinyjs::info("Expression data is NULL!")
                                return(NULL)
                        }else if(nrow(expr.data)==0){
                                shinyjs::info("Expression data is empty!")
                                return(NULL)
                        }

                        agg.data <- NULL
                        map <- NULL

                        if(grepl("ag_exp.*", input$arrType)){
                                updateProgress(detail="Aggregating Probes by ID...", value=2/3)
				annDF <- gVars$annDF
				idIDX <- which(colnames(annDF)==input$ID)
				map <- annDF[,c(1,idIDX)]
                                
                                agg.data <- aggreg.probes.2(expr.data, map)
                                agg.data <- as.matrix(agg.data)
                                expr.data <- agg.data
                                correctionLvl <- 2
                        }
                        updateProgress(detail="Completed!", value=3/3)
			gVars$agg.data <- agg.data
			gVars$expr.data <- expr.data
                        gVars$map <- map

                        tmpChoices <- rownames(gVars$expr.data)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

			shinyBS::toggleModal(session, "importAnnModal", toggle="close")
			shinyBS::updateButton(session, "launch_ann_modal", style="success", icon=icon("check-circle"))
                        shinyBS::updateCollapse(session, "bsSidebar", open="DIFFERENTIAL ANALYSIS", style=list("ANNOTATION"="success", "DIFFERENTIAL ANALYSIS"="warning"))
                        shiny::updateTabsetPanel(session, "display", selected="tvTab")
                        shiny::updateTabsetPanel(session, "mtv", selected="mdsTab")
                        shiny::updateTabsetPanel(session, "mdsBox", selected="postAggTab")
                })

                observeEvent(input$submit_ann, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "ANNOTATION & AGGREGATION")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

                        if(is.null(gVars$expr.data))
                        return(NULL)

                        #Clear differential analysis result reactive variables
			gVars$deg.list <- NULL
                        gVars$deComps <- NULL

                        shinyjs::html(id="loadingText", "ANNOTATION & AGGREGATION")

                        progress$set(message="Annotation:", value=0)

                        correctionLvl <- gVars$correctionLvl
			arrType <- input$arrType
                        if(is.null(correctionLvl)){
                                expr.data <- gVars$norm.data
                        }else if(correctionLvl==1){
                                expr.data <- gVars$comb.data
                        }else{
                                shinyjs::info("Could not locate expression data!")
                                return(NULL)
                        }
                        if(is.null(expr.data)){
                                shinyjs::info("Expression data is NULL!")
                                return(NULL)
                        }else if(nrow(expr.data)==0){
                                shinyjs::info("Expression data is empty!")
                                return(NULL)
                        }

                        agg.data <- NULL

                        if(arrType=="ag_exp2"){
                                rgList <- gVars$rgList
                                print(str(rgList$genes))
                                ncIdx <- which(rgList$genes$ControlType==0)
                                map <- data.frame("ProbeName"=rgList$genes$ProbeName[ncIdx], "SystematicName"=rgList$genes$SystematicName[ncIdx], stringsAsFactors=FALSE)
                        }else if(arrType=="ag_exp1"){
                                map <- gVars$map
                        }

                        dups <- which(duplicated(map[,1]))
                        if(length(dups)>0){
                                map <- map[-dups,]
                        }
                        print(str(map))

                        agg.data <- aggreg.probes.2(expr.data, map)
                        agg.data <- as.matrix(agg.data)
                        expr.data <- agg.data
                        correctionLvl <- 2

                        updateProgress(detail="Completed!", value=2/2)

                        gVars$map <- map
                        gVars$agg.data <- agg.data
                        gVars$expr.data <- expr.data ## exprt.data here is agg.data, see above.

                        tmpChoices <- rownames(gVars$expr.data)
                        updateSelectizeInput(session=session, inputId='expGenes', choices=tmpChoices, server=TRUE)
                        gVars$expGenes <- NULL

                        shinyBS::updateCollapse(session, "bsSidebar", open="DIFFERENTIAL ANALYSIS", style=list("ANNOTATION"="success", "DIFFERENTIAL ANALYSIS"="warning"))
                        shiny::updateTabsetPanel(session, "display", selected="tvTab")
                        shiny::updateTabsetPanel(session, "mtv", selected="mdsTab")
                        shiny::updateTabsetPanel(session, "mdsBox", selected="postAggTab")
                })

		observeEvent(input$add_comp_submit, {
			if(is.null(gVars$condChoices())){
				gVars$comps <- list()
				return(NULL)
			}
			treat <- input$treatment
			control <- input$control
			comps <- paste0(treat, "-", control)
			compsSel <- comps
                        if(!is.null(input$comps)){
                                compsSel <- c(input$comps, compsSel)
                        }
			gVars$compsSel <- compsSel
			gVars$comps[[comps]] <- 1
		})

		gVars$compChoices <- reactive({
			if(length(gVars$comps)==0)
			return("NA")

			return(names(gVars$comps))
		})

		observeEvent(input$de_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "DIFFERENTIAL ANALYSIS")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

			arrType <- input$arrType
                        correctionLvl <- gVars$correctionLvl

                        if(is.null(gVars$expr.data))
                        return(NULL)

                        shinyjs::html(id="loadingText", "DIFFERENTIAL ANALYSIS")

                        progress$set(message="Differential Analysis:", value=0)

                        pvAdjMethod <- input$pvAdjMethod
			phTable <- gVars$phTable
                        sampleColName <- gVars$sampleColName
                        rownames(phTable) <- phTable[,sampleColName]
			phFactor <- gVars$phFactor
                        factorCols <- gVars$factorCols
                        varI <- input$varILimma
                        coVar <- input$coVarLimma
                        varCols <- append(varI, coVar)
                        varColsIdx <- NULL
                        if(length(varCols)>0){
                                varColsIdx <- which(colnames(phTable) %in% varCols)
                                print("varCols:")
                                print(varCols)
                                print(varColsIdx)
                        }

                        remSampleName <- NULL
                        remSampleIdx <- NULL
                        nrlevels.singular <- NULL
                        if(input$chkKeepNA && length(varColsIdx)>0){
                                remSampleName <- rownames(phTable)[apply(phTable[,varColsIdx, drop=FALSE], 1, function(x){res<-any(is.na(x));res})]
                                print("Sample names to remove:")
                                print(remSampleName)
                                if(length(remSampleName)>0){
                                        remSampleIdx <- which(rownames(phTable) %in% remSampleName)
                                        print("Sample Index to remove:")
                                        print(remSampleIdx)
                                        print("str(phTable) before removing:")
                                        print(str(phTable))
                                        phTable <- phTable[-remSampleIdx,]
                                        print("str(phTable) after removing:")
                                        print(str(phTable))
                                        print(phTable)
                                        shinyjs::alert(paste0("Following sample were removed before limma analysis :\n",paste(remSampleName,collapse=", "),"\n\nOne or more covariate was found to have NA values for these samples!!!."))

                                        #Remove columns with single level data
                                        nrlevels <- apply(phTable, 2, function(x){length(levels(factor(x)))})
                                        nrlevels.singular <- which(nrlevels<=1)
                                        print("nrlevels:")
                                        print(nrlevels)

                                        remInfo <- 0
                                        remStr <- ""
                                        if(length(nrlevels.singular)>0){
                                                remInfo <- 1
                                                remStr <- paste0(remStr, "Following columns are removed because they contain only single repeated value:\n[",paste0(names(nrlevels.singular), collapse=", "), "]")
                                                if(length(nrlevels.singular)==ncol(phTable)){
                                                        remStr <- paste0(remStr, "\n\nNo column survived filtering!!! Please define phenotype data columns with singular data type.")
                                                        shinyjs::info(remStr)
                                                        return(NULL)
                                                }
                                                col2rem <- which(colnames(phTable) %in% names(nrlevels.singular))
                                                phTable <- phTable[,-col2rem, drop=FALSE]
                                        }

                                        #Inform user with the columns removed from the data frame
                                        if(remInfo==1){
                                                shinyjs::info(remStr)
                                        }

                                        if(is.null(factorCols)){
                                                phFactor <- phTable
                                        }else{
                                                tmpIdx <- which(factorCols %in% colnames(phTable))
                                                if(length(tmpIdx)>0){
                                                        factorCols <- factorCols[tmpIdx]
                                                        factorIdx <- which(colnames(phTable) %in% factorCols)
                                                        phFactor <- factorize_cols(phTable, factorIdx)
                                                }else{
                                                        phFactor <- phTable
                                                }
                                        }
                                }
                        }

                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        if(length(remSampleIdx)>0){
                                                svaSV <- gVars$svaSV[-remSampleIdx,]
                                                svaSV <- as.data.frame(apply(svaSV, 2, factor))
                                                svaSVc <- gVars$svaSVc[-remSampleIdx,]
                                        }else{
                                                svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                                svaSVc <- gVars$svaSVc
                                        }
                                                phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }

                        if(is.null(coVar)){
                                coVar <- NULL
                        }else{
                                coVar2rem <- NULL
                                if(length(nrlevels.singular)>0){
                                        coVar2rem <- which(coVar %in% names(nrlevels.singular))
                                }

                                if(length(coVar2rem)>0){
                                        coVar2remNames <- coVar[coVar2rem]
                                        coVar <- as.list(coVar[-coVar2rem])
                                        shinyjs::alert(paste0("Following User specified CoVariate(s) were removed :\n",paste(coVar2remNames,collapse=", "),"\n\nFiltered in previous check!!!."))
                                }else{
                                    coVar <- as.list(coVar)
                                }
                        }

			batchCorVar <- list(var.int=varI, covariates=coVar) ## No batch for limma model
			print(str(batchCorVar))

                        data <- gVars$expr.data
                        if(length(remSampleName)>0){
                                data <- data[,-remSampleIdx]
                        }

			updateProgress(detail="Building Model...", value=1/3)
                        print(str(batchCorVar))
                        print("str(phTable) : ")
                        print(str(phTable))
                        print("str(phFactor) : ")
                        print(str(phFactor))
                        print(phFactor)
			des <- build.model.matrix(phFactor, intercept=-1, batchCorVar$var.int, batchCorVar$covariates, verbose=T)
                        ne <- limma::nonEstimable(des)
                        if(!is.null(ne)){
                                shinyjs::alert(paste0("Coefficients not estimable : ",paste(ne,collapse=", "),"\n\nPlease check your limma model definition!!! Possibly contains confounders."))
                                return(NULL)
                        }

			annDF <- gVars$map
                        dupIdx <- which(duplicated(annDF[,2])==T)
                        if(length(dupIdx)>0){
                                annDF <- annDF[-which(duplicated(annDF[,2])==T),]
                        }
			rownames(annDF) <- annDF[,2]

			updateProgress(detail="Evaluating Expression...", value=2/3)

                        comps <- as.vector(input$comps)

                        #Check and filter comps
                        compCheck <- sapply(strsplit(comps, "-"), function(x){all(x %in% phTable[,varI])})
                        if(sum(compCheck)==0){
                                shinyjs::alert(paste0("Could not find any valid comparisons:\n\nFollowing comparisons were removed:\n",paste(comps[!compCheck],collapse=", "),"\n\nSelected Variable of Interest - '", varI,"'\n\nPlease check selected varaiable of interest and selected comparisons."))
                                return(NULL)
                        }else if(sum(compCheck)<length(comps)){
                                shinyjs::alert(paste0("Following comparisons were removed:\n", paste(comps[!compCheck],collapse=", "),"\n\nRemoved comparison does not match varaiable of interest."))
                        }
                        comps <- comps[compCheck]

                        comps <- sapply(strsplit(comps, "-"), function(x){res<-make.names(x);res<-paste0(res, collapse="-")})
			deg.list <- diff.gene.expr(data, des, contrasts=comps, pvalue=1, fcvalue=0, p.adjust.method=pvAdjMethod, annot=annDF, plot=F, verbose=T)

                        print("Filtering Differential Expresssion Table...")
			comp <- comps[1]
			print(paste0("COMP SEL: ", comp))
			degDF <- deg.list[[comp]]

			lfc <- as.numeric(input$lfcThr)
			pvThr <- -log10(as.numeric(input$pvThr))
                        pvType <- input$pvType
			selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)
			if(length(selVec)==0){
                                degDF <- NULL
                        }else{
                                degDF <- degDF[selVec,]
                        }
			gVars$filteredDeTable <- degDF
                        gVars$compFilt <- comp
			gVars$pvFilt <- pvThr
			gVars$pvType <- pvType
			gVars$lfcFilt <- input$lfcThr

                        updateProgress(detail="Completed!", value=3/3)

			gVars$deg.list <- deg.list
                        gVars$deComps <- comps
                        gVars$varI <- varI
                        gVars$conditions <- levels(phFactor[,varI])
                        gVars$conditionsBoxplot <- levels(phFactor[,varI])

                        gVars$intersectComps <- comps
                        gVars$heatComps <- comps[1]
                        gVars$conditionsHeat <- levels(phFactor[,varI])
                        gVars$percDE <- input$percDE

                        shiny::updateTabsetPanel(session, "display", selected="diffTab")
                        shiny::updateTabsetPanel(session, "diffTBox", selected="diffSummTableTab")
                        shinyBS::updateCollapse(session, "bsSidebar", style=list("DIFFERENTIAL ANALYSIS"="success"))
		})

		output$preBoxPlot <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
			filt.data <- gVars$filt.data
                        if(input$arrType=="il_methyl"){
                                filt.data <- minfi::preprocessRaw(filt.data)
                                filt.data <- minfi::getMethSignal(filt.data, "M")
			        boxplot(filt.data, las=2, cex=0.7, main="Before Normalization")
                        }else{
			        boxplot(log2(filt.data), las=2, cex=0.7, main="Before Normalization")
                        }
		})

		output$postBoxPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
			norm.data <- gVars$norm.data
			boxplot(norm.data, las=2, cex=0.7, main="After Normalization")
		})

		output$preDensityPlot <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(input$arrType!="ag_exp1", "No plot for Agilent single color data!"))
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
                        if(input$arrType=="il_methyl"){
                                RGset <- gVars$RGset
                                minfi::densityPlot(RGset, main="Before Normalization", xlab="Beta")
                        }else{
                                rgList <- gVars$rgList
                                limma::plotDensities(rgList, main="Before Normalization")
                        }
		})

		output$postDensityPlot <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(input$arrType!="ag_exp1", "No plot for Agilent single color data!"))
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
                        if(input$arrType=="il_methyl"){
                                Mset <- gVars$Mset
                                shiny::validate(
                                        need(any(is(Mset, "RGChannelSet"), is(Mset, "MethylSet")), "No plot for this nomalization method!")
                                )
                                minfi::densityPlot(Mset, main="After Normalization", xlab="Beta")
                        }else{
				shiny::validate(need(!is.null(gVars$rgList.norm), "RG channel object is null for normalized data..."))
                                rgList.norm <- gVars$rgList.norm
                                limma::plotDensities(rgList.norm, main="After Normalization")
                        }
		})

		output$preMDPlot <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(input$arrType!="ag_exp1", "No plot for Agilent single color data!"))
                        shiny::validate(need(input$arrType!="il_methyl", "No plot for Illumina methylation data!"))
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
			rgList <- gVars$rgList
			limma::plotMD(rgList, main="Before Normalization")
		})

		output$postMDPlot <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(input$arrType!="ag_exp1", "No plot for Agilent single color data!"))
                        shiny::validate(need(input$arrType!="il_methyl", "No plot for Illumina methylation data!"))
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
                        shiny::validate(need(!is.null(gVars$rgList.norm), "RG channel object is null for normalized data..."))
			rgList.norm <- gVars$rgList.norm
			limma::plotMD(rgList.norm, main="After Normalization")
		})

		output$confPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))

			ph <- gVars$phTable
			phFactor <- gVars$phFactor
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
                        phFactor <- phFactor[,test]

                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }
                        print("Printing str(phFactor):")
                        print(str(phFactor))
			swamp::confounding(phFactor, margins=c(10,10))
		},
		height=function(){
			if(is.null(input$confPlotDiv_height) || input$confPlotDiv_height=="" || input$confPlotDiv_height==0){
				return((40*(session$clientData$output_confPlot_width*100/60))/100)
			}else{
				input$confPlotDiv_height
			}
		},
		width = function(){
			if(is.null(input$confPlotDiv_width) || input$confPlotDiv_width=="" || input$confPlotDiv_width==0){
				return(session$clientData$output_confPlot_width)
			}else{
				input$confPlotDiv_width
			}
		})

		output$princePlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
			#print("In Prince Plot...")
			#print(paste0("Plot Width: ", session$clientData$output_princePlot_width))
			#print(paste0("Plot Height: ", session$clientData$output_princePlot_height))
			#print(paste0("Param Height: ", input$princePlotDiv_height))
			#print(paste0("Param Width: ", input$princePlotDiv_width))

			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$norm.data
			npc = 10
			test <- sapply(colnames(ph), function(b){length(table(ph[,b]))>1 && length(table(ph[,b]))!=length(ph[,b])})
                        print("test cols:")
                        print(test)
                        phFactor <- phFactor[,test]
                        print("str(data):")
                        print(str(data))
                        print("str(ph):")
                        print(str(ph))
                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }
                        if(npc>ncol(data)){
                                npc <- ncol(data)
                        }
			pr <- swamp::prince(data, phFactor, top=npc)
			# generate the prince plot
			swamp::prince.plot(prince=pr, margins=c(15,15), note=TRUE)
		},
		height=function(){
			if(is.null(input$princePlotDiv_height) || input$princePlotDiv_height=="" || input$princePlotDiv_height==0){
				#print("Ht1")
				#print(session$clientData$output_princePlot_width)
				return((40*(session$clientData$output_princePlot_width*100/60))/100)
			}else{
				#print("Ht2")
				input$princePlotDiv_height
			}
		},
		width = function(){
			if(is.null(input$princePlotDiv_width) || input$princePlotDiv_width=="" || input$princePlotDiv_width==0){
				#print("Wd1")
				#print(session$clientData$output_princePlot_width)
				return(session$clientData$output_princePlot_width)
			}else{
				#print("Wd2")
				input$princePlotDiv_width
			}
		}
		)

		output$hcPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
			cat("Accessing PH from extracted data...")
			cat("\n")

			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$norm.data
                        colnames(data) <- rownames(ph)
			print(dim(data))
			print(colnames(data))
			print(dim(ph))
			print(rownames(ph))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
                        phFactor <- phFactor[,test]
                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }
                        print("phFactor:")
                        print(phFactor)
                        print("Str phFactor:")
                        print(str(phFactor))
			# hc plot
			swamp::hca.plot(data, phFactor, method="correlation")
		},
		height=function(){
			if(is.null(input$hcPlotDiv_height) || input$hcPlotDiv_height==""){
				return((40*(session$clientData$output_hcPlot_width*100/60))/100)
			}else{
				input$hcPlotDiv_height
			}
		},
		width=function(){
			if(is.null(input$hcPlotDiv_width) || input$hcPlotDiv_width==""){
				return(session$clientData$output_hcPlot_width)
			}else{
				input$hcPlotDiv_width
			}
		})

                gVars$mdsLblCols <- reactive({
                        if(is.null(input$mdsColor) || is.null(gVars$phFactor)){
                                return(NULL)
                        }
			phFactor <- gVars$phFactor
			colVec <- get_color_palette(iVec=phFactor[,input$mdsColor], asFactor=TRUE)
                        return(colVec)
                })

		output$preCorMDS <- renderPlot({
                        shiny::validate(need(!is.null(gVars$norm.data), "Waiting for normalization..."))
                        shiny::validate(need(!is.null(input$mdsLabel), "Waiting for label input..."))
                        shiny::validate(need(!is.null(input$mdsColor), "Waiting for color input..."))

			norm.data <- gVars$norm.data
			phTable <- gVars$phTable
			phFactor <- gVars$phFactor

                        print(colnames(norm.data))
                        print(phTable[,input$mdsLabel])
                        print(as.character(phTable[,input$mdsLabel]))

                        percTop <- 1
                        topNum <- floor((nrow(norm.data)*percTop)/100)

                        cat("Total Features: ", nrow(norm.data), "\n")
                        cat("Top Num: ", topNum, "\n")

                        #Get unique vector of colors for MDS plot
			colVec <- gVars$mdsLblCols()

			limma:::plotMDS(norm.data, top=topNum, labels=as.character(phFactor[,input$mdsLabel]), col=colVec, gene.selection="common", main="Before removing any batch", font=2)
                },
		height = function(){
			if(is.null(input$preCorMDSDiv_height) || input$preCorMDSDiv_height==""){
				return((40*(session$clientData$output_preCorMDS_width*100/60))/100)
			}else{
				input$preCorMDSDiv_height
			}
		},
		width = function(){
			if(is.null(input$preCorMDSDiv_width) || input$preCorMDSDiv_width==""){
				return(session$clientData$output_preCorMDS_width)
			}else{
				input$preCorMDSDiv_width
			}
		})

		output$postCorMDS <- renderPlot({
                        shiny::validate(need(!is.null(gVars$comb.data), "Waiting for batch correction..."))
                        shiny::validate(need(!is.null(input$mdsLabel), "Waiting for label input..."))
                        shiny::validate(need(!is.null(input$mdsColor), "Waiting for color input..."))
			comb.data <- gVars$comb.data
			phTable <- gVars$phTable
			phFactor <- gVars$phFactor
                        percTop <- 1
                        topNum <- floor((nrow(comb.data)*percTop)/100)

                        cat("Total Features: ", nrow(comb.data), "\n")
                        cat("Top Num: ", topNum, "\n")

                        #Get unique vector of colors for MDS plot
			colVec <- gVars$mdsLblCols()

			limma:::plotMDS(comb.data, top=topNum, labels=phFactor[,input$mdsLabel], col=colVec, gene.selection="common", main="After removing batches", font=2)
		},
		height=function(){
			if(is.null(input$postCorMDSDiv_height) || input$postCorMDSDiv_height==""){
				return((40*(session$clientData$output_postCorMDS_width*100/60))/100)
			}else{
				input$postCorMDSDiv_height
			}
		},
		width=function(){
			if(is.null(input$postCorMDSDiv_width) || input$postCorMDSDiv_width==""){
				return(session$clientData$output_postCorMDS_width)
			}else{
				input$postCorMDSDiv_width
			}
		})

		output$postAggMDS <- renderPlot({
                        shiny::validate(need(input$arrType!="af_exp", "No plot for Affymetrix data!"))
                        shiny::validate(need(input$arrType!="il_methyl", "No plot for Illumina methylation data!"))
                        shiny::validate(need(!is.null(gVars$agg.data), "Waiting for batch correction..."))
			agg.data <- gVars$agg.data
			phTable <- gVars$phTable
			phFactor <- gVars$phFactor

                        percTop <- 1
                        topNum <- floor((nrow(agg.data)*percTop)/100)

                        cat("Total Features: ", nrow(agg.data), "\n")
                        cat("Top Num: ", topNum, "\n")

                        #Get unique vector of colors for MDS plot
			colVec <- gVars$mdsLblCols()

			limma:::plotMDS(agg.data, top=topNum, labels=phFactor[,input$mdsLabel], col=colVec, gene.selection="common", main="After aggregation", font=2)
		},
		height=function(){
			if(is.null(input$postAggMDSDiv_height) || input$postAggMDSDiv_height==""){
				return((40*(session$clientData$output_postAggMDS_width*100/60))/100)
			}else{
				input$postAggMDSDiv_height
			}
		},
		width=function(){
			if(is.null(input$postAggMDSDiv_width) || input$postAggMDSDiv_width==""){
				return(session$clientData$output_postAggMDS_width)
			}else{
				input$postAggMDSDiv_width
			}
		})

		output$postHcPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$comb.data), "Waiting for batch correction..."))
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$comb.data
                        colnames(data) <- rownames(ph)
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))

                        phFactor <- phFactor[,test]
                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }

			# hc plot
			swamp::hca.plot(data, phFactor, method = "correlation")
		},
		height=function(){
			if(is.null(input$postHcPlotDiv_height) || input$postHcPlotDiv_height==""){
				return((40*(session$clientData$output_postHcPlot_width*100/60))/100)
			}else{
				input$postHcPlotDiv_height
			}
		},
		width=function(){
			if(is.null(input$postHcPlotDiv_width) || input$postHcPlotDiv_width==""){
				return(session$clientData$output_postHcPlot_width)
			}else{
				input$postHcPlotDiv_width
			}
		})

                output$postSvaHcPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$comb.sva.data), "No unknown batches removed..."))
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$comb.sva.data
                        colnames(data) <- rownames(ph)
			print(head(data))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))

                        phFactor <- phFactor[,test]
                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }

			# hc plot
			swamp::hca.plot(data, phFactor, method = "correlation")
		},
		height=function(){
			if(is.null(input$postSvaHcPlotDiv_height) || input$postSvaHcPlotDiv_height==""){
				return((40*(session$clientData$output_postSvaHcPlot_width*100/60))/100)
			}else{
				input$postSvaHcPlotDiv_height
			}
		},
		width=function(){
			if(is.null(input$postSvaHcPlotDiv_width) || input$postSvaHcPlotDiv_width==""){
				return(session$clientData$output_postSvaHcPlot_width)
			}else{
				input$postSvaHcPlotDiv_width
			}
		})

		output$postPrincePlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$comb.data), "Waiting for normalization..."))
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$comb.data
			npc = 10
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
                        phFactor <- phFactor[,test]
                        if(!is.null(gVars$svaStep)){
                                if(gVars$svaStep==1){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phFactor <- cbind(phFactor, svaSV, svaSVc)
                                }
                        }
                        if(npc>ncol(data)){
                                npc <- ncol(data)
                        }
			pr <- swamp::prince(data, phFactor, top=npc)
			# generate the prince plot
			swamp::prince.plot(prince=pr, margins=c(15,15), note=TRUE)
		},
		height=function(){
			if(is.null(input$postPrincePlotDiv_height) || input$postPrincePlotDiv_height==""){
				return((40*(session$clientData$output_postPrincePlot_width*100/60))/100)
			}else{
				input$postPrincePlotDiv_height
			}
		},
		width=function(){
			if(is.null(input$postPrincePlotDiv_width) || input$postPrincePlotDiv_width==""){
				return(session$clientData$output_postPrincePlot_width)
			}else{
				input$postPrincePlotDiv_width
			}
		})

                output$postSvaPrincePlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$comb.sva.data), "No unknown batches removed..."))
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			phFactor <- gVars$phFactor
                        rownames(ph) <- ph[,gVars$sampleColName]
                        rownames(phFactor) <- ph[,gVars$sampleColName]
			data <- gVars$comb.sva.data
			npc = 10
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
                        phFactor <- phFactor[,test]
                        if(npc>ncol(data)){
                                npc <- ncol(data)
                        }
			pr <- swamp::prince(data, phFactor, top=npc)
			# generate the prince plot
			swamp::prince.plot(prince=pr, margins=c(15,15), note=TRUE)
		},
		height=function(){
			if(is.null(input$postSvaPrincePlotDiv_height) || input$postSvaPrincePlotDiv_height==""){
				return((40*(session$clientData$output_postSvaPrincePlot_width*100/60))/100)
			}else{
				input$postSvaPrincePlotDiv_height
			}
		},
		width=function(){
			if(is.null(input$postSvaPrincePlotDiv_width) || input$postSvaPrincePlotDiv_width==""){
				return(session$clientData$output_postSvaPrincePlot_width)
			}else{
				input$postSvaPrincePlotDiv_width
			}
		})

                observeEvent(input$filterDE_submit, {
                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}

			#start loading screen
			shinyjs::html(id="loadingText", "FILTERING DIFFERENTIAL TABLE")
			shinyjs::show(id="loading-content")

                        on.exit({
                                progress$close()
                                shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                        })

			if(is.null(gVars$deg.list) || is.null(input$compDE))
			return(NULL)

                        shinyjs::html(id="loadingText", "FILTERING DIFFERENTIAL TABLE")

                        progress$set(message="Filtering DE:", value=0)

                        updateProgress(detail="...", value=1/2)
                	deg.list <- gVars$deg.list
			comp <- input$compDE
			print(paste0("COMP SEL: ", input$compDE))
			degDF <- deg.list[[comp]]

			lfc <- as.numeric(input$lfcThr)
			pvThr <- -log10(as.numeric(input$pvThr))
                        pvType <- input$pvType
			selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)
			if(length(selVec)==0){
                                degDF <- NULL
                        }else{
                                degDF <- degDF[selVec,]
                        }

			gVars$filteredDeTable <- degDF
                        gVars$compFilt <- comp
			gVars$pvFilt <- pvThr
			gVars$pvType <- pvType
			gVars$lfcFilt <- input$lfcThr
		})

                output$deSummTable <- DT::renderDataTable({
			shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))
			shiny::validate(need(length(gVars$deg.list)>0, "Waiting for Differential Analysis Results..."))

			lfc <- as.numeric(input$lfcThr)
                	deg.list <- gVars$deg.list

                        summDF <- get_deg_summary(deg_list=deg.list, names=names(deg.list), lfc=lfc)
			DT::datatable(summDF, filter=list(position='top', clear=FALSE), options=list(search=list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		},server=TRUE)

		output$deTable <- DT::renderDataTable({
			shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))
			shiny::validate(need(length(gVars$deg.list)>0, "Waiting for Differential Analysis Results..."))
			shiny::validate(need(!is.null(gVars$filteredDeTable), "Filtered Differential Table is Empty!..."))

			arrType <- input$arrType
			degDF <- gVars$filteredDeTable
                        if(arrType=="il_methyl"){
				if(gVars$RGarray=="EPIC"){
                                        data(Other, package="IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
                                }else if(gVars$RGarray=="450k"){
                                        data(Other, package="IlluminaHumanMethylation450kanno.ilmn12.hg19")
                                }
                                degDF <- as.data.frame(cbind(degDF, Other[degDF$ID, c(3:ncol(Other))]), stringsAsFactors=FALSE)
				rm(Other)
                                gc()
                        }
			DT::datatable(degDF, filter=list(position='top', clear=FALSE), options=list(search=list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		},server=TRUE)

                output$exportDE <- shiny::downloadHandler(
			filename = function(){
				paste("Differential_Expression_Tables_", Sys.Date(), '.xlsx', sep='')
			},
			content = function(con){
				#screen loading screen
				shinyjs::html(id="loadingText", "EXPORTING DIFFERENTIAL TABLE")
				shinyjs::show(id="loading-content")

                                on.exit({
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                })
                                deg.list <- gVars$deg.list
				lfc <- as.numeric(input$lfcThr)
                                pvThr <- gVars$pvFilt
                                pvType <- gVars$pvType

				#Update Names
				names(deg.list) <- strtrim(paste0(c(1:length(deg.list)), ".", names(deg.list)), 30)
                                if(input$chkExportDE==TRUE){
                                        deg.list.filt <- deg.list
                                        summDF <- get_deg_summary(deg_list=deg.list, names=names(deg.list), lfc=0)
                                }else{
                                        deg.list.filt <- lapply(deg.list, function(degDF){
                                                selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)
                                                if(length(selVec>0)){
                                                        degDF[selVec,]
                                                }else{NULL}
                                        })
                                        summDF <- get_deg_summary(deg_list=deg.list, names=names(deg.list), lfc=lfc)
                                }
                                print("str(summDF)")
                                print(str(summDF))
                                print("str(deg.list.filt) before filtering")
                                print(str(deg.list.filt))
                                deg.list.filt <- c("Summary"=list(summDF), deg.list.filt)
				deg.list.filt <- deg.list.filt[unlist(lapply(deg.list.filt, function(degDF) !is.null(degDF)))]
				#Add annotation columns for Methylation
				if(input$arrType=="il_methyl"){
					if(gVars$RGarray=="EPIC"){
						data(Other, package="IlluminaHumanMethylationEPICanno.ilm10b2.hg19")
					}else if(gVars$RGarray=="450k"){
						data(Other, package="IlluminaHumanMethylation450kanno.ilmn12.hg19")
					}
					deg.list.filt <- lapply(deg.list.filt, function(degDF){as.data.frame(cbind(degDF, Other[degDF$ID, c(3:ncol(Other))]), stringsAsFactors=FALSE)})
					rm(Other)
                                        gc()
				}
                                print("str(deg.list.filt) after filtering")
                                print(str(deg.list.filt))
				print("DE GeneList Names:")
				print(names(deg.list.filt))
                                WriteXLS::WriteXLS(deg.list.filt, ExcelFileName=con, col.names=T, AdjWidth=T, BoldHeaderRow=T)
			}
		)

                output$exportNormMat <- shiny::downloadHandler(
			filename = function(){
				paste("Expression_Matrix_Normalized_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				#screen loading screen
				shinyjs::html(id="loadingText", "EXPORTING NORAMLIZED MATRIX")
				shinyjs::show(id="loading-content")

                                on.exit({
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                })
                                data <- gVars$norm.data
				write.table(data, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

                output$exportFiltMat <- shiny::downloadHandler(
			filename = function(){
				paste("Expression_Matrix_Filtered_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				#screen loading screen
				shinyjs::html(id="loadingText", "EXPORTING FILTERED MATRIX")
				shinyjs::show(id="loading-content")

                                on.exit({
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                })
				if(input$arrType=="il_methyl"){
					if(is.null(gVars$methFilt.data)){
						data <- "Something went wrong while preparing filtered object..."
					}else{
						data <- gVars$methFilt.data
					}
				}else{
					data <- gVars$filt.data
				}
				write.table(data, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

                output$exportCorrMat <- shiny::downloadHandler(
			filename = function(){
				paste("Expression_Matrix_Corrected_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				#screen loading screen
				shinyjs::html(id="loadingText", "EXPORTING CORRECTED MATRIX")
				shinyjs::show(id="loading-content")

                                on.exit({
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                })
                                data <- gVars$comb.data
				write.table(data, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

                output$exportAggMat <- shiny::downloadHandler(
			filename = function(){
				paste("Expression_Matrix_Aggregated_", Sys.Date(), '.txt', sep='')
			},
			content = function(con){
				#screen loading screen
				shinyjs::html(id="loadingText", "EXPORTING AGGREGATED FEATURES MARIX")
				shinyjs::show(id="loading-content")

                                on.exit({
                                        shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
                                })
                                data <- gVars$agg.data
				write.table(data, con, row.names=TRUE, col.names=TRUE, quote=FALSE, sep="\t")
			}
		)

		observeEvent(input$generateRpt, {
			#start loading screen
			shinyjs::html(id="loadingText", "CREATING ANALYSIS REPORT")
			shinyjs::show(id="loading-content")

			on.exit({
				shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
			})
			#Disable Warning
			oldw <- getOption("warn")
			options(warn = -1)

			tempReport <- file.path(tempdir(), "report.Rmd")
			file.copy("report.Rmd", tempReport, overwrite=TRUE)

			exn <- ".pdf"
			con <- tempfile(pattern="eUTOPIA_Analysis_Report_", fileext=exn)

			message("generateRpt srcDir:")
			message(srcDir)

			#params <- list(gVars=gVars, input=input)
			params <- list(gVars=gVars, input=input, srcDir=srcDir)
			rmarkdown::render(tempReport, output_file=con,
				params=params,
				envir=new.env(parent=globalenv())
			)

			#Enable Warning
			options(warn = oldw)

			gVars$analysisReport_created <- TRUE
			gVars$analysisReport_con <- con

			shinyBS::toggleModal(session, "analysisReportModal", toggle="open")
		})

                output$exportRpt <- shiny::downloadHandler(
			filename = function(){
				paste("eUTOPIA_Analysis_Report_", Sys.Date(), '.pdf', sep='')
			},
			content = function(con){
				if(gVars$analysisReport_created==TRUE){
					file.copy(gVars$analysisReport_con, con)
				}
			}
		)

                gVars$volGGplot <- reactive({
                        if(is.null(gVars$deg.list)){
                               return(NULL)
                        }
                        pvThr <- gVars$pvFilt
                        pvType <- gVars$pvType
                        deg.list <- gVars$deg.list
			comp <- gVars$compFilt
			deg <- deg.list[[comp]]
			deg <- data.frame(x=as.numeric(deg$logFC), y=-log10(as.numeric(deg[,pvType])), ID=rownames(deg))
			p <- ggplot(deg, aes(x, y, label= ID)) + geom_point() +
			geom_vline(xintercept = input$lfcThr, color = "blue") +
			geom_vline(xintercept = -input$lfcThr, color = "blue") +
			geom_hline(yintercept = pvThr, color = "red") +
			labs(x="log2(Fold-change)", y="-log10(P.Value)")
                        gVars$volDeg <- deg
                        return(p)
                })

		output$volcanoPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))

                        input$vol_submit

                        pvThr <- gVars$pvFilt
                        lfc <- gVars$lfcFilt

                        p <- gVars$volGGplot()
                        deg <- gVars$volDeg
                        print("str(deg)")
                        print(str(deg))
                        degSel <- deg[deg$y>=pvThr & abs(deg$x)>=lfc,]

                        isolate({
                                p <- p +
                                geom_point(data = degSel, color = "red") +
                                scale_x_continuous("log2(Fold-change)", limits = range(-input$xAxis,input$xAxis)) +
                                scale_y_continuous("-log10(P.Value)", limits = range(0,input$yAxis)) + theme_bw()

                                if(input$gName){
                                    p <- p + geom_text_repel(
                                            data = degSel,
                                            aes(label = ID),
                                            size = 5
                                    )
                                }
                        })
			print(p)
		}, height=900, width=900)

                observeEvent(input$venn_submit, {
                        gVars$intersectComps <- input$intersectComps
                })

                output$intersectPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))

                        input$venn_submit

			deg.list <- gVars$deg.list
			lfc <- as.numeric(gVars$lfcFilt)
			pvThr <- gVars$pvFilt
			pvType <- gVars$pvType

                        isolate({
                                intersectComps <- input$intersectComps

                                deg.item.list <- list()
                                for(comp in intersectComps){
                                        degDF <- deg.list[[comp]]
                                        selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)

                                        if(length(selVec)==0)
                                        next

                                        deg.item.list[[comp]] <- rownames(degDF[selVec,])
                                }

                                shiny::validate(need(length(deg.item.list)>1, "Not Enough Sets to Intersect!"))

                                print(str(deg.item.list))
                                if(length(deg.item.list)<=4){
                                        fillCol <- c("blue", "red", "yellow", "green")
                                        cat_label_size <- 1
                                        area_label_size <- 2
                                        cat_dist <- c(0.22,0.22,0.1,0.1)
                                        printMode <- "raw"
                                        marginSize <- 0.12
                                        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                                        venn.plot <- VennDiagram::venn.diagram(
                                                  x=deg.item.list,
                                                  filename=NULL,
                                                  cat.cex=cat_label_size,
                                                  cat.dist=cat_dist[1:length(deg.item.list)],
                                                  cex=area_label_size,
                                                  fill=fillCol[1:length(deg.item.list)],
                                                  margin=marginSize,
                                                  print.mode=printMode
                                        )
                                        grid.draw(venn.plot)
                                }else{
                                        deg.item.DF <- data.frame(genes=unique(unlist(deg.item.list)), stringsAsFactors=F)
                                        for(itm in names(deg.item.list)){
                                                deg.item.DF[,itm] <- 0
                                                deg.item.DF[which(deg.item.DF$genes %in% deg.item.list[[itm]]),itm] <- 1
                                        }
                                        setNames <- names(deg.item.list)
                                        setColors <- rep("gray23", length(deg.item.list))
                                        keepingOrder <- F
                                        textScale <- 1
                                        intersectSize <- 10
                                        xLabel <- "Set Size"
                                        UpSetR::upset(deg.item.DF, nintersects=intersectSize, sets=setNames, order.by="freq", decreasing=T, line.size=0.1, sets.x.label=xLabel, keep.order=keepingOrder, text.scale=textScale, sets.bar.color=setColors)
                                }
                        })
		}, height=900, width=900)

                #Venn diagram for set intersection
                output$vennPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))

                        input$venn_submit

			deg.list <- gVars$deg.list
			lfc <- as.numeric(gVars$lfcFilt)
			pvThr <- gVars$pvFilt
			pvType <- gVars$pvType

                        intersectComps <- gVars$intersectComps

                        deg.item.list <- list()
                        for(comp in intersectComps){
                                degDF <- deg.list[[comp]]
                                selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)

                                if(length(selVec)==0)
                                next

                                deg.item.list[[comp]] <- rownames(degDF[selVec,])
                        }

                        shiny::validate(need(length(deg.item.list)>1, "Not Enough Sets to Intersect!"))
                        shiny::validate(need(length(deg.item.list)<=4, "Too Many Sets for Venn Representation!"))

                        print(str(deg.item.list))
                        #Create venn diagram
                        fillCol <- c("blue", "red", "yellow", "green")
                        cat_label_size <- 1
                        area_label_size <- 2
                        cat_dist <- c(0.22,0.22,0.1,0.1)
                        printMode <- "raw"
                        marginSize <- 0.12
                        futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")
                        venn.plot <- VennDiagram::venn.diagram(
                                  x=deg.item.list,
                                  filename=NULL,
                                  cat.cex=cat_label_size,
                                  cat.dist=cat_dist[1:length(deg.item.list)],
                                  cex=area_label_size,
                                  fill=fillCol[1:length(deg.item.list)],
                                  margin=marginSize,
                                  print.mode=printMode
                        )

                        grid.draw(venn.plot)
		}, height=900, width=900)

                #UpSet plot for set intersection
                output$upsetPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))

                        input$venn_submit

			deg.list <- gVars$deg.list
			lfc <- as.numeric(gVars$lfcFilt)
			pvThr <- gVars$pvFilt
			pvType <- gVars$pvType

                        intersectComps <- gVars$intersectComps

                        deg.item.list <- list()
                        for(comp in intersectComps){
                                degDF <- deg.list[[comp]]
                                selVec <- which(-log10(degDF[,pvType])>pvThr & abs(degDF$logFC)>lfc)

                                if(length(selVec)==0)
                                next

                                deg.item.list[[comp]] <- rownames(degDF[selVec,])
                        }

                        shiny::validate(need(length(deg.item.list)>1, "Not Enough Sets to Intersect!"))

                        print(str(deg.item.list))
                        #Create UpSet plot
                        deg.item.DF <- data.frame(genes=unique(unlist(deg.item.list)), stringsAsFactors=F)
                        for(itm in names(deg.item.list)){
                                deg.item.DF[,itm] <- 0
                                deg.item.DF[which(deg.item.DF$genes %in% deg.item.list[[itm]]),itm] <- 1
                        }
                        setNames <- names(deg.item.list)
                        setColors <- rep("gray23", length(deg.item.list))
                        keepingOrder <- F
                        textScale <- 1
                        intersectSize <- 10
                        xLabel <- "Set Size"

                        UpSetR::upset(deg.item.DF, nintersects=intersectSize, sets=setNames, order.by="freq", decreasing=T, line.size=0.1, sets.x.label=xLabel, keep.order=keepingOrder, text.scale=textScale, sets.bar.color=setColors)
		}, height=900, width=900)

                observeEvent(input$expBoxplot_submit, {
                        gVars$varIBoxPlot <- input$varIBoxPlot
                        gVars$expGenes <- input$expGenes
                        gVars$conditionsBoxplot <- input$conditionsBoxplot
                })

                output$expressionBoxPlot <- renderPlot({
                        shiny::validate(need(!is.null(gVars$expr.data), "Missing expression data..."))

                        data <- gVars$expr.data
                        phTable <- gVars$phTable

                        varI <- gVars$varIBoxPlot
                        genes <- gVars$expGenes
                        conditions <- gVars$conditionsBoxplot

                        isCond <- FALSE
                        if(!is.null(conditions)){
                                if(!is.na(conditions) && conditions!="NA")
                                isCond <-TRUE
                        }

                        shiny::validate(need(isCond, "Need at least 1 condition to create boxplot!"))
                        shiny::validate(need(length(genes)>=1, "Need at least 1 gene to create boxplot!"))

                        dataCols <- as.integer(gsub("X", "", colnames(data)))
                        print(dataCols)
                        dataColsVarI <- phTable[,varI]

                        tmp.data <- data
                        colnames(tmp.data) <- dataColsVarI

                        selCols <- which(colnames(tmp.data) %in% conditions)
                        selRows <- which(rownames(tmp.data) %in% genes)
                        print("In expression boxplot $@#@#!!")
                        print(colnames(tmp.data))
                        print(conditions)
                        print(selCols)
                        print(selRows)
                        melted.tmp.data <- reshape2::melt(tmp.data[selRows, selCols])
                        if(ncol(melted.tmp.data)==1){
                                melted.tmp.data <- data.frame(Var1=rownames(tmp.data)[selRows], Var2=colnames(tmp.data)[selCols], value=melted.tmp.data, stringsAsFactors=F)
                        }
                        p <- ggplot(melted.tmp.data, aes(x=Var2, y=value, fill=Var2)) +
                        geom_boxplot() +
                        geom_jitter(aes(Var2,value),
                                position=position_jitter(width=0.1,height=0),
                                alpha=0.6,
                                size=3,
                                show.legend=FALSE
                        ) +
                        facet_grid(.~Var1, scales="free") +
                        labs(x="Condition", y="Expression Values") +
                        scale_fill_discrete(name="Condition") +
                        theme(panel.grid.minor=element_blank(), panel.grid.major.y=element_blank(), panel.grid.major.x=element_blank(), panel.background=element_rect(fill='#f5f5f5', color='#e8e8e8'), axis.text.x=element_text(angle=90, hjust=1))
                        p
                }, height=900, width=900)

                observeEvent(input$expHeatmap_submit, {
                        gVars$heatComps <- input$heatComps
                        gVars$conditionsHeat <- input$conditionsHeat
                        gVars$percDE <- input$percDE
                })

                output$expressionHeatmap <- renderPlot({
                        shiny::validate(need(!is.null(gVars$deg.list), "Waiting for Differential Analysis Results..."))
                        print("DEG LIST:")
                        print(gVars$deg.list)
                        print("DEG LIST STR:")
                        print(str(gVars$deg.list))

                        comps <- gVars$heatComps
                        conditions <- gVars$conditionsHeat

                        shiny::validate(need(!is.null(comps), "Waiting for comparisons input..."))
                        shiny::validate(need(!is.null(conditions), "Waiting for conditions input..."))

                        deg.list <- gVars$deg.list
                        data <- gVars$expr.data
                        phTable <- gVars$phTable

                        varI <- gVars$varI
                        percDE <- gVars$percDE

                        uq <- (100-(percDE/2))/100
                        lq <- (percDE/2)/100
                        print(paste0("LQ: ", lq, ", UQ: ", uq))

                        deList <- list()
                        for(comp in comps){
                                deg.list.comp <- deg.list[[comp]]
                                qVals <- quantile(deg.list.comp$score, probs=c(lq, uq))
                                lc <- qVals[[1]]
                                uc <- qVals[[2]]
                                print(paste0("COMP : ", comp))
                                print(paste0("LC: ", lc, ", UC: ", uc))
                                deList[[comp]] <- rownames(deg.list.comp)[which(deg.list.comp$score <= lc | deg.list.comp$score >= uc)]
                        }
                        deGenes <- unique(unlist(deList))
                        print("head(deGenes)")
                        print(head(deGenes))
                        shiny::validate(need(length(deGenes)>0, "No Differential Features to Plot!"))

                        dataColsVarI <- phTable[,varI]
                        print("dataColsVarI : ")
                        print(dataColsVarI)

                        tmp.data <- data
                        colnames(tmp.data) <- phTable[, gVars$sampleColName]

                        selCols <- which(dataColsVarI %in% conditions)
                        selRows <- which(rownames(tmp.data) %in% deGenes)
                        print("In expression Heatmap. Pogchamp!!")
                        print(colnames(tmp.data))
                        print(conditions)
                        print(length(selCols))
                        print(length(selRows))

                        tmp.data.sel <- tmp.data[selRows, selCols]
                        print(dim(tmp.data.sel))
                        dataColsVarI.sel <- dataColsVarI[which(dataColsVarI %in% conditions)]
                        print("dataColsVarI.sel : ")
                        print(dataColsVarI.sel)
                        classPalette <- get_color_palette(iVec=sort(conditions), asFactor=FALSE)
                        gVars$classPalette <- classPalette

                        made4::heatplot(tmp.data.sel,
                                dend="column",
                                scale="none",
                                cols.default=F,
                                lowcol="yellow",
                                highcol="red",
                                cexCol=1,
                                margins=c(8,16),
                                labRow=rownames(tmp.data.sel),
                                classvec=dataColsVarI.sel,
                                classvecCol=classPalette
                        )
                        par(lend = 1)           # square line ends for the color legend
                        legend("bottomleft",      # location of the legend on the heatmap plot
                                legend = sort(conditions), # category labels
                                col = classPalette,
                                bty = "n",
                                bg = "white",
                                border = "white",
                                cex = .7,
                                lty= 1,            # line style
                                lwd = 10           # line width
                        )
                }, height=900, width=900)

		output$phCounts <- renderPrint({
			phRows <- input$ph_rows_selected
			if (length(phRows)) {
				cat(phRows, sep = ', ')
			}
		})

                output$htmlAnnMaps <- renderUI({
			if(is.null(gVars$annLoaded)){
				dispStr <- "<p 'style color:red'><b>Waiting for upload!</b></p>"
			}else{
                                dispStr <- gVars$annMapDispStr
                        }
			HTML(dispStr)
		})

                output$totalSampleBox <- shinydashboard::renderInfoBox({
			if(is.null(gVars$phTable)){
				totalCount <- "NA"
			}else{
				totalCount <- gVars$totalSamples
			}

			shinydashboard::infoBox(
				"Total Samples", totalCount, icon=icon("list")
			)
		})

                output$filteredSampleBox <- shinydashboard::renderInfoBox({
			if(is.null(gVars$phTable)){
				filteredCount <- "NA"
			}else{
				filteredCount <- gVars$filteredSamples
			}

			shinydashboard::infoBox(
				"Remaining Samples", filteredCount, icon=icon("check"), color="green"
			)
		})

                output$removedSampleBox <- shinydashboard::renderInfoBox({
			if(is.null(gVars$phTable)){
				removedCount <- "NA"
			}else{
				removedCount <- gVars$removedSamples
			}

			shinydashboard::infoBox(
				"Samples Removed", removedCount, icon=icon("close"), color="red"
			)
		})

		output$varIValueBox <- shinydashboard::renderValueBox({
			if(is.null(input$varI)){
				varI <- "NA"
			}else{
				varI <- input$varI
			}

			shinydashboard::valueBox(
				varI, "Variable of Interest", color="blue", icon=NULL
			)
		})

		output$coVarValueBox <- shinydashboard::renderValueBox({
			if(is.null(input$coVar)){
				coVar <- "NA"
			}else{
				coVar <- paste0(input$coVar, collapse=", ")
			}

			shinydashboard::valueBox(
				coVar, "Co-Variates", color="light-blue", icon=NULL
			)
		})

		output$batchValueBox <- shinydashboard::renderValueBox({
			if(is.null(input$batch)){
				batch <- "NA"
			}else{
				batch <- paste0(input$batch, collapse=", ")
			}

			shinydashboard::valueBox(
				batch, "Batches", color="orange", icon=NULL
			)
		})

		observeEvent(input$generateQC, {
			#Progress bar callback function
			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if(is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
			
			#Start loading screen
			shinyjs::html(id="loadingText", "CREATING QC REPORT")
			shinyjs::show(id="loading-content")

			on.exit({
				progress$close()
				shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")
			})

			if(input$arrType=="af_exp"){
				exn <- ".pdf"
				con <- tempfile(pattern="eUTOPIA_Affymetrix_QC_Report_", fileext=exn)

				#Generate affy_QC_report
				phTable <- gVars$phTable
				fileNameColName <- gVars$fileNameColName
				sampleColName <- gVars$sampleColName
				#celDir <- gVars$celDir()
				rawDir <- gVars$rawDir
				affCDF <- input$affCDF
				affyPool <- as.integer(input$affyPool)
				cores <- NULL
				chkParallel <- input$chkParallel
				if(isFALSE(is.null(chkParallel)) && isTRUE(chkParallel)){
					cores <- as.integer(input$cores)
				}

				message("Running affy_QC_report...")
				affy_QC_report(fileNamesCol=fileNameColName, samplesCol=sampleColName, phTable=phTable, celDir=rawDir, nCores=cores, cdfname=affCDF, outputFile=con, isParallel=chkParallel, minFiles=affyPool, updateProgress=updateProgress)
				message("After affy_QC_report!")

				if(!gVars$QC_passed){
					shinyBS::updateCollapse(session, "bsSidebar", open="BATCH CORRECTION", style=list("QUALITY CONTROL"="success", "BATCH CORRECTION"="warning"))
				}
			}else if(input$arrType=="ag_exp2"){
				exn <- ".tar.gz"
				con <- tempfile(pattern="eUTOPIA_Agilent_QC_Report_", fileext=exn)

				tempReportDir <- file.path(tempdir(), Sys.Date())
				arrayQualityMetrics::arrayQualityMetrics(expressionset=gVars$rgList, outdir=tempReportDir, force=TRUE, reporttitle=paste0("Array Quality Metrics"))
				tar(tarfile=con, files=tempReportDir, compression="gzip", tar="tar")
				if(!gVars$QC_passed){
					shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
				}
			}else if(input$arrType=="ag_exp1"){
				exn <- ".tar.gz"
				con <- tempfile(pattern="eUTOPIA_Agilent_QC_Report_", fileext=exn)

				eset <- Biobase::ExpressionSet(assayData=Biobase::assayDataNew(exprs=gVars$eListRaw$E))
				tempReportDir <- file.path(tempdir(), Sys.Date())
				arrayQualityMetrics::arrayQualityMetrics(expressionset=eset, outdir=tempReportDir, force=TRUE, reporttitle=paste0("Array Quality Metrics"))
				tar(tarfile=con, files=tempReportDir, compression="gzip", tar="tar")
				if(!gVars$QC_passed){
					shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
				}
			}else if(input$arrType=="il_methyl"){
				exn <- ".html"
				con <- tempfile(pattern="eUTOPIA_Illumina_Methylation_QC_Report_", fileext=exn)

				if(!gVars$QC_passed){
					shinyBS::updateCollapse(session, "bsSidebar", open="PROBE FILTERING", style=list("QUALITY CONTROL"="success", "PROBE FILTERING"="warning"))
				}
			}
			gVars$QC_passed <- TRUE
			gVars$QC_con <- con

			shinyBS::toggleModal(session, "qcReportModal", toggle="open")
		})

                output$exportQC <- shiny::downloadHandler(
			filename = function(){
                                if(input$arrType=="af_exp"){
					exn <- ".pdf"
                                        resStr <- paste0("eUTOPIA_Affymetrix_QC_Report_", Sys.Date(), exn)
                                }else if(grepl("ag_exp.*", input$arrType)){
                                        exn <- ".tar.gz"
                                        resStr <- paste0("eUTOPIA_Agilent_QC_Report_", Sys.Date(), exn)
                                }else if(input$arrType=="il_methyl"){
                                        exn <- ".html"
                                        resStr <- paste0("eUTOPIA_Illumina_Methylation_QC_Report_", Sys.Date(), exn)
                                }
                                return(resStr)
			},
			content = function(con){
				if(gVars$QC_passed==TRUE){
					file.copy(gVars$QC_con, con)
				}
			}
		)

		##DYNAMIC UI WIDGETS

                output$loading_text <- renderText({
                        loadingText <- "LOADING"
                        return(loadingText)
                })

		output$selFileNameCol <- renderUI({
			selectInput("fileNameCol", "Filename Variable", choices=gVars$phColChoices())
		})

		output$selDyeCol <- renderUI({
			selectInput("dyeCol", "Dye Variable", choices=gVars$phColChoices())
		})

		output$selSampleIDCol <- renderUI({
			selectInput("sampleIDCol", "Sample ID Variable", choices=gVars$phColChoices())
		})

		output$selSep <- renderUI({
			selectInput("sepS", "Field Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
		})

		output$selQuote <- renderUI({
			selectInput("quote", "Quotes", choices=gVars$quoteChoices, selected=gVars$quoteChoices[1])
		})

                output$selSepAnno <- renderUI({
			selectInput("sepSAnno", "Field Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
		})

		output$selAffCDF <- renderUI({
			selectInput("affCDF", "Available Annotation CDF", choices=gVars$cdfChoices(), selected=gVars$cdfChoices()[1])
		})

		output$selVarI <- renderUI({
			selectInput("varI", label="Variable of Interest", choices=gVars$pcChoices())
		})

		output$selCoVar <- renderUI({
			selectInput("coVar", label="Co-Variates", choices=gVars$pcChoices(), multiple=TRUE)
		})

		output$selBatch <- renderUI({
			selectInput("batch", label="Batches", choices=gVars$pcChoicesSV(), multiple=TRUE)
		})

		output$selVarISva <- renderUI({
			selectInput("varISva", label="Variable of Interest", choices=gVars$pcChoices())
		})

		output$selCoVarSva <- renderUI({
			selectInput("coVarSva", label="Co-Variates", choices=gVars$pcChoices(), multiple=TRUE)
		})

		output$selBatchSva <- renderUI({
			selectInput("batchSva", label="Batches", choices=gVars$pcChoices(), multiple=TRUE)
		})

                output$selVarICombat <- renderUI({
			selectInput("varICombat", label="Variable of Interest", choices=gVars$pcChoices())
		})

		output$selCoVarCombat <- renderUI({
			selectInput("coVarCombat", label="Co-Variates", choices=gVars$pcChoices(), multiple=TRUE)
		})

		output$selBatchCombat <- renderUI({
                        tmpChoices <- gVars$pcChoices()
                        if(input$corrType=="sc"){
                                if(gVars$svaStep==1){
                                        tmpChoices <- names(gVars$svaSV)
                                }
                        }
			selectInput("batchCombat", label="Batches", choices=tmpChoices, multiple=TRUE)
		})

                output$selVarILimma <- renderUI({
                        phTable <- gVars$phTable
			if(is.null(phTable)){
                                tmpChoices <- c("NA")
                                sel <- tmpChoices[1]
                        }else{
                                tmpChoices <- colnames(phTable)
                                if(!input$chkKeepNA){
                                        tmpChoices <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
                                        if(length(tmpChoices)==0){
                                                shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
                                                tmpChoices <- c("NA")
                                        }
                                }
                                if(input$corrType=="sc" || input$corrType=="c"){
                                        sel <- gVars$varICombat
                                }else if(input$corrType=="s"){
                                        sel <- gVars$varISva
                                }else if(input$corrType=="n"){
                                        #sel <- gVars$pcChoices()[1]
                                        sel <- tmpChoices[1]
                                }
                        }
			selectInput("varILimma", label="Variable of Interest", choices=tmpChoices, selected=sel)
		})

		output$selCoVarLimma <- renderUI({
                        phTable <- gVars$phTable

			if(is.null(phTable)){
                                tmpChoices <- c("NA")
                                sel <- tmpChoices[1]
                        }else{
                                if(!is.null(gVars$svaSV)){
                                        svaSV <- as.data.frame(apply(gVars$svaSV, 2, factor))
                                        svaSVc <- gVars$svaSVc
                                        phTable <- cbind(phTable, svaSV, svaSVc)
                                }

                                tmpChoices <- colnames(phTable)
                                if(!input$chkKeepNA){
                                        tmpChoices <- colnames(phTable)[apply(phTable, 2, function(x){res<-any(is.na(x));!res})]
                                        if(length(tmpChoices)==0){
                                                shinyjs::info("No available phenotype variables.\nAll variables contain NAs\nCheck 'Keep NA Variables' to include variables with NA values")
                                                tmpChoices <- c("NA")
                                        }
                                }
                                if(input$corrType=="sc" || input$corrType=="c"){
                                        sel <- gVars$batchCombat
                                }else if(input$corrType=="s"){
                                        sel <- names(gVars$svaSVc)
                                }else if(input$corrType=="n"){
                                        sel <- NULL
                                }
                        }
			selectInput("coVarLimma", label="Co-Variates", choices=tmpChoices, multiple=TRUE, selected=sel)
		})

		output$selMdsLabel <- renderUI({
                        if(is.null(gVars$varICombat)){
                                selIdx <- 1
                        }else{
                                selIdx <- which(gVars$pcChoices() %in% gVars$varICombat)
                        }
			selectInput("mdsLabel", "MDS Label Name", choices=gVars$pcChoices(), selected=gVars$pcChoices()[selIdx])
		})

		output$selMdsColor <- renderUI({
                        if(is.null(gVars$varICombat)){
                                selIdx <- 1
                        }else{
                                selIdx <- which(gVars$pcChoices() %in% gVars$varICombat)
                        }
			selectInput("mdsColor", "MDS Label Color", choices=gVars$pcChoices(), selected=gVars$pcChoices()[selIdx])
		})

		output$selVarIBoxPlot <- renderUI({
                        if(input$corrType=="sc" || input$corrType=="c"){
                                sel <- gVars$varICombat
                        }else if(input$corrType=="s"){
                                sel <- gVars$varISva
                        }else if(input$corrType=="n"){
                                sel <- gVars$pcChoices()[1]
                        }
			selectInput("varIBoxPlot", label="Variable of Interest", choices=gVars$pcChoices(), selected=sel)
		})

		output$selNormMethod <- renderUI({
                        if(input$arrType=="il_methyl"){
			        selectInput("normMethod", "Normalization Method",
                                        choices=c(
                                                "SWAN"="SWAN",
                                                "Raw"="Raw",
                                                "Quantile"="Quantile",
                                                "Normal-exponent out-of-band"="Noob",
                                                "Illumina"="Illumina",
                                                "Functional Normalization"="FunNorm"
                                        ),
                                        selected="SWAN",
                                        multiple=FALSE
                                )
                        }else{
			        selectInput("normMethod", "Normalization Type", choices=gVars$normChoices, selected="BA", multiple=FALSE)
                        }
		})

		output$selNormMethod2 <- renderUI({
			selectInput("normMethod2", "Method", choices=gVars$baChoices, multiple=FALSE)
		})

		output$selID <- renderUI({
			selectInput("ID", "Select ID", choices=gVars$idChoices()[-1], multiple=FALSE)
		})

		output$selTreatment <- renderUI({
			selectInput("treatment", "Condition 1", choices=gVars$condChoices(), multiple=FALSE)
		})

		output$selControl <- renderUI({
			selectInput("control", "Condition 2", choices=gVars$condChoices(), multiple=FALSE)
		})

		output$selComps <- renderUI({
                        compsSel <- gVars$compsSel
                        if(length(compsSel)==0){
                                compsSel <- "NA"
                        }
			selectInput("comps", "Comparisons", choices=gVars$compChoices(), multiple=TRUE, selected=compsSel)
		})

		output$selIntersectComps <- renderUI({
                        if(is.null(gVars$deComps)){
                                deChoices <- c("NA")
                        }else{
                                deChoices <- gVars$deComps
                        }
			selectInput("intersectComps", "Sets to Intersect", choices=deChoices, multiple=TRUE, selected=deChoices)
		})

		output$selPvAdjMethod <- renderUI({
			selectInput("pvAdjMethod", "P.Value Adjustment Method", choices=gVars$pvAdjChoices, multiple=FALSE, selected="none")
		})

		output$selCompDeTable <- renderUI({
                        if(is.null(gVars$deComps)){
                                deChoices <- c("NA")
                        }else{
                                deChoices <- gVars$deComps
                        }
			selectInput("compDE", "Select Comparison", choices=deChoices, multiple=FALSE)
		})

		output$selConditionsBoxplot <- renderUI({
                        phTable <- gVars$phTable
                        if(is.null(phTable)){
                                tmpChoices <- c("NA")
                        }else{
                                rownames(phTable) <- phTable[,gVars$sampleColName]
                                phTable <- as.data.frame(apply(phTable, 2, factor))
                                varI <- input$varIBoxPlot
                                conditions <- levels(phTable[,varI])
                                if(is.null(conditions)){
                                        tmpChoices <- c("NA")
                                }else{
                                        tmpChoices <- conditions
                                }
                        }
			selectInput("conditionsBoxplot", "Select Conditions", choices=tmpChoices, multiple=TRUE, selected=tmpChoices)
		})

		output$selBoxplotComp <- renderUI({
                        if(is.null(gVars$deComps)){
                                tmpChoices <- c("NA")
                        }else{
                                tmpChoices <- gVars$deComps
                        }
			selectInput("boxplotComp", "Select Comparison", choices=tmpChoices, multiple=FALSE, selected=tmpChoices[1])
		})

		output$selConditionsHeat <- renderUI({
                        if(is.null(gVars$conditions)){
                                tmpChoices <- c("NA")
                        }else{
                                tmpChoices <- gVars$conditions
                        }
			selectInput("conditionsHeat", "Select Conditions", choices=tmpChoices, multiple=TRUE, selected=tmpChoices)
		})

                output$selHeatComps <- renderUI({
                        if(is.null(gVars$deComps)){
                                deChoices <- c("NA")
                        }else{
                                deChoices <- gVars$deComps
                        }
			selectInput("heatComps", "For Comparison(s)", choices=deChoices, multiple=TRUE, selected=deChoices[1])
		})

		output$selCores <- shiny::renderUI({
                        shinyjs::disabled(shiny::selectInput("cores", "No. of Cores", choices=(1:detectCores()), multiple=FALSE, selected=2))
                })

                output$selPoolSize <- shiny::renderUI({
                        if(gVars$loadedRaw){
                                phTable <- gVars$phTable
                                fileNameColName <- gVars$fileNameColName
                                fileNames <- unique(phTable[,fileNameColName])
                                min <- min(length(fileNames), 20)
                                max <- min(length(fileNames), 50)
                                poolChoices <- c(min:max)
                                poolSel <- min
                        }else{
                                poolChoices <- c("NA")
                                poolSel <- c("NA")
                        }
                        shiny::selectInput("affyPool", "Pool Size of Raw Data", choices=poolChoices, multiple=FALSE, selected=poolSel)
                })

		observe({
                        disable_ann <- FALSE
                        disable_de <- FALSE
                        disable_sva <- FALSE
                        disable_combat <- FALSE
                        disable_corr_skip <- FALSE
                        disable_limma <- FALSE

			if(is.null(input$fPheno)){
				shinyjs::disable("load_pheno_submit")
                        }else{
				shinyjs::enable("load_pheno_submit")
                        }

			if(is.null(gVars$phTable)){
				shinyjs::disable("dirButton")
			}else{
				shinyjs::enable("dirButton")
			}

                        if(is.null(gVars$phLoaded)){
			        shinyjs::hide("phenoPreviewDiv")
                        }else{
				shinyjs::show("phenoPreviewDiv")
                        }

                        if(is.null(gVars$loadedRaw) || !isTRUE(gVars$loadedRaw)){
                                print("Disable QC buttons...")
			        shinyjs::disable("exportQC")
			        shinyjs::disable("qc_methyl_submit")
			        shinyjs::disable("qc_skip_submit")
                        }else{
                                print("Enable QC buttons...")
				shinyjs::enable("exportQC")
			        shinyjs::enable("qc_methyl_submit")
			        shinyjs::enable("qc_skip_submit")
                        }

                        if(is.null(gVars$normalized) || !isTRUE(gVars$normalized)){
			        shinyjs::disable("exportRpt")
                        }else{
				shinyjs::enable("exportRpt")
                        }

                        if(is.null(gVars$filt.data)){
			        shinyjs::disable("exportFiltMat")
				shinyjs::disable("norm_submit")
                        }else{
				shinyjs::enable("exportFiltMat")
				shinyjs::enable("norm_submit")
                        }

                        if(is.null(gVars$norm.data)){
			        shinyjs::disable("exportNormMat")
                        }else{
				shinyjs::enable("exportNormMat")
                        }

                        if(is.null(gVars$comb.data)){
			        shinyjs::disable("exportCorrMat")
                        }else{
				shinyjs::enable("exportCorrMat")
                        }

                        if(is.null(gVars$agg.data)){
			        shinyjs::disable("exportAggMat")
                        }else{
				shinyjs::enable("exportAggMat")
                        }

			if(input$arrType=="ag_exp2"){
				shinyjs::show("selDyeCol")
			}else if(input$arrType=="af_exp"){
				if(is.null(input$affCDF) || is.na(input$affCDF)){
					shinyjs::disable("upload_raw_submit")
				}else{
					shinyjs::enable("upload_raw_submit")
				}
				shinyjs::hide("selDyeCol")
                        }else if(input$arrType=="ag_exp1"){
				shinyjs::hide("selDyeCol")
                        }else if(input$arrType=="il_methyl"){
				shinyjs::hide("selDyeCol")
                        }

			if(is.null(gVars$rawDir)){
				shinyjs::disable("upload_raw_submit")
			}else{
				shinyjs::enable("upload_raw_submit")
			}

			if(is.null(input$cdfURL) || input$cdfURL==""){
				shinyjs::disable("install_cdf_submit")
				shinyBS::addTooltip(session, id="install_cdf_submit", title="Please provide a valid URL for the CDF file!", placement="bottom")
			}else{
				shinyjs::enable("install_cdf_submit")
				shinyBS::removeTooltip(session, id="install_cdf_submit")
			}

                        if(gVars$QC_passed){
                                shinyjs::enable("filt_submit")
                                shinyjs::enable("filtDist")
                                shinyjs::enable("perSamples")
                        }else{
                                shinyjs::disable("filt_submit")
                                shinyjs::disable("filtDist")
                                shinyjs::disable("perSamples")
                                if(input$arrType=="af_exp"){
                                        #Disable sva
                                        disable_sva <- TRUE

                                        #Disable combat
                                        disable_combat <- TRUE

                                        #Disable corr_skip
                                        disable_corr_skip <- TRUE
                                }
                        }

                        if(input$annType=="raw"){
			        shinyjs::hide("launchAnnModalDiv")
				shinyjs::show("submitAnnDiv")
                        }else if(input$annType=="file"){
			        shinyjs::hide("submitAnnDiv")
				shinyjs::show("launchAnnModalDiv")
                        }

			if(!is.null(input$sepS) && input$sepS=="OTHER"){
				shinyjs::enable("sepT")
			}else{
				shinyjs::disable("sepT")
			}

			if(!is.null(input$sepSAnno) && input$sepSAnno=="OTHER"){
				shinyjs::enable("sepTAnno")
			}else{
				shinyjs::disable("sepTAnno")
			}

                        if(!is.null(input$fileNameCol) && !is.null(input$sampleIDCol)){
                                disable <- FALSE
				if(input$fileNameCol==input$sampleIDCol){
                                        disable <- TRUE
				}else{
                                        if(!is.null(input$dyeCol)){
                                                if(input$fileNameCol==input$dyeCol || input$sampleIDCol==input$dyeCol){
                                                        disable <- TRUE
                                                }
                                        }
				}

                                if(disable){
					shinyjs::disable("upload_pheno_submit")
                                        shinyBS::addTooltip(session, id="upload_pheno_submit", title="Specified variables cannot be represented by the same column", placement="bottom")
                                }else{
					shinyjs::enable("upload_pheno_submit")
                                        shinyBS::removeTooltip(session, id="upload_pheno_submit")
                                }
			}

                        if(gVars$normalized){
                                if(input$corrType=="sc"){
                                        #if(input$arrType=="ag_exp2"){
					if(grepl("ag_exp.*", input$arrType)){
                                                if(is.null(gVars$comb.data)){
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Please peform batch correction before annotation import!", placement="bottom")
                                                        disable_ann <- TRUE
                                                        disable_de <- TRUE
                                                }else{
                                                        shinyBS::removeTooltip(session, id="launch_ann_modal")
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Launch a graphical window, to configure and import annotation!", placement="bottom")
                                                        if(is.null(gVars$agg.data)){
                                                                disable_de <- TRUE
                                                        }else{
                                                                if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                        disable_de <- TRUE
                                                                }else{
                                                                }
                                                        }
                                                }
                                        }else{
                                                if(is.null(gVars$comb.data)){
                                                        disable_de <- TRUE
                                                }else{
                                                        if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                disable_de <- TRUE
                                                        }else{
                                                        }
                                                }
                                        }

                                        if(is.null(gVars$svaStep)){
                                                disable_combat <- TRUE
                                        }else{
                                                shinyBS::addTooltip(session, id="launch_combat_modal", title="Launch a graphical window, to configure and run batch correction with ComBat!", placement="bottom")
                                        }
                                }else if(input$corrType=="s"){
                                        #if(input$arrType=="ag_exp2"){
					if(grepl("ag_exp.*", input$arrType)){
                                                if(is.null(gVars$svaStep)){
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Please peform batch correction before annotation import!", placement="bottom")
                                                        disable_ann <- TRUE
                                                        disable_de <- TRUE
                                                }else{
                                                        shinyBS::removeTooltip(session, id="launch_ann_modal")
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Launch a graphical window, to configure and import annotation!", placement="bottom")
                                                        if(is.null(gVars$agg.data)){
                                                                disable_de <- TRUE
                                                        }else{
                                                                if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                        disable_de <- TRUE
                                                                }else{
                                                                }
                                                        }
                                                }
                                        }else{
                                                if(is.null(gVars$svaStep)){
                                                        disable_de <- TRUE
                                                }else{
                                                        if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                disable_de <- TRUE
                                                        }else{
                                                        }
                                                }
                                        }
                                }else if(input$corrType=="c"){
                                        #if(input$arrType=="ag_exp2"){
					if(grepl("ag_exp.*", input$arrType)){
                                                if(is.null(gVars$comb.data)){
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Please peform batch correction before annotation import!", placement="bottom")
                                                        disable_ann <- TRUE
                                                        disable_de <- TRUE
                                                }else{
                                                        shinyBS::removeTooltip(session, id="launch_ann_modal")
                                                        shinyBS::addTooltip(session, id="launch_ann_modal", title="Launch a graphical window, to configure and import annotation!", placement="bottom")
                                                        if(is.null(gVars$agg.data)){
                                                                disable_de <- TRUE
                                                        }else{
                                                                if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                        disable_de <- TRUE
                                                                }else{
                                                                }
                                                        }
                                                }
                                        }else{
                                                if(is.null(gVars$comb.data)){
                                                        disable_de <- TRUE
                                                }else{
                                                        if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                disable_de <- TRUE
                                                        }else{
                                                        }
                                                }
                                        }
                                        shinyBS::addTooltip(session, id="launch_combat_modal", title="Launch a graphical window, to configure and run batch correction with ComBat!", placement="bottom")
                                        shinyjs::enable("launch_combat_modal")
                                }else if(input$corrType=="n"){
                                        #if(input$arrType=="ag_exp2"){
					if(grepl("ag_exp.*", input$arrType)){
                                                shinyBS::removeTooltip(session, id="launch_ann_modal")
                                                shinyBS::addTooltip(session, id="launch_ann_modal", title="Launch a graphical window, to configure and import annotation!", placement="bottom")

                                                if(!gVars$corrected){
                                                        disable_ann <- TRUE
                                                }

                                                if(is.null(gVars$agg.data)){
                                                        disable_de <- TRUE
                                                }else{
                                                        if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                                disable_de <- TRUE
                                                        }else{
                                                        }
                                                }
                                        }else{
                                                if(is.null(input$comps[1]) || is.na(input$comps[1]) || input$comps[1]=="NA"){
                                                        disable_de <- TRUE
                                                }else{
                                                }
                                        }
                                }
                                
                        }else{
                                disable_ann <- TRUE

                                #Disable differential
                                disable_de <- TRUE

                                #Disable sva
                                disable_sva <- TRUE

                                #Disable combat
                                disable_combat <- TRUE

                                #Disable corr_skip
                                disable_corr_skip <- TRUE
                        }

			if(!is.null(input$treatment) && !is.null(input$control)){
				if(input$treatment==input$control){
					shinyjs::disable("add_comp_submit")
				}else{
					shinyjs::enable("add_comp_submit")
				}
			}

                        if(is.null(gVars$deg.list)){
                                shinyjs::disable("exportDE")
                        }else{
                                shinyjs::enable("exportDE")
                        }

                        if(is.null(input$batchCombat)){
                                shinyjs::disable("combat_submit")
                        }else{
                                shinyjs::enable("combat_submit")
                        }

                        #Disable or enable condition resolution
                        if(disable_ann){
                                shinyjs::disable("launch_ann_modal")
                                shinyjs::disable("submit_ann")
                        }else{
                                shinyjs::enable("launch_ann_modal")
                                shinyjs::enable("submit_ann")
                        }

                        if(disable_de){
                                shinyjs::disable("de_submit")
                        }else{
                                shinyjs::enable("de_submit")
                        }

                        if(disable_sva){
                                shinyjs::disable("launch_sva_modal")
                        }else{
                                shinyjs::enable("launch_sva_modal")
                        }

                        if(disable_combat){
                                shinyjs::disable("launch_combat_modal")
                        }else{
                                shinyjs::enable("launch_combat_modal")
                        }

                        if(disable_corr_skip){
                                shinyjs::disable("corr_skip_submit")
                        }else{
                                shinyjs::enable("corr_skip_submit")
                        }

			if(is.null(input$chkParallel) || isFALSE(input$chkParallel)){
                                print("Disable selCores!")
                                shinyjs::disable("selCores")
                        }else{
                                print("Enable selCores!")
                                print("input$chkParallel")
                                print(input$chkParallel)
                                shinyjs::enable("selCores")
                        }
		})

		shinyBS::addTooltip(session, id="selVarI", title="Variable (Column of data) from the phenotype file that represents the biological variation of primary interest!", placement="top", trigger="focus")
		shinyBS::addTooltip(session, id="selCoVar", title="Variables that represents the biological variation not associated with the primary variable of interest (Should not be confounded with other variables)!", placement="top", trigger="focus")
		shinyBS::addTooltip(session, id="selBatch", title="Variables that represents the technical variation, such as 'slide', 'dye', 'area', 'array', etc. (Should not be confounded with other variables)!", placement="top", trigger="focus")
		shinyBS::addTooltip(session, id="selTreatment", title="Select a treated condition from the 'Variable of Interest'!", placement="top", trigger="focus")
		shinyBS::addTooltip(session, id="selControl", title="Select the corresponding control condition from the 'Variable of Interest'!", placement="top", trigger="focus")
		shinyBS::addTooltip(session, id="selComps", title="Add/Remove the composed comparisons for differential analysis!", placement="top", trigger="focus")

		##Hide the loading message when the rest of the server function has executed
		Sys.sleep(1)
		shinyjs::hide(id="loading-content", anim=TRUE, animType="fade")

		##Check the array type from the main screen and adjust widgets accordingly
		observeEvent(input$submit_array_type, {
			shinyjs::hide(id="sel-array-type", anim=TRUE, animType="fade")
                        arrType <- input$arrType
			print(paste0("Array Type: ", arrType))
                        if(arrType=="ag_exp2"){
				shinyjs::hide("affAnnDiv")
				shinyjs::hide("detectPV")
				shinyjs::show("agSourceDiv")
                        }else if(arrType=="af_exp"){
				shinyjs::show("affAnnDiv")
				shinyjs::hide("launch_ann_modal")
				shinyjs::hide("detectPV")
				shinyjs::hide("agSourceDiv")

				#Test bsCollapsePanel hiding function
				hideBSCollapsePanel(session, panel.name="PROBE FILTERING")
				hideBSCollapsePanel(session, panel.name="NORMALIZATION")
				hideBSCollapsePanel(session, panel.name="ANNOTATION")

				#For affymetyrix QC with Affymetrix outlier QC script
				shinyjs::show("affy_QC_div")
                        }else if(arrType=="ag_exp1"){
				shinyjs::hide("affAnnDiv")
				shinyjs::hide("detectPV")
				shinyjs::show("agSourceDiv")
                        }else if(arrType=="il_methyl"){
				shinyjs::show("detectPV")
				shinyjs::hide("selNormMethod2")
				shinyjs::hide("filtDist")
				shinyjs::hide("affAnnDiv")
				shinyjs::hide("launch_ann_modal")
				shinyjs::hide("agSourceDiv")

                                #For methylation QC with shinyMethyl
                                shinyjs::hide("download_QC_div")
                                shinyjs::show("methyl_QC_div")
                                shinyjs::show("ilChkDiv")

				hideBSCollapsePanel(session, panel.name="ANNOTATION")
                        }
		})

		observeEvent(input$corrType, {
                        corrType <- input$corrType
                        if(is.null(corrType))
                        return(NULL)

                        if(corrType=="c"){
                                shinyjs::show("combatDiv")
                                shinyjs::hide("svaDiv")
                                shinyjs::hide("skipDiv")
                        }else if(corrType=="sc"){
                                shinyjs::show("combatDiv")
                                shinyjs::show("svaDiv")
                                shinyjs::hide("skipDiv")
                        }else if(corrType=="s"){
                                shinyjs::hide("combatDiv")
                                shinyjs::show("svaDiv")
                                shinyjs::hide("skipDiv")
                        }else if(corrType=="n"){
                                shinyjs::hide("combatDiv")
                                shinyjs::hide("svaDiv")
                                shinyjs::show("skipDiv")
                        }
		})
	}
)
