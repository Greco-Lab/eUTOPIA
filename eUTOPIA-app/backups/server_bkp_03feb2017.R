library(shiny)
library(DT)
library(shinyjs)
library(shinyBS)
library(shinydashboard)
library(shinyFiles)
library(limma)
library(sva)
library(swamp)
library(ggplot2)
library(biomaRt)
library(infotheo)
source("/home/veer/git/AROMA/agilent_twocolors_array_analysis_pipe_fixed.R")

options(shiny.maxRequestSize=500*1024^2)

phRows <- NULL
celTable <- NULL
phTable <- NULL
extractedList <- NULL
pcChoices <- list("RIN", "array", "area", "slide", "dye", "group")

shinyServer(
	function(input, output, session){		
                gVars <- shiny::reactiveValues(phTable=NULL, rgList=NULL, celDir=NULL, totalSamples=NULL, filteredSamples=NULL, removedSamples=NULL, norm.data=NULL, celDir=NULL, pcChoices=NULL, comb.data=NULL, agg.data=NULL, comps=list())

		gVars$sepChoices <- c("TAB", ",", ";", "SPACE", "OTHER")
		gVars$pvAdjChoices <- c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
		gVars$normChoices <- c("quantile", "within arrays", "between arrays", "median absolute values", "robust spline", "cyclic loess")
		
		phRows <- reactive(as.numeric(input$pheno_rows_selected))

                roots <- c(wd="/")
                shinyDirChoose(input, "dirButton", roots=roots)
                gVars$celDir <- reactive({
                        if(is.null(input$dirButton))
                        return(NULL)

                        roots <- c(wd="/")
                        celDir <- parseDirPath(roots, input$dirButton)
                        return(celDir)
                })

		gVars$annDF <- reactive({
			annFile <- input$fAnno
			if (is.null(annFile))
			return(NULL)
			
			annDF <- read.csv(annFile$datapath, header=TRUE, sep="\t", stringsAsFactors=FALSE, quote="")
			return(annDF)
		})

		output$annDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$annDF()), "No phenotype file!")
			)

                        annDF <- gVars$annDF()
			DT::datatable(annDF, filter="none", 
				options = list(
					ch = list(regex=TRUE, caseInsensitive=FALSE), 
					scrollX=TRUE, 
					pageLength=3,
					lengthMenu=c(1,2,3),
					ordering=F
				)
			)
		},server=TRUE)

		gVars$idChoices <- reactive({
                        annDF <- gVars$annDF()
			if (is.null(annDF))
                        return("Waiting for annotation file!")

			return(colnames(annDF))
		})

                output$dirText <- renderText({
                        if(is.null(gVars$celDir()))
                        return("Select Directory")

                        return(gVars$celDir())
                })

		gVars$inputPh <- eventReactive(input$load_pheno_submit, {
			if(is.null(input$fPheno))
			return(NULL)

                        phFile <- input$fPheno
			sepS <- input$sepS
			sepT <- input$sepT
			sepChar=NULL
			if(sepS=="OTHER"){
				sepChar <- sepT
			}else{
				if(sepS=="TAB"){
					sepChar="\t"
				}else if(sepS=="SPACE"){
					sepChar=" "
				}else{
					sepChar=sepS
				}
			}
			print(sepChar)

                        phTable <- read.csv(phFile$datapath, header=TRUE, sep=sepChar, stringsAsFactors=FALSE)
			return(phTable)
		})

		output$phRowsText <- renderText({
			if(is.null(gVars$inputPh())){
				nRow <- "NA"
			}else{
				nRow <- nrow(gVars$inputPh())
			}

			return(paste0("Rows: ", nRow))
		})

		output$phColsText <- renderText({
			if(is.null(gVars$inputPh())){
				nCol <- "NA"
			}else{
				nCol <- ncol(gVars$inputPh())
			}

			return(paste0("Columns: ", nCol))
		})

		gVars$phColChoices <- reactive({
			if(is.null(gVars$inputPh()))
			return("NA")

			choicesVec <- seq(1,ncol(gVars$inputPh()))
			choicesNames <- paste0("Column ", choicesVec)
			names(choicesVec) <- choicesNames
			return(choicesVec)
		})

		output$phenoDT <- DT::renderDataTable({
			shiny::validate(
				need(!is.null(gVars$inputPh()), "No phenotype file!")
			)

                        phTable <- gVars$inputPh()
			DT::datatable(phTable, filter="none", 
				options = list(
					ch = list(regex=TRUE, caseInsensitive=FALSE), 
					scrollX=TRUE, 
					pageLength=3,
					lengthMenu=c(1,2,3),
					ordering=F
				)
			)
		},server=TRUE)

                observeEvent(input$upload_pheno_submit, {
                        #phFile <- input$fPheno
			#shiny::validate(
			#	need(!is.null(phFile), "No Phenotype File Provided!")
			#)
			shiny::validate(
				need(!is.null(gVars$inputPh()), "No Phenotype File Provided!")
			)
			
                        #phTable <- read.csv(phFile$datapath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                        #dyeColID <- grep("Cy3|Cy5", phTable[1,], ignore.case=T)
                        #id_file_names <- grep("file", colnames(phTable), ignore.case = TRUE)
                        phTable <- gVars$inputPh()
			dyeColID <- as.integer(input$dyeCol)
			fileNameColID <- as.integer(input$fileNameCol)
                        phTable <- phTable[order(phTable[,dyeColID], phTable[,fileNameColID]),]
                        gVars$phTable <- phTable
			gVars$dyeColID <- dyeColID
			gVars$fileNameColID <- fileNameColID
                        gVars$totalSamples <- nrow(phTable)
                        gVars$filteredSamples <- nrow(phTable)
                        gVars$removedSamples <- 0
			gVars$removedSamplesInfo <- NULL
                })

		gVars$pcChoices <- reactive({
                        phTable <- gVars$phTable
			if (is.null(phTable))
                        return("No Options!")

			return(colnames(phTable))
		})

		gVars$condChoices <- reactive({
                        phFile <- input$fPheno
			if (is.null(phFile))
                        return("No Options!")
                        
			varI <- input$varI
			if (is.null(varI))
                        return("No Options!")

                        phTable <- read.csv(phFile$datapath, header=TRUE, sep="\t", stringsAsFactors=FALSE)
                        conds <- levels(factor(phTable[,varI]))
			return(conds)
		})

                observeEvent(input$filterPh, {
                        if(is.null(gVars$phTable) || is.null(input$filtered_rows_selected))
                        return(NULL)

			dyeColID <- gVars$dyeColID
			fileNameColID <- gVars$fileNameColID
                        phTable <- gVars$phTable
                        phRows2Remove <- as.integer(input$filtered_rows_selected)
			if (length(phRows2Remove)>0){
				phRows2Remove <- which(rownames(phTable) %in% phRows2Remove)
				print(phRows2Remove)
				removedSamplesInfo <- phTable[phRows2Remove,c(fileNameColID, dyeColID)]
				phTable <- phTable[-phRows2Remove,]
			}
                        gVars$phTable <- phTable
                        gVars$filteredSamples <- nrow(phTable)
                        gVars$removedSamples <- gVars$totalSamples - nrow(phTable)
			gVars$removedSamplesInfo <- rbind(gVars$removedSamplesInfo, removedSamplesInfo)
			#print(gVars$removedSamplesInfo)
                })

		output$filtered <- DT::renderDataTable({
                        if(is.null(gVars$phTable))
                        return(NULL)

                        phTable <- gVars$phTable
			#DT::datatable(phTable, filter = list(position='top', clear=FALSE), options = list(search = list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
			DT::datatable(phTable, filter="none", 
				options = list(
					search = list(regex=TRUE, caseInsensitive=FALSE), 
					scrollX=TRUE, 
					ordering=F
				)
			)
		},server=TRUE)

                observeEvent(input$upload_raw_submit, {
			shiny::validate(
				need(!is.null(gVars$celDir()), "No directory selected!")
			)

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())
                        progress$set(message="Raw Data:", value=0)

                        phTable <- gVars$phTable
			removedSamplesInfo <- gVars$removedSamplesInfo
			fileNameColID <- gVars$fileNameColID
			celDir <- gVars$celDir()

			smp2rmG <- NULL
			smp2rmR <- NULL
			if(!is.null(removedSamplesInfo)){
                                print("Checking samples to remove...")
				tmpVal <- table(removedSamplesInfo[,1])
				smp2rm <- names(tmpVal[tmpVal==1])
				smp2rmDF <- removedSamplesInfo[which(removedSamplesInfo[,1] %in% smp2rm),]
				smp2rmID <- grep("Cy3", smp2rmDF[,2], ignore.case=T)
				if(length(smp2rmID)>0){
					smp2rmG <- smp2rmDF[smp2rmID,1]
					smp2rmG <- gsub("\\.[a-zA-Z]*$", "", smp2rmG)
					print("Green Samples to Remove")
					print(smp2rmG)
				}
				smp2rmID <- grep("Cy5", smp2rmDF[,2], ignore.case=T)
				if(length(smp2rmID)>0){
					smp2rmR <- smp2rmDF[smp2rmID,1]
					smp2rmR <- gsub("\\.[a-zA-Z]*$", "", smp2rmR)
					print("Red Samples to Remove")
					print(smp2rmR)
				}
			}

                        updateProgress(detail="Reading CEL Files...", value=1/3)
			fileNames <- unique(phTable[,fileNameColID])
                        rgList <- read.maimages(files=fileNames, source="agilent.median", path=celDir, verbose=TRUE)

                        updateProgress(detail="Preprocessing & Filtering...", value=2/3)
			gCh <- rgList$G
			rCh <- rgList$R
			if(!is.null(smp2rmG)){
				print(str(gCh))
				colID <- which(colnames(gCh) %in% smp2rmG)
				gCh <- gCh[,-c(colID)]
				print(str(gCh))
			}
			if(!is.null(smp2rmR)){
				print(str(rCh))
				colID <- which(colnames(rCh) %in% smp2rmR)
				rCh <- rCh[,-c(colID)]
				print(str(rCh))
			}

			#data <- cbind(rgList$G, rgList$R)
			data <- cbind(gCh, rCh)
			rownames(data) <- rgList$genes$ProbeName
			colnames(data) <- rownames(phTable)
			nc.data <- data[which(rgList$genes$ControlType==0),]
			c.data  <- data[which(rgList$genes$ControlType==-1),]
                        updateProgress(detail="Completed!", value=2/3)

                        gVars$rgList <- rgList
                        gVars$nc.data <- nc.data
                        gVars$c.data <- c.data

                        shinyjs::show(id="filtOptions", anim=TRUE)
                        shinyjs::hide(id="rawOptions", anim=TRUE)
		})
		
		gVars$filt.data <- reactive({
			if(is.null(gVars$rgList))
			return(NULL)

			#rgList <- gVars$rgList
			nc.data <- gVars$nc.data
			c.data <- gVars$c.data
			qfilt <- as.numeric(input$filtDist)
			perc <- as.numeric(input$perSamples)

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())
                        progress$set(message="Filtering:", value=0)

                        updateProgress(detail="...", value=1/2)
			#data <- cbind(rgList$G, rgList$R)
			#rownames(data) <- rgList$genes$ProbeName
			#colnames(data) <- rownames(pd)
			#nc.data <- data[which(rgList$genes$ControlType==0),]
			#c.data  <- data[which(rgList$genes$ControlType==-1),]
			#gVars$nc.data <- nc.data
			filt.data <- filt.by.neg.probes(nc.data, c.data, qdist=qfilt, perc=perc, verbose=T)
                        updateProgress(detail="Completed!", value=2/2)

			return(filt.data)
		})

		output$percProbesText <- renderText({
			if(is.null(gVars$filt.data()))
			return("NA")

			filt.data <- gVars$filt.data()
			nc.data <- gVars$nc.data
			countProbes <- nrow(nc.data)-nrow(filt.data)
			percProbes <- round((countProbes/nrow(nc.data))*100, digits=2)

			percProbesText <- paste0("Probes Removed: ", countProbes, " (", percProbes, ")")
			return(percProbesText)
		})

		output$numProbesText <- renderText({
			if(is.null(gVars$filt.data()))
			return("NA")
			
			filt.data <- gVars$filt.data()
			countProbes <- nrow(filt.data)

			numProbesText <- paste0("Probes Remaining: ", countProbes)
			return(numProbesText)
		})

                observeEvent(input$norm_submit, {
			shiny::validate(
				need(!is.null(gVars$rgList), "CEL files not loaded!"),
				need(!is.null(gVars$filt.data()), "Waiting for the filteration step!")
			)

                        shiny::updateTabsetPanel(session, "display", selected="normTab")

                        progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())
                        progress$set(message="Normalization:", value=0)

			filt.data <- gVars$filt.data()

                        updateProgress(detail="Processing...", value=1/2)
			norm.data <- run.norm(filt.data=filt.data, plot=FALSE)
                        updateProgress(detail="Completed!", value=2/2)

                        gVars$norm.data <- norm.data
                })

                gVars$mart <- reactive({
                        useMart("ENSEMBL_MART_ENSEMBL", dataset=input$org, host="www.ensembl.org")
                })

                gVars$annChoices <- reactive({
                        if(is.null(gVars$mart()))
                        return(NULL)

                        att <- listAttributes(gVars$mart())
                        choicesVec <- att[grep("agilent", att$name),1]
                        choicesNames <- att[grep("agilent", att$name),2]
                        names(choicesVec) <- choicesNames
                        return(choicesVec)
                })

		observeEvent(input$batch_submit, {
			shiny::validate(
				need(!is.null(gVars$celDir()), "No directory selected!")
			)
			shiny::validate(
				need(!is.null(gVars$norm.data), "Normalization is not performed!")
			)

                        shiny::updateTabsetPanel(session, "display", selected="batchTab")

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())
                        progress$set(message="Batch Correction:", value=0)

			varI <- input$varI
                        if(is.null(input$coVar)){
			    coVar <- NULL
                        }else{
			    coVar <- as.list(input$coVar)
                        }
			batch <- as.list(input$batch)
                        #norm.data <- as.matrix(gVars$norm.data)
                        norm.data <- gVars$norm.data
			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
                        
			batchCorVar <- list(var.int=varI, covariates=coVar, batches=batch)
			print(str(batchCorVar))

                        updateProgress(detail="Removing Known Batches...", value=1/4)
			comb.data <- remove.batch.effects(norm.data, phTable, 10, batchCorVar, method="Combat", plot=F, verbose=T)

                        comb.data.cont <- comb.data
                        comb.sva.data <- NULL
                        if(input$rmSva){
                                updateProgress(detail="Removing Unknown Batches...", value=2/4)
                                batches.sva <- get.sva.batch.effects(comb.data=comb.data, pd=phTable, vars=batchCorVar, cmd.ui=F)
                                assoc.cutoff <- 0.05
                                sv.filt.logic <- apply(batches.sva$pd, 2, function(x) x[batchCorVar$var.int]<assoc.cutoff)
                                sv.filt.names <- names(sv.filt.logic[sv.filt.logic==F])
                                if(length(sv.filt.names>0)){
                                        #Check for confounded sva variables
                                        sva.assoc.mat <- assoc.var.int(batches.sva$sv[sv.filt.names], batches.sva$sv[sv.filt.names])
                                        colnames(sva.assoc.mat) <- rownames(sva.assoc.mat)
                                        sva.assoc.mat[lower.tri(sva.assoc.mat, diag=T)] <- NA
                                        sva.assoc.DF <- as.data.frame(as.table(sva.assoc.mat), stringsAsFactors=F)
                                        sva.assoc.DF <- sva.assoc.DF[-which(is.na(sva.assoc.DF$Freq)),]
                                        rowSel <- which(sva.assoc.DF$Freq<0.01)
                                        if(length(rowSel)>0){
                                                sva.assoc.DF.conf <- sva.assoc.DF[rowSel,]
                                                sva.rm.names <- unique(sva.assoc.DF.conf[,2])
                                                sv.filt.names <- sv.filt.names[-which(sv.filt.names %in% sva.rm.names)]
                                        }
                                        #Add selected sva sv and perform combat
                                        sva.sv.filt <- batches.sva$sv[,sv.filt.names]
                                        batchCorVar.sva <- batchCorVar
                                        batchCorVar.sva$batches <- as.list(names(sva.sv.filt))
                                        sva.sv.filt.factor <- as.data.frame(apply(sva.sv.filt, 2, factor))
                                        phTable.sva <- cbind(phTable, sva.sv.filt.factor)
                                        comb.sva.data <- remove.batch.effects(comb.data, phTable.sva, 10, batchCorVar.sva, method="Combat", plot=F, verbose=T)
                                        comb.data.cont <- comb.sva.data
                                }
                        }

			updateProgress(detail="Aggregrating Probes by ID...", value=3/4)
                        if(input$annType=="file"){
			        annDF <- gVars$annDF()
                                idIDX <- which(colnames(annDF)==input$selID)
                                print("Mapping IDs")
                                print(idIDX)
                                #map <- annDF[,c(1,idIDX)]
                                map <- annDF[,c(1:2)]
                        }else if(input$annType=="mart"){
                                mart <- gVars$mart()
                                map <- getBM(mart, attributes=c(input$probeID, input$mapID), filters=c(input$probeID), values=rownames(comb.data.cont))
                                
                                if(nrow(map)==0){
                                    shinyjs::alert("Ensembl query returned 0 mapping!")
                                    return(NULL)
                                }
                                
                                if(length(which(complete.cases(map)==FALSE))>0)
                                map <- map[complete.cases(map),]
                        }
                        colnames(map)[2] <- "TargetID"
			agg.data <- aggreg.probes.2(comb.data.cont, map)
			agg.data <- as.matrix(agg.data)
                        updateProgress(detail="Completed!", value=4/4)

			gVars$batchCorVar <- batchCorVar
			gVars$comb.data <- comb.data
                        gVars$comb.sva.data <- comb.sva.data
			gVars$agg.data <- agg.data
                        gVars$map <- map
		})

		observeEvent(input$add_comp_submit, {
			if(is.null(gVars$condChoices())){
				gVars$comps <- list()
				return(NULL)
			}

			treat <- input$treatment
			control <- input$control
			comps <- paste0(treat, "-", control)
			gVars$comps[[comps]] <- 1
		})

		gVars$compChoices <- reactive({
			if(length(gVars$comps)==0)
			return("No Comparisons Specified!")

			return(names(gVars$comps))
		})

		observeEvent(input$de_submit, {
			shiny::validate(
				need(!is.null(gVars$agg.data), "Batch correction is not performed!")
			)

                        shiny::updateTabsetPanel(session, "display", selected="diffTab")

			progress <- shiny::Progress$new()
			updateProgress <- function(value=NULL, detail=NULL){
				if (is.null(value)) {
					value <- progress$getValue()
					value <- value + (progress$getMax() - value) / 5
				}
				progress$set(value = value, detail = detail)
			}
                        on.exit(progress$close())
                        progress$set(message="Differential Analysis:", value=0)

			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
			batchCorVar <- gVars$batchCorVar
			agg.data <- gVars$agg.data

			updateProgress(detail="Building Model...", value=1/3)
                        print(str(batchCorVar))
			des <- build.model.matrix(phTable, -1, batchCorVar$var.int, batchCorVar$covariates, verbose=T) 
                        annDF <- NULL
                        if(input$annType=="file"){
                                annDF <- gVars$annDF()
			        
                        }else if(input$annType=="mart"){
                                annDF <- gVars$map
                        }
                        annDF <- annDF[-which(duplicated(annDF[,2]) ==T),]
                        rownames(annDF) <- annDF[,2]

			updateProgress(detail="Evaluating Expression...", value=2/3)
			deg.list <- diff.gene.expr(agg.data, des, contrasts=input$comps, pvalue=1, fcvalue=0, p.adjust.method="none", annot=annDF, plot=F, verbose=T)
                        updateProgress(detail="Completed!", value=3/3)

			gVars$deg.list <- deg.list
		})

		output$preBoxPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			filt.data <- gVars$filt.data()
			boxplot(log2(filt.data), las=2, cex=0.7, main="Before Normalization")
		})

		output$postBoxPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			norm.data <- gVars$norm.data
			boxplot(norm.data, las=2, cex=0.7, main="After Normalization")
		})

		output$preDensityPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			rgList <- gVars$rgList
			plotDensities(rgList, main="Before Normalization")
		})

		output$postDensityPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			phTable <- gVars$phTable
			dyeColID <- gVars$dyeColID
			norm.data <- gVars$norm.data

			gSampleCount <- table(phTable[,dyeColID])[[1]]
			rSampleCount <- table(phTable[,dyeColID])[[2]]
			ttlSampleCount <- gSampleCount+rSampleCount

			rgList.norm <- new("RGList", list(G=norm.data[,c(1:gSampleCount)], R=norm.data[,c((gSampleCount+1):ttlSampleCount)]))
			plotDensities(rgList.norm, main="After Normalization")
		})

		output$preMDPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			rgList <- gVars$rgList
			plotMD(rgList, main="Before Normalization")
		})

		output$postMDPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			phTable <- gVars$phTable
			dyeColID <- gVars$dyeColID
			norm.data <- gVars$norm.data

			gSampleCount <- table(phTable[,dyeColID])[[1]]
			rSampleCount <- table(phTable[,dyeColID])[[2]]
			ttlSampleCount <- gSampleCount+rSampleCount
			rgList.norm <- new("RGList", list(G=norm.data[,c(1:gSampleCount)], R=norm.data[,c((gSampleCount+1):ttlSampleCount)]))
			plotMD(rgList.norm, main="After Normalization")
		})

                gVars$confPlotHt <- eventReactive(input$confHeight,{
                        percWd <- input$confWidth
                        percHt <- input$confHeight
                        nWidth <- session$clientData$output_confPlot_width
                        nWidthFull <- (nWidth/percWd)*100
                        nHeight <- (percHt*nWidthFull)/100
                        return(nHeight)
                })

                gVars$confPlotWd <- eventReactive(input$confWidth,{
                        percWd <- input$confWidth
                        nWidth <- session$clientData$output_confPlot_width
                        nWidthFull <- (nWidth/percWd)*100
                        nWidth <- (percWd*nWidthFull)/100
                        return(nWidth)
                })

		output$confPlot <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			print(paste0("Plot Width: ", session$clientData$output_confPlot_width))
			print(paste0("Plot Width: ", session$clientData[["output_confPlot_width"]]))
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			# generate the counfounding plot
			#print(test)
			confounding(ph[,test], margins = c(10,10))
		}, 
                height = function(){
                        percWd <- input$confWidth
                        percHt <- input$confHeight
                        if(is.null(gVars$confDim)){
                                nWidth <- session$clientData$output_confPlot_width
                                nWidthFull <- (nWidth/percWd)*100
                                gVars$confDim <- nWidthFull
                        }else{
                                nWidthFull <- gVars$confDim
                        }
                        nHeight <- (percHt*nWidthFull)/100
                        return(nHeight)
		},
                width = function(){
                        if(is.null(gVars$confDim)){
                                nWidth <- session$clientData$output_confPlot_width
                        }else{
                                percWd <- input$confWidth
                                nWidthFull <- gVars$confDim
                                nWidth <- (percWd*nWidthFull)/100
                        }
                        return(nWidth)
                })

		output$princePlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			print("Accessing PH from extracted data...")
			print(paste0("Plot Width: ", session$clientData$output_princePlot_width))
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$norm.data
			npc = 10
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			pr <- prince(data, ph[,test], top=npc)
			# generate the prince plot
			prince.plot(prince=pr, margins = c(15,15))
		},
                height = function() {session$clientData$output_princePlot_width}#,
                )

		output$hcPlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$norm.data
			print(head(data))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			# hc plot
			hca.plot(data, ph[,test], method = "correlation")
			#list(conf = res, pr = pr)
		})

		output$preCorMDS <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$norm.data), "Waiting for normalization...")
                        )

			norm.data <- gVars$norm.data
			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
			limma:::plotMDS(norm.data, top=500, labels=phTable[,input$varI], col=as.numeric(phTable[,input$varI]), gene.selection="common", main = "Before removing any batch.")
                }, 
                height = function(){
                        percWd <- input$mdsWidth
                        percHt <- input$mdsHeight
                        if(is.null(gVars$mdsDim)){
                                nWidth <- session$clientData$output_confPlot_width
                                nWidthFull <- (nWidth/percWd)*100
                                gVars$mdsDim <- nWidthFull
                        }else{
                                nWidthFull <- gVars$confDim
                        }
                        nHeight <- (percHt*nWidthFull)/100
                        return(nHeight)
		},
                width = function(){
                        if(is.null(gVars$mdsDim)){
                                nWidth <- session$clientData$output_confPlot_width
                        }else{
                                percWd <- input$mdsWidth
                                nWidthFull <- gVars$mdsDim
                                nWidth <- (percWd*nWidthFull)/100
                        }
                        return(nWidth)
		})

		output$postCorMDS <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$comb.data), "Waiting for batch correction...")
                        )

			comb.data <- gVars$comb.data
			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
			limma:::plotMDS(comb.data, top=500, labels=phTable[,input$varI], col=as.numeric(phTable[,input$varI]), gene.selection="common", main = "After removing batches.")
		})

                output$postSvaCorMDS <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$comb.sva.data), "No Unknown Batches Removed...")
                        )

			comb.data <- gVars$comb.sva.data
			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
			limma:::plotMDS(comb.data, top=500, labels=phTable[,input$varI], col=as.numeric(phTable[,input$varI]), gene.selection="common", main = "After removing unknown batches.")
		})

		output$postAggMDS <- renderPlot({
                        shiny::validate(
                                need(!is.null(gVars$agg.data), "Waiting for batch correction...")
                        )

			agg.data <- gVars$agg.data
			phTable <- gVars$phTable
			phTable <- as.data.frame(apply(phTable, 2, factor))
			limma:::plotMDS(agg.data, top=500, labels=phTable[,input$varI], col=as.numeric(phTable[,input$varI]), gene.selection="common", main = "After removing batches.")
		})

		output$postHcPlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$comb.data), "Waiting for batch correction...")
                        )
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$comb.data
			print(head(data))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			# hc plot
			hca.plot(data, ph[,test], method = "correlation")
			#list(conf = res, pr = pr)
		})

                output$postSvaHcPlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$comb.sva.data), "No unknown batches removed...")
                        )
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$comb.sva.data
			print(head(data))
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			# hc plot
			hca.plot(data, ph[,test], method = "correlation")
			#list(conf = res, pr = pr)
		})

		output$postPrincePlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$comb.data), "Waiting for normalization...")
                        )
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$comb.data
			npc = 10
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			pr <- prince(data, ph[,test], top=npc)
			# generate the prince plot
			prince.plot(prince=pr, margins = c(15,15))
		})

                output$postSvaPrincePlot <- renderPlot({			
                        shiny::validate(
                                need(!is.null(gVars$comb.sva.data), "No unknown batches removed...")
                        )
			cat("Accessing PH from extracted data...")
			cat("\n")
			ph <- gVars$phTable
			ph <- as.data.frame(apply(ph, 2, factor))
			data <- gVars$comb.sva.data
			npc = 10
			test <- sapply(colnames(ph), function(b) length(table(ph[,b])) > 1 & length(table(ph[,b])) != length(ph[,b]))
			pr <- prince(data, ph[,test], top=npc)
			# generate the prince plot
			prince.plot(prince=pr, margins = c(15,15))
		})

		output$deTable <- DT::renderDataTable({
			shiny::validate(
                                need(!is.null(gVars$deg.list), "Waiting for annotation enrichment results..."),
                                need(length(gVars$deg.list)>0, "Waiting for annotation enrichment results...")
                        )

                        deg.list <- gVars$deg.list
			comp <- input$compDE
			print(comp)
			print(names(deg.list))
			print(head(deg.list[[comp]]))
			print("Check Complete###############")

			DT::datatable(deg.list[[comp]], filter=list(position='top', clear=FALSE), options=list(search=list(regex=TRUE, caseInsensitive=FALSE), scrollX=TRUE))
		},server=TRUE)

		output$volcanoPlot <- renderPlot({
			deg.list <- gVars$deg.list
			comp <- input$compDE
			deg <- deg.list[[comp]]
			deg <- data.frame(x=as.numeric(deg$logFC), y=-log10(as.numeric(deg$P.Value)), ID=rownames(deg))
			p <- ggplot(deg, aes(x, y, label= ID)) + geom_point() +
			geom_vline(xintercept = input$lfcThr, color = "blue") + 
			geom_vline(xintercept = -input$lfcThr, color = "blue") + 
			geom_hline(yintercept = input$pvThr, color = "red") +  
			labs(x="log2(Fold-change)", y="-log10(P.Value)") + 
			scale_x_continuous("log2(Fold-change)", limits = range(-input$xAxis,input$xAxis)) +
			scale_y_continuous("-log10(P.Value)", limits = range(0,input$yAxis)) + theme_bw()

			if(input$gName){
			    degSel <- deg[deg$y>input$pvThr & abs(deg$x)>input$lfcThr,]
			    p <- p + annotate("text", x=degSel$x, y=degSel$y, 
			    label=degSel$ID, size=degSel$y, 
			    vjust=-0.1, hjust=-0.1)
			}
			print(p)
		})

		output$phCounts <- renderPrint({
			phRows <- input$ph_rows_selected
			if (length(phRows)) {
				cat(phRows, sep = ', ')
			}
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

		output$selFileNameCol <- renderUI({
			selectInput("fileNameCol", "Select Filename Column", choices=gVars$phColChoices())
		})

		output$selDyeCol <- renderUI({
			selectInput("dyeCol", "Select Dye Column", choices=gVars$phColChoices())
		})

		output$selSep <- renderUI({
			selectInput("sepS", "Select Column Seperator", choices=gVars$sepChoices, selected=gVars$sepChoices[1])
		})

		output$selVarI <- renderUI({
			selectInput("varI", "Select Variable of Interest", choices=gVars$pcChoices())
		})
		
		output$selCoVar <- renderUI({
			selectInput("coVar", label="Select Co-Variates", choices=gVars$pcChoices(), multiple=TRUE)
		})
		
		output$selBatch <- renderUI({
			selectInput("batch", "Select Batches", choices=gVars$pcChoices(), multiple=TRUE)
		})
		
                output$selProbeID <- renderUI({
			selectInput("probeID", "Select Probe ID", choices=gVars$annChoices(), multiple=FALSE)
		})
		
		output$slideFiltDist <- renderUI({
			sliderInput("filtDist", "Qunatile Based Cutoff", min=0.10, max=1.0, value=0.75, step=0.05, round=2)
		})
		
		output$slidePerSamples <- renderUI({
			sliderInput("perSamples", "Select Percentage of Samples", min=1, max=100, value=75, step=1)
		})

		output$selNormMethod <- renderUI({
			selectInput("normMethod", "Select Normalization Method", choices=gVars$normChoices, multiple=FALSE)
		})
		
		output$selID <- renderUI({
			selectInput("ID", "Select ID", choices=gVars$idChoices(), multiple=FALSE)
		})
		
		output$selTreatment <- renderUI({
			selectInput("treatment", "Condition 1", choices=gVars$condChoices(), multiple=FALSE)
		})
		
		output$selControl <- renderUI({
			selectInput("control", "Condition 2", choices=gVars$condChoices(), multiple=FALSE)
		})
		
		output$selComps <- renderUI({
			selectInput("comps", "Comparisons", choices=gVars$compChoices(), multiple=TRUE, selected=gVars$compChoices())
		})

		output$selPvAdjMethod <- renderUI({
			selectInput("pvAdjMethod", "P.Value Adjustment Method", choices=gVars$pvAdjChoices, multiple=FALSE, selected="none")
		})

		output$selCompDeTable <- renderUI({
			selectInput("compDE", "Select Comparison", choices=gVars$compChoices(), multiple=FALSE)
		})
		
		shinyjs::onclick("showPhenoOptions", shinyjs::toggle(id="phenoOptions", anim=TRUE))
		shinyjs::onclick("showPhenoOptions", shinyjs::hide(id="rawOptions", anim=TRUE), add=T)
		shinyjs::onclick("showPhenoOptions", shinyjs::hide(id="filtOptions", anim=TRUE), add=T)
		shinyjs::onclick("showPhenoOptions", shinyjs::hide(id="normOptions", anim=TRUE), add=T)
		shinyjs::onclick("showPhenoOptions", shinyjs::hide(id="batchOptions", anim=TRUE), add=T)
		shinyjs::onclick("showPhenoOptions", shinyjs::hide(id="diffOptions", anim=TRUE), add=T)

		shinyjs::onclick("showRawOptions", shinyjs::toggle(id="rawOptions", anim=TRUE))
		shinyjs::onclick("showRawOptions", shinyjs::hide(id="phenoOptions", anim=TRUE), add=T)
		shinyjs::onclick("showRawOptions", shinyjs::hide(id="filtOptions", anim=TRUE), add=T)
		shinyjs::onclick("showRawOptions", shinyjs::hide(id="normOptions", anim=TRUE), add=T)
		shinyjs::onclick("showRawOptions", shinyjs::hide(id="batchOptions", anim=TRUE), add=T)
		shinyjs::onclick("showRawOptions", shinyjs::hide(id="diffOptions", anim=TRUE), add=T)

		shinyjs::onclick("showFiltOptions", shinyjs::toggle(id="filtOptions", anim=TRUE))
		shinyjs::onclick("showFiltOptions", shinyjs::hide(id="phenoOptions", anim=TRUE), add=T)
		shinyjs::onclick("showFiltOptions", shinyjs::hide(id="rawOptions", anim=TRUE), add=T)
		shinyjs::onclick("showFiltOptions", shinyjs::hide(id="normOptions", anim=TRUE), add=T)
		shinyjs::onclick("showFiltOptions", shinyjs::hide(id="batchOptions", anim=TRUE), add=T)
		shinyjs::onclick("showFiltOptions", shinyjs::hide(id="diffOptions", anim=TRUE), add=T)

		shinyjs::onclick("showNormOptions", shinyjs::toggle(id="normOptions", anim=TRUE))
		shinyjs::onclick("showNormOptions", shinyjs::hide(id="phenoOptions", anim=TRUE), add=T)
		shinyjs::onclick("showNormOptions", shinyjs::hide(id="rawOptions", anim=TRUE), add=T)
		shinyjs::onclick("showNormOptions", shinyjs::hide(id="filtOptions", anim=TRUE), add=T)
		shinyjs::onclick("showNormOptions", shinyjs::hide(id="batchOptions", anim=TRUE), add=T)
		shinyjs::onclick("showNormOptions", shinyjs::hide(id="diffOptions", anim=TRUE), add=T)

		shinyjs::onclick("showBatchOptions", shinyjs::toggle(id="batchOptions", anim=TRUE))
		shinyjs::onclick("showBatchOptions", shinyjs::hide(id="rawOptions", anim=TRUE), add=T)
		shinyjs::onclick("showBatchOptions", shinyjs::hide(id="phenoOptions", anim=TRUE), add=T)
		shinyjs::onclick("showBatchOptions", shinyjs::hide(id="filtOptions", anim=TRUE), add=T)
		shinyjs::onclick("showBatchOptions", shinyjs::hide(id="normOptions", anim=TRUE), add=T)
		shinyjs::onclick("showBatchOptions", shinyjs::hide(id="diffOptions", anim=TRUE), add=T)

		shinyjs::onclick("showDiffOptions", shinyjs::toggle(id="diffOptions", anim=TRUE))
		shinyjs::onclick("showDiffOptions", shinyjs::hide(id="phenoOptions", anim=TRUE), add=T)
		shinyjs::onclick("showDiffOptions", shinyjs::hide(id="rawOptions", anim=TRUE), add=T)
		shinyjs::onclick("showDiffOptions", shinyjs::hide(id="filtOptions", anim=TRUE), add=T)
		shinyjs::onclick("showDiffOptions", shinyjs::hide(id="normOptions", anim=TRUE), add=T)
		shinyjs::onclick("showDiffOptions", shinyjs::hide(id="batchOptions", anim=TRUE), add=T)

		observe({
                        if(input$annType=="mart"){
			        shinyjs::hide("annFile")
				shinyjs::show("annQuery")
                        }else if(input$annType=="file"){
			        shinyjs::hide("annQuery")
				shinyjs::show("annFile")
                        }

			if(!is.null(input$sepS) && input$sepS=="OTHER"){
				shinyjs::enable("sepT")
			}else{
				shinyjs::disable("sepT")
			}

			if(!is.null(input$treatment) && !is.null(input$control)){
				if(input$treatment==input$control){
					shinyjs::disable("add_comp_submit")
				}else{
					shinyjs::enable("add_comp_submit")
				}
			} 
		})
	}
)
