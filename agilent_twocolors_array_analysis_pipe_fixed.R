## suppose to have one file for each slot, 
## for a two-color experiment we have two samples associated with one slot
agilent.twocolors.analysis <- function(pd.file, id.rem.samp = NULL, type.data = "mRNA",
                                       qdist = c(.9), perc = c(50), num.pc = 5, gr = TRUE,
                                       pdf.file = NULL, save.excel = NULL, 
                                       verbose = TRUE, plot = TRUE, ...) {  
  require(limma)
  require(sva)
  require(swamp)
  require(infotheo)
  require(biomaRt)
  require(gplots)
  require(Biobase)
  #require(BACA)
  require(RColorBrewer)
  
  # loading methada, raw data, annotations and pheno data 
  pd <- read.table(pd.file, sep="\t", header=T)
  id_file_names <- grep("file", colnames(pd), ignore.case = TRUE)
  if(length(id_file_names) == 0) stop("The <file.name> field is missing!")
  rgList <- read.maimages(files=unique(pd[,id_file_names]), 
                          source="agilent.median", verbose = verbose)
  # generate nc.data, c.data and pd
  st <- fix.data(pd, id_file_names, rgList, GR = T, id.rem.samp)
  pd <- st$pd
  # opening a pdf file to save the generated plot
  if(!is.null(pdf.file))  pdf(pdf.file)
  if(verbose) cat("Analysis of the sample annotations.. \n")
  monitor.technical.variation(st$nc.data, st$pd)
  # extract from the pd file the variables of interest, the covariates and the technical variables
  if(verbose)  cat("Extracting phenotype data and ..\n")
  vct <- get.info(pd)
  if(verbose)  print(vct)
  if(verbose)  cat("Pre-processing data <",nrow(st$nc.data),",",ncol(st$nc.data),"> \n", sep="")
  norm.data <- pre.proc(nc.data = st$nc.data, c.data = st$c.data, pd[,vct$var.int], qdist = qdist, perc = perc, 
                        verbose = verbose, plot = plot)
  if(verbose) cat("Removing batch effects..<iterative process>\n")
  removing.batches <- c("Combat", "Limma")
  print(data.frame(rb.methods = removing.batches))
  idm <- read.input.number(" - indicate the method to use to remove batch effects:")
  if(idm == 1 | idm == 2)
    comb.data <- remove.batch.effects(norm.data, pd, num.pc, vct, method = removing.batches[idm], plot = plot, verbose = verbose)
  else 
    stop("invalid method")
  if(verbose) cat("Analysis of the sample annotations .after removing the batch effects.. \n")
  monitor.technical.variation(comb.data, st$pd)
  if(verbose) cat("Connecting to BioMart..\n")
  con <- connect.to.biomart.2()
  if(verbose) cat("Annotating..\n")
  biomart.data <- query.biomart(con$mart, con$organism, comb.data)
  if(verbose) cat("Combining probe expression values..\n")
  fdata <- aggreg.probes(comb.data, biomart.data$map, var.mds = pd[,vct$var.int], plot = plot)
  plot.heatmap(fdata, pd[,vct$var.int], 500)
  # return(list(data = fdata, biomart.data = biomart.data, vct = vct, pd = pd)) ### -------------------------> ATTENTION
  # return(list(norm.data = norm.data, combat = comb.data, data = fdata, biomart.data = biomart.data, vct = vct, pd = pd))
  if(verbose) cat("Differential gene expression analysis of ", type.data, "..\n", sep="")
  des <- build.model.matrix(pd, -1, vct$var.int, vct$covariates, verbose = verbose) 
  deg <- diff.gene.expr(fdata, des, vct$contrasts, biomart.data$annot, 
                        save.file = save.excel, plot = plot, verbose = verbose)
  if(!is.null(pdf.file)) dev.off()
  list(norm = norm.data, 
       combat = comb.data, 
       data = fdata, 
       BM = biomart.data, 
       deg = deg,
       mod = vct)
}

## PRE-PROCESSING

get.info <- function(pd, verbose = TRUE, ...) {
  vars <- colnames(pd)
  print(data.frame(Variables = vars))
  n1 <- NULL
  while(is.null(n1)) {  
    n1 <- read.input.number(" - indicate the variable of interest:")
    if(!is.numeric(n1)) n1 <- NULL
    if(n1 > length(vars)) n1 <- NULL
  }
  var.int <- vars[n1]
  vars <- vars[-n1]
  # n1 <- read.input.number(" - how many contrasts:")
  levs <- names(table(pd[,var.int]))
  print(data.frame(Levels = levs))
  contrasts <- c("")
  flag <- TRUE
  while(flag) {
    cont <- read.input.word(message = " - indicate a contrast:")
    cont <- as.numeric(unlist(strsplit(cont, ",")))
    contrasts <- c(contrasts, paste(levs[cont[1]],levs[cont[2]],sep="-")) 
    flag <- continue.read.inputs()
  }
  contrasts <- contrasts[-1]
  n2 <- NULL
  # indicating the covariates
  covariates <- NULL
  while(!is.numeric(n2)) {
    print(data.frame(Variables = vars))
    n2 <- read.input.number(" - add a covariate <enter '0' to go to the next step>:")
    if(n2 == 0) {
      n2 <- NULL
      break
    }
    if(is.null(covariates))  covariates <- n2
    else  covariates <- c(covariates, n2)
    n2 <- NULL
  }
  if(!is.null(covariates)) {
    covs <- vars[covariates]
    vars <- vars[-covariates]
  }
  else  covs <- NULL
  # indicating the batch variables 
  batches <- NULL
  while(!is.numeric(n2)) {
    print(data.frame(Variables = vars))
    n2 <- read.input.number(" - add a batch <enter '0' to go to the next step>:")
    if(n2 == 0) break
    if(is.null(batches)) batches <- n2
    else batches <- c(batches, n2)
    n2 <- NULL
    print(batches)
  }
  if(!is.null(batches)) bats <- vars[batches]
  else bats <- NULL
  list(var.int = var.int,
       covariates = covs,
       batches = bats,
       contrasts = contrasts)
}
read.contrasts <- function(group) {
  levs <- names(table(group))
  print(data.frame(Levels = levs))
  contrasts <- c("")
  flag <- TRUE
  while(flag) {
    cont <- read.input.word(message = " - indicate a contrast:")
    cont <- as.numeric(unlist(strsplit(cont, ",")))
    contrasts <- c(contrasts, paste(levs[cont[1]],levs[cont[2]],sep="-")) 
    flag <- continue.read.inputs()
  }
  contrasts <- contrasts[-1]
  contrasts
}
run.norm <- function(filt.data, rgList, method, method2, plot=FALSE){
  norm.data <- NULL

  if(plot)  {
   boxplot(log2(filt.data), las=2, cex=0.7, main="Before normalization.")
  }
  norm_switch <- function(method){
   switch(method,
    "quantile" = normalizeQuantiles(log2(filt.data)),
    "vsn" = normalizeVSN(log2(filt.data)),
    "cl" = normalizeCyclicLoess(log2(filt.data), method="fast"),
    "BA" = normalizeBetweenArrays(log2(filt.data), method=method2)
   )
  }
  #norm.data <- normalizeQuantiles(log2(filt.data)) 
  norm.data <- norm_switch(method) 
  if(plot)  {
   boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
  }

  return(norm.data)
} 
pre.proc <- function(pd, rgList, var.mds, qdist = c(.75, .9), perc = c(75, 50), verbose = TRUE, plot = TRUE, ...) {
  norm.data <- NULL
  data <- cbind(rgList$G, rgList$R)
  rownames(data) <- rgList$genes$ProbeName
  # defining the colnames for data
  print(head(colnames(data)))
  print(head(rownames(pd)))
  colnames(data) <- rownames(pd)
  nc.data <- data[which(rgList$genes$ControlType==0),]
  c.data  <- data[which(rgList$genes$ControlType==-1),]
  # filter with different combinations of qdist and perc
  if(length(qdist) != 1 & length(perc) != 1) {
    grid <- expand.grid(qdist, perc)
    filt.data <- apply(grid, 1, FUN = function(x) {
                          if(verbose) cat("*** (qidst:",x[1],", perc.samples:",x[2], ")\n", sep="")
                          #fd <- 
                          filt.by.neg.probes(nc.data, c.data, qdist = x[1], perc = x[2], verbose)
                          #filt.by.var.dup.probes(fd, qdist = x[1], perc = x[2], verbose)
                          })
    # print out some information
    grid.info <- apply(grid, 1, function(x) paste(x[1], x[2], sep="-"))
    perc.probes <- unlist(lapply(filt.data, function(f) {
      round(100 - (((nrow(nc.data) - nrow(f))/nrow(nc.data))*100), digits = 2)}))
    print(data.frame(grid.info, perc.probes, num.probes = unlist(lapply(filt.data, nrow))))
    #n1 <- read.input.number(message = paste("Insert a number between 1 and", length(grid.info)), more.info=NULL)
    n1 <- 1 
    # normalize the data
    if(verbose) cat(" - quantile normalization method..\n")
    if(plot)  {
      boxplot(log2(filt.data[[n1]]), las=2, cex=0.7, main="Before normalization.")
    }
    norm.data <- normalizeQuantiles(log2(filt.data[[n1]]))
    if(plot)  {
      boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
    }
  }
  else {
    if(verbose) cat("*** (qdist:",qdist,", perc.samples:",perc, ")\n", sep="")
    filt.data <- filt.by.neg.probes(nc.data, c.data, qdist = qdist, perc = perc, verbose)
    #filt.data <- filt.by.var.dup.probes(filt.data, qdist = qdist, perc = perc, verbose)
    cat("perc.probes", round(100 - (((nrow(nc.data) - nrow(filt.data))/nrow(nc.data))*100), digits = 2), "\n")
    cat("num.probes", nrow(filt.data), "\n")
    if(plot)  {
     boxplot(log2(filt.data), las=2, cex=0.7, main="Before normalization.")
    }
    norm.data <- normalizeQuantiles(log2(filt.data)) 
    if(plot)  {
     boxplot(norm.data, las=2, cex=0.7, main="After normalization.")
    }
  }
  return(norm.data)
} 
filt.by.neg.probes <- function(data, cdata, qdist = .9, perc = 50, verbose = TRUE, ...) {
  if(verbose) cat(" - filtering by control-negative probes..\n")
  if(verbose) print(paste0("Qdist - ", qdist, ", Perc - ", perc))
  print("nrow(data)")
  print(nrow(data))
  # compile a threshold value from the negative control probes for each sample
  thres.vec <- apply(cdata, 2, function(x) quantile(x, qdist, na.rm=TRUE))
  # filter out the probes whose...
  score.neg <- matrix(0, nrow=nrow(data), ncol=ncol(data))
  for(i in 1:dim(data)[2])  score.neg[which(data[,i] > thres.vec[i]),i] <- 1
  neg.sumrow <- apply(score.neg, 1, sum)
  if(verbose) print(paste0("Rounded - ", round(ncol(data)*perc/100)))
  filt.data <- data[neg.sumrow >= round(ncol(data)*perc/100),]
  print("nrow(filt.data)")
  print(nrow(filt.data))
  if(verbose) cat("     * number of probes to remove:", (nrow(data) - nrow(filt.data)), "\n")
  return(filt.data)
}
filt.by.var.dup.probes <- function(data, qdist = .9, perc = 50, verbose = TRUE, plot = TRUE, ...) {
  if(verbose) cat(" - filtering by variance duplicated probes..\n")
  # build the list of indexes for each duplicated non control probe
  un.probes <- unique(rownames(data))
  list.id.dup.probes <- lapply(un.probes, function(p) {which(rownames(data) == p)})
  un.probes <- un.probes[-which(sapply(list.id.dup.probes, length) < 2)]
  list.id.dup.probes <- list.id.dup.probes[-which(sapply(list.id.dup.probes, length) < 2)]
  score.neg <- matrix(0, nrow=length(list.id.dup.probes), ncol=ncol(data))
  for(j in 1:ncol(data)) {
    dist.var <- unlist(lapply(list.id.dup.probes, function(i) var(data[i,j])))
    thres.vec <- quantile(dist.var, qdist, na.rm=TRUE)
    score.neg[which(dist.var > thres.vec), j] <- 1
  }
  neg.sumrow <- apply(score.neg, 1, sum)
  probs.to.rem <- un.probes[neg.sumrow >= round(ncol(data)*perc/100)]
  probs.to.rem <- which((rownames(data) %in% probs.to.rem))
  if(verbose) cat("     * number of probes to remove:", length(probs.to.rem), "\n")
  return(data[-probs.to.rem,])
}
## VISUALIZE/IDENTIFY & REMOVE BATCH EFFECTS
monitor.technical.variation <- function(data, pd, npc = 10, verbose = TRUE, ...) {
  test <- sapply(colnames(pd), function(b) length(table(pd[,b])) > 1 & length(table(pd[,b])) != length(pd[,b]))
  # generate the counfounding plot
  print(test)
  res <- confounding(pd[,test], margins = c(10,10))
  pr <- prince(data, pd[,test], top=npc)
  print(pr)
  # generate the prince plot
  prince.plot(prince=pr, margins = c(15,15))
  # hc plot
  hca.plot(data, pd[,test], method = "correlation")
  list(conf = res, pr = pr)
}
pc.anaylsis <- function(data, pd, npc = 10, verbose = TRUE, ...) {
  print("PC-analysis")
  # run the principal-component analysis
  var.names <- colnames(pd)
  # remove variables with one level
  test <- sapply(var.names, function(v) { 
                 if(is.factor(pd[,v]))
                   if(nlevels(pd[,v]) == 1) TRUE
                   else FALSE
                 else FALSE}) 
  ids <- which(test == TRUE)
  if(length(ids) > 0) var.names <- var.names[-ids]
  # run princomp
  pc <- princomp(data, cor=T)
  eig <- with(pc, unclass(loadings))
  if(ncol(eig) > npc) eig <- eig[,1:npc]
  perc.importance <- round(pc$sdev^2/sum(pc$sdev^2)*100, 2)
  # compile assoc matrix
  assoc <- matrix(NA, nrow=length(var.names), ncol=ncol(eig))
  rownames(assoc) <- var.names
  colnames(assoc) <- colnames(eig)
  for(i in 1:npc) {
    assoc[,i] <- sapply(var.names, function(a) {
      if(is.factor(pd[,a])) 
        anova(lm(eig[,i]~as.factor(pd[,a])))[1,5]
      else if(is.numeric(pd[,a]))
        anova(lm(eig[,i]~pd[,a]))[1,5]
    })
  }
  colnames(assoc) <- paste(rep("PC", npc), 1:npc, " (", perc.importance[1:npc], "%) ", sep="")
  return(assoc)
}
pc.anaylsis.2 <- function (g, o, npc = 10, center = T){
  test <- sapply(colnames(o), function(b) length(table(o[,b])) > 1 & length(table(o[,b])) != length(o[,b]))
  if(nrow(o) == 1) {
    o <- as.data.frame(o)
    colnames(o) <- "sva.1"
  }
  print(test)
  if (is.matrix(g) != T) {
    stop("g is not a matrix")
  }
  if (is.data.frame(o) != T) {
    stop("o is not a data.frame")
  }
  classes <- unlist(lapply(unclass(o), class))
  if (any(classes == "character")) {
    stop("o contains characters")
  }
  nrlevels <- unlist(lapply(unclass(o), function(x) length(levels(x))))
  if (any(nrlevels == 1)) {
    stop("o contains factors with only one level")
  }
  if (npc > ncol(g)) {
    stop("npc is larger than ncol(g)")
  }
  if (npc > nrow(g)) {
    stop("npc is larger than nrow(g)")
  }
  if (!identical(rownames(o), colnames(g))) {
    warning("Colnames of g are not the same as rownames of o")
  }
  if (center == T) {
    pr <- prcomp(t(g))
  }
  if (center == F) {
    pr <- prcomp(t(g), center = F)
  }
  linp <- matrix(ncol = npc, nrow = ncol(o))
  rownames(linp) <- colnames(o)
  rsquared <- matrix(ncol = npc, nrow = ncol(o))
  rownames(rsquared) <- colnames(o)
  for (i in 1:ncol(o)) {
    for (j in 1:npc) {
      fit <- lm(pr$x[, j] ~ o[, i])
      s <- summary(fit)
      linp[i, j] <- pf(s$fstatistic[1], s$fstatistic[2], 
                       s$fstatistic[3], lower.tail = FALSE)
      rsquared[i, j] <- s$r.squared[1]
    }
  }
  prop <- (pr$sdev[1:npc]^2/sum(pr$sdev^2)) * 100
  assoc <- linp
  colnames(assoc) <- paste(rep("PC", npc), 1:npc, " (", round(prop[1:npc],2), "%) ", sep="")
  return(assoc)
}
assoc.var.int <- function(X, pd, verbose = TRUE, ...) {
  # run the principal-component analysis
  var.names <- colnames(pd) 
  # remove variables with one level
  test <- sapply(var.names, function(v) { 
    if(is.factor(pd[,v]))
      if(nlevels(pd[,v]) == 1) TRUE
    else FALSE
    else FALSE}) 
  ids <- which(test == TRUE)
  if(length(ids) > 0) var.names <- var.names[-ids]
  # compile assoc matrix
  assoc <- matrix(NA, nrow=length(var.names), ncol=ncol(X))
  rownames(assoc) <- var.names
  colnames(assoc) <- colnames(colnames(X))
  print("assoc:")
  print(assoc)
  print("str(assoc):")
  print(str(assoc))
  for(i in 1:ncol(assoc)) {
    cat("assoc col : ", i, "\n")
    assoc[,i] <- sapply(var.names, function(a) {
      cat("pd var name : ", a, "\n")
      if(is.factor(pd[,a])) {
        cat("Is a factor\n")
        if(length(levels(pd[,a])) > 0)
          anova(lm(X[,i]~as.factor(pd[,a])))[1,5]
      }
      else if(is.numeric(pd[,a]))
        cat("Is numeric\n")
        anova(lm(X[,i]~pd[,a]))[1,5]
    })
  }
  return(assoc)
}
remove.batch.effects <- function(comb.data, pd, num.pc, vars, inter = "", method = "Combat", plot = TRUE, verbose = TRUE, ...) {
  batches <- vars$batches
  covs <- vars$covariates
  if(plot){
      cat(" - Running Principal Component Analysis..\n")
      assoc <- pc.anaylsis.2(g = comb.data, o = pd, npc = num.pc)
      print(assoc)
      limma:::plotMDS(comb.data, top=500, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), gene.selection="common", main = "Before removing any batch.")
  }
  for(batch in batches){
      id.var <- grep(batch, batches)[1]
      print(id.var)
      pre.batches <- batches
      batches <- batches[-id.var]
      string.formula <- paste("~", inter, "pd$", vars$var.int, sep="")
      if(length(batches) > 0)
        string.formula <- paste(string.formula, "+pd$", paste(batches, collapse="+pd$"), sep="")
      if(length(covs) > 0)
        string.formula <- paste(string.formula, "+pd$", paste(covs, collapse="+pd$"), sep="")
      form <- formula(string.formula)
      print(form)
      prev.comb.data <- comb.data
      if(method == "Combat") 
        res <- try(comb.data <- ComBat(dat=comb.data, batch=pd[,batch], 
                                     mod=model.matrix(form), 
                                     par.prior=T,  prior.plots=F))
      if(method == "Limma")
        res <- try(comb.data <- removeBatchEffect(x = comb.data,
                                                  batch = pd[,batch],
                                                  design = model.matrix(form)))
      
      if(inherits(res, "try-error")) {
        ##error handling: restore the previous state
        #comb.data <- prev.comb.data 
        #batches <- pre.batches
        #print(c(vars$var.int,batches))
        #next

        #error handling: return the custom error statement, with variable names
        error.str <- paste0("Variables are confounded!\n\nPlease remove one or more variable, so the design is not confunded!\n\nError encountered while correcting for BATCH:", batch, "\n\nPlease check ComBat resource for more information!")
        return(error.str)
      }
      if(plot){
          cat(" - Running Principal Component Analysis..\n")
          assoc <- pc.anaylsis.2(g=comb.data, o=pd, npc=num.pc)
          print(assoc)
          limma:::plotMDS(comb.data, top=500, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]),gene.selection="common", main = paste("After removing", batch))
      }
  }
  if(plot)  {
    boxplot(comb.data, las=2, cex=0.7, main="After removing batch effects.")
    #plot.heatmap(comb.data, pd[,vars$var.int], 1000)
  }
  return(comb.data)
}
remove.batch.effects.2 <- function(comb.data, pd, num.pc, vars, inter = "", method = "Combat", plot = TRUE, verbose = TRUE, ...) {
  batches <- vars$batches
  covs <- vars$covariates
  cat(" - Running Principal Component Analysis..\n")
  assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
  print(assoc)
  if(plot) limma:::plotMDS(comb.data, top=1000, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), 
                           gene.selection="common", main = "Before removing any batch.")
  #limma:::plotMDS(comb.data, top=500, labels=pd[,"TimeExposition"], col=as.numeric(pd[,"TimeExposition"]), 
  #                gene.selection="common", main = "Before removing any batch.")
  flag <- continue.read.inputs()
  while(flag) {
    batch  <- read.input.word(message = " - Please enter a variable <name> to use as batch..\n", 
                              more.info = paste("batches:", paste(batches, collapse=",")))
    id.var <- which(batches == batch)
    print(id.var)
    if(length(id.var) > 0) {
      prev.comb.data <- comb.data
      pre.batches <- batches
      print(c(covs, batches)[(length(covs) + id.var)])
      comb.data <- swamp:::combat(comb.data, pd[,c(covs, batches)], length(covs) + id.var, par.prior=T,  prior.plots=F)
      batches <- batches[-id.var]
      cat(" - Running Principal Component Analysis..\n")
      assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
      print(assoc)
      if(plot) limma:::plotMDS(comb.data, top=1000, labels=pd[,vars$var.int], col=as.numeric(pd[,vars$var.int]), 
                               gene.selection="common", main = paste("After removing", batch))
      # roll back if the user wants to
      rb <- continue.read.inputs(message = "Do you want to roll back? (y/n)")
      if(rb == TRUE) {
        comb.data <- prev.comb.data 
        batches <- pre.batches
        cat(" - Previous Principal Component Analysis..\n")
        assoc <- pc.anaylsis(data = comb.data, pd = pd, npc = num.pc)
        print(assoc)
      }
    }
    else {
      cat(" - Invalid variable name.\n")
      next
    }
    flag <- continue.read.inputs()
  }
  if(plot)  boxplot(comb.data, las=2, cex=0.7, main="After removing batch effects.")
  return(comb.data)
}
get.sva.batch.effects <- function(comb.data, pd, vars, npc = 10, verbose = T, cmd.ui = T) {
  print("Searching for unknown batch effects:")
  covs <- vars$covariates
  string.formula <- paste("~pd$", vars$var.int, sep="")
  if(length(covs) > 0)  string.formula <- paste(string.formula, "+pd$", paste(covs, collapse="+pd$"), sep="")
  print(string.formula)
  form <- formula(string.formula)
  print(form)
  X <- sva(dat=comb.data, mod=model.matrix(form), method = "two-step")$sv
  cat("class(X) - ", class(X), "\n")
  cat("dim(X) - ", dim(X), "\n")
  X.c <- X
  cat("class(X.c) - ", class(X.c), "\n")
  cat("dim(X.c) - ", dim(X.c), "\n")
  X <- discretize(as.matrix(X), disc="equalfreq", nbins=NROW(X)^(1/3))
  X <- as.data.frame(X)
  if(verbose) print(X)
  colnames(X) <- paste("sva",c(1:ncol(X)),sep=".")
  #colnames(X) <- paste("svaD",c(1:ncol(X)),sep=".")
  rownames(X) <- colnames(comb.data)
  #colnames(X.c) <- paste("svaC",c(1:ncol(X)),sep=".")
  colnames(X.c) <- colnames(X)
  rownames(X.c) <- rownames(X)
  assoc.pca <- pc.anaylsis.2(g = comb.data, o = X, npc = npc)
  if(verbose) print(assoc.pca)
  assoc.pd <- assoc.var.int(X, pd)
  colnames(assoc.pd) <- colnames(X)
  if(verbose) print(assoc.pd)
  # indicating the batch variables found by uysing SVA
  if(cmd.ui){
    batches <- read.input.word(message = " - indicate batches <n to exit>: ")
    if(batches != "n") {
      ids <- as.numeric(unlist(strsplit(batches, ",")))
      X[,ids]
    }else NULL
  }else list(pd=assoc.pd, sv=X, svc=X.c)
}
## ANNOT-PROCESSING functions
connect.to.biomart <- function(verbose = TRUE, ...) {
  marts <- listMarts(host="www.ensembl.org")
  print(as.data.frame(marts[,1]))
  n1 <- read.input.number("Please enter a number to indicate the BioMart database to use.")
  if(verbose) cat(" - selected BioMart database:", as.character(marts[n1,1]), "\n", sep="")
  mart <- useMart(as.character(marts[n1,1]))
  datasets <- listDatasets(mart)
  organism <- read.input.word(1, "mmusculus or hsapiens")
  print(class(organism))
  ids <- grep(organism, datasets[,1])
  if(length(ids) > 0) cat(" - selected organism:", as.character(datasets[ids,1]), "\n", sep="")
  else stop("dataset not found")
  mart <- useMart(biomart=as.character(marts[n1,1]), dataset=as.character(datasets[ids,1]))
  mart
}
connect.to.biomart.2 <- function(verbose = TRUE, ...) {
  organisms <- c("hsapiens_gene_ensembl","mmusculus_gene_ensembl") 
  print(as.data.frame(organisms, stringsAsFactors = F))
  cat("Please select the BioMart database to use...\n")
  n1 <- scan(n=1)     
  if(is.numeric(n1)) 
    list(mart = useMart("ENSEMBL_MART_ENSEMBL", dataset=organisms[n1], host="www.ensembl.org"), organisms[n1])
  else
    NULL
}
query.biomart <- function(mart, organism, data, verbose = TRUE, ...) {
  probes <- as.character(unique(rownames(data)))
  if(interactive()){
    agilent.annot <- c("efg_agilent_wholegenome_4x44k_v1",  "efg_agilent_wholegenome_4x44k_v2",
                       "efg_agilent_sureprint_g3_ge_8x60k", "efg_agilent_sureprint_g3_ge_8x60k_v2")
    print(as.data.frame(agilent.annot, stringsAsFactors = F))
    cat("Please enter a number to indicate the filter to use as input to the query...\n")
    n1 <- scan(n=1)
    attributes <- listAttributes(mart)[c(1:100),1]
    print(as.data.frame(attributes, stringsAsFactors = F))
    cat("Please enter a number to indicate the annotation to use...\n")
    n2 <- scan(n=1)
    resq <- NULL       
    if(is.numeric(n1) & is.numeric(n2)) {
      # querying BioMart
      print(c(agilent.annot[n1], attributes[n2]))
      print(agilent.annot[n1])
      print(head(probes))
      resq <- getBM(attributes = c(agilent.annot[n1], attributes[n2]), 
                    filters = agilent.annot[n1], 
                    values = probes, 
                    mart = mart)
    }
    else stop("Invalid inputs.")
    # keep only complete cases
    if(verbose) cat(" - number of probes referring to 'NA' values:", length(which(complete.cases(resq) == FALSE)),"\n", sep="")
    if(length(which(complete.cases(resq) == FALSE)) > 0) resq <- resq[complete.cases(resq),]
    # remove probes referring to different gene ids
    u.probes <- unique(resq[,1])
    tests <- unlist(lapply(u.probes, function(x) { length(resq[which(x == resq[,1]), 2]) > 1}))
    if(verbose) cat(" - number of shared probes:", length(which(tests == TRUE)),"\n", sep="")
    if(length(which(tests == TRUE)) > 0) {
      u.probes <- u.probes[which(tests == FALSE)]
      ids.to.keep <- unlist(lapply(u.probes, function(p) which(resq[,1] == p)))
      resq <- resq[ids.to.keep,]
    }
    if(verbose) cat(" - number of unique probes:", length(unique(resq[,1])), "\n", sep="")
    # check species
    testh <- grep("hsap", organism)
    if(length(testh) > 0) spec <- "Homo sapiens"
    else spec <- "##"
    print(head(resq[,2]))
    print(attributes[n2])
    annotations <- get.genes(geneids = resq[,2], mart = mart, species=spec, type = attributes[n2])
    return(list(map = resq, annot = annotations))
  }
  NULL
}
getGene <- function (id, type, mart) {
  #biomaRt:::martCheck(mart, "ensembl")
  biomaRt:::checkWrapperArgs(id, type, mart)
  symbolAttrib = switch(strsplit(biomaRt:::martDataset(mart), "_", fixed = TRUE, 
                                 useBytes = TRUE)[[1]][1], hsapiens = "hgnc_symbol", mmusculus = "mgi_symbol", "external_gene_id")
  typeAttrib = switch(type, affy_hg_u133a_2 = "affy_hg_u133a_v2", type)
  attrib = c(typeAttrib, symbolAttrib, "description", "chromosome_name", 
             "band", "strand", "start_position", "end_position", "ensembl_gene_id")
  table = biomaRt:::getBM(attributes = attrib, filters = type, values = id, mart = mart)
  return(table)
}
get.genes <- function(geneids, mart, species="Homo sapiens", type) {
  require(biomaRt)
  if(species=="Homo sapiens")
    gannot <- getGene(id = geneids, type = type, mart = mart)[,c(type, "hgnc_symbol","description","ensembl_gene_id")]
  else {
    gannot <- getGene(id = geneids, type = type, mart = mart)  #[,c(type, "symbol","description","ensembl_gene_id")]
    print(colnames(gannot))
  }
  id.dups <- which(duplicated(gannot[,type]))
  if(length(id.dups) > 0) gannot <- gannot[-id.dups,]
  print(head(gannot)) 
  gannot[,"description"] <- gsub(" \\[.*\\]", "", gannot[,"description"])
  rownames(gannot) <- gannot[,type]
  gannot
}
## GENE EXPRESSION ANALYSIS
aggreg.probes <- function(data, map, var.mds = NULL, plot = TRUE, verbose = TRUE, ...) {
  genes <- unique(map[,2])
  gene.to.probes <- lapply(genes, function(g) map[which(map[,2] == g),1])
  names(gene.to.probes) <- genes
  probes.to.dupl <- lapply(unlist(gene.to.probes), function(probe) which(rownames(data) == probe))
  names(probes.to.dupl) <- unlist(gene.to.probes)
  mat <- apply(data, 2, function(s) {
               sapply(gene.to.probes, 
                      function(probes) {
                          #print(unlist(probes.to.dupl[probes]))
                          median(s[unlist(probes.to.dupl[probes])])
                      })
               })
  if(verbose) cat("Number of annotated genes: ", nrow(mat), "\n", sep="")
  if(plot & !is.null(var.mds)) 
    limma:::plotMDS(mat, top=1000, labels=var.mds, col=as.numeric(var.mds), 
            gene.selection="common", main = "After mapping probes to genes")
  mat
}
aggreg.probes.2 <- function(data, map, var.mds = NULL, plot = TRUE, verbose = TRUE, ...) {
  #Update map column 2 to TargetID
  colnames(map)[2] <- "TargetID"
  #Aggregate dup probes before aggregating by annotation IDs
  tmpRowNames <- rownames(data)
  dups <- which(duplicated(tmpRowNames))
  if(length(dups)>0){
    dupRowNames <- unique(tmpRowNames[dups])
    unique.data <- data[-which(tmpRowNames %in% dupRowNames),]
    colnames(unique.data) <- make.names(colnames(unique.data))
    dup.data <- data[which(tmpRowNames %in% dupRowNames),]
    dup.data <- data.frame(ProbeID=rownames(dup.data), dup.data, row.names=NULL)
    dup.datag <- aggregate(.~ProbeID, data = dup.data, median)
    rownames(dup.datag) <- dup.datag[,1]
    dup.datag <- dup.datag[,-1]
    data <- rbind(unique.data, dup.datag)
  }
  rownames(map) <- map[,1]
  map <- map[rownames(data),]
  map <- cbind(map,data)
  map <- map[,-1]
  blankIdx <- which(map$TargetID=="")
  if(length(blankIdx)>0){
    map <- map[-blankIdx,]
  }
  if(verbose)  print(str(map))
  datag = aggregate(.~TargetID, data = map, median)
  rownames(datag) <- datag[,1]
  datag <- datag[,-1]
  if(plot & !is.null(var.mds))  
    limma:::plotMDS(datag, top=1000, labels=var.mds, col=as.numeric(var.mds), 
                    gene.selection="common", main = "After mapping probes to genes")
  datag
}
diff.gene.expr <- function(data, des, contrasts, pvalue, fcvalue, p.adjust.method, annot = NULL, max.ngenes=Inf, save.file = NULL, plot = TRUE, verbose = TRUE, ...) {
  # make the matrix of contrasts
  colnames(des) <- make.names(colnames(des))
  cont <- makeContrasts(contrasts=contrasts, levels=des)
  if(verbose) {
    print(des)
    print(cont)
  }
  fit <- lmFit(data, des)
  ##Added the print for contrast matrix
  #print(contrasts.fit(fit, cont))
  fit2 <- eBayes(contrasts.fit(fit, cont))
  #max.ngenes <- read.input.number(" - indicate the max number of genes:")
  if(plot & length(contrasts) < 6) {
    vennDiagram(decideTests(fit2, p.value=pvalue, lfc=log2(fcvalue), adjust.method=p.adjust.method),
                cex=c(.8,.8,0.7),
                main = paste("Venn Diagram <p.value=", pvalue, ",FC=", fcvalue, ">",sep=""))
    #vennDiagram(decideTests(fit2, p.value=pvalue, lfc=log2(fcvalue), adjust.method=methods[n1]), include=c("up","down"), counts.col=c("red","green"),
    #            main = paste("Venn Diagram <p.value=", pvalue, ",FC=", fcvalue, ">",sep=""))
    ##for(i in 1:length(vars$contrasts))  volcanoplot(fit2, coef=i, highlight=50)
  }
  # to output a table of the differentially expressed genes
  list.top.tables <- list()
  if(verbose)  cat(" - generating top tables.. \n", sep="")
  for(i in 1:length(contrasts)) {
    list.top.tables[[i]] <- topTable(fit2, coef=i, p.value=pvalue, lfc=log2(fcvalue), 
                                     adjust.method=p.adjust.method,
                                     number=max.ngenes, sort.by="logFC")
    if(verbose) {
      cat("     > ", contrasts[i], " - number of sig. genes:",  nrow(list.top.tables[[i]]), "\n", sep="")
    }
  }
  names(list.top.tables) <- contrasts

  # how to order the gene expression data
  for(i in 1:length(list.top.tables)) {
    ids <- NULL
    if(dim(list.top.tables[[i]])[1] == 0) next
    if(!is.null(annot)) {
      list.top.tables[[i]] <- base:::merge(list.top.tables[[i]], annot, by = "row.names")
      rownames(list.top.tables[[i]]) <- list.top.tables[[i]]$Row.names
      list.top.tables[[i]] <- list.top.tables[[i]][,-1]
      #list.top.tables[[i]]$score <- abs(list.top.tables[[i]]$logFC) * -log(list.top.tables[[i]]$P.Value)
    }else{
      ids <- rownames(list.top.tables[[i]])
    }
    list.top.tables[[i]]$score <- list.top.tables[[i]]$logFC * -log(list.top.tables[[i]]$P.Value)
    if(!is.null(ids)) list.top.tables[[i]]$ID <- rownames(list.top.tables[[i]])
    colnames(list.top.tables[[i]])[which(colnames(list.top.tables[[i]]) == "t")] <- "t-statistic" 
    colnames(list.top.tables[[i]])[which(colnames(list.top.tables[[i]]) == "B")] <- "B-statistic"
  }
  list.top.tables
}
volcano.plot <- function(data, cutoff.pvalue = .05, cutoff.lfc = 1) {
  require(calibrate)
  plot(data$logFC, -log10(data$P.Value), pch=20, main="Volcano plot", xlim=c(-2.5,2))
  pdata <- subset(data, P.Value < cutoff.pvalue) 
  points(pdata$logFC, -log10(pdata$P.Value), pch=20, col="gray")
  fdata <- subset(data, abs(logFC) > cutoff.lfc)
  points(fdata$logFC, -log10(fdata$P.Value), pch=20, col="orange")
  gdata <- subset(data, P.Value < cutoff.pvalue & abs(logFC) > cutoff.lfc)
  points(gdata$logFC, -log10(gdata$P.Value), pch=20, col="blue")
  textxy(gdata$logFC, -log10(gdata$P.Value), labs=gdata$hgnc_symbol, cex=.6)
}
save.xls.file <- function(lists.diff.genes = NULL, file.name = "diff_genes.xlsx") {
  #require(xlsx)
  require(XLConnect)
  options(java.parameters = "-Xmx4g")
  if(!is.null(lists.diff.genes)) {
    gc()
    # creating work book
    # wb <- createWorkbook()
    xlcFreeMemory()
    wb <- loadWorkbook(file.name, create=TRUE)
    for(i in 1:length(lists.diff.genes)) {
      # making each as a sheet
      #sheet <- createSheet(wb, sheetName=names(lists.diff.genes)[i])
      sheet <- names(lists.diff.genes)[i]
      createSheet(wb, sheet)
      #addDataFrame(lists.diff.genes[[i]], sheet)
      writeWorksheet(wb, lists.diff.genes[[i]], sheet, startRow=1, startCol=1, header=TRUE, rownames=rownames(lists.diff.genes[[i]]))
    }
    # saving the workbook
    #saveWorkbook(wb, file.name)
    saveWorkbook(wb)
  }
  NULL
}
david.annot <- function(list.top.tables) {
  require(BACA)
  require(ggplot2)
  split <- read.input.number("Enter [1] to split each list in two sublists: down- and up-regulated genes.")
  # building the gene list
  if(split == 1) {
    gene.lists <- unlist(lapply(list.top.tables, function(d) { 
      down.genes <- rownames(d)[which(d$logFC < 0)]
      if(length(down.genes)==0) down.genes <- "1"
      up.genes <- rownames(d)[which(d$logFC >= 0)]  
      if(length(up.genes)==0) up.genes <- "1"
      list(down.genes, up.genes)}), recursive = FALSE)
    # naming the gene lists 
    names(gene.lists) <- unlist(lapply(names(list.top.tables), function(x) paste(x,c(1,2),sep="_")))
    print(length(gene.lists))
  }
  else {
    gene.lists <- lapply(list.top.tables, rownames)
    names(gene.lists) <- names(list.top.tables)
  }
  return(gene.lists)
  # querying DAVID for KEGG.PATHWAYS
  kegg.pathways <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                               idType = "ENTREZ_GENE_ID", annotation = "KEGG_PATHWAY")

  bp.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                               idType = "ENTREZ_GENE_ID", annotation = "GOTERM_BP_ALL")
  
  mf.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                               idType = "ENTREZ_GENE_ID", annotation = "GOTERM_MF_ALL")
  
  cc.terms <- DAVIDsearch(gene.lists, david.user = "vittorio.fortino@ttl.fi", 
                          idType = "ENTREZ_GENE_ID", annotation = "GOTERM_CC_ALL")
  
  list(show.bbplot(kegg.pathways, "BBplot - KEGG", names(list.top.tables)),
       show.bbplot(bp.terms, "BBplot - GO terms (BP)", names(list.top.tables)),
       show.bbplot(mf.terms, "BBplot - GO terms (MF)", names(list.top.tables)),
       show.bbplot(cc.terms, "BBplot - GO terms (CC)", names(list.top.tables)))
  
} 
show.bbplot <- function(res.david.query, title.query, col.names, verbose = TRUE) {
  cont.it <- TRUE
  if(verbose)  cat("Building the bbplot..", title.query, "\n", sep="")
  while(cont.it) {
    max.pval <- read.input.number("Indicate the p.value:")
    min.ngenes <- read.input.number("Indicate the min number of genes:")
    cat("Do you want to specify the max number of genes?", "\n")
    if(continue.read.inputs()) {
      max.ngenes <- read.input.number("Indicate the max number of genes:")
      bbplot <- BBplot(res.david.query, max.pval = max.pval, 
                       min.ngenes = min.ngenes, 
                       max.ngenes = max.ngenes,
                       adj.method = "",
                       name.com = col.names, 
                       title = paste(title.query, "<p.value=", max.pval,
                                     ";min.num.genes=", min.ngenes,
                                     ";max.num.genes=", max.ngenes,
                                     ">", sep=""), print.term = "full") 
      print(bbplot)
    }
    else {
      # plotting
      bbplot <- BBplot(res.david.query, max.pval = max.pval, min.ngenes = min.ngenes, name.com = col.names, 
                       title = paste(title.query, "<p.value=", max.pval,";min.num.genes=", min.ngenes,">", sep=""), print.term = "full") 
      print(bbplot)
    }
    cont.it <- continue.read.inputs()
  }
  bbplot
}
## USER INPUTS FUNCTIONS
read.input.word <- function(n=1, message = "Indicate a number", more.info=NULL) {
  cat(message, "\n", sep="")
  if(!is.null(more.info)) cat(" *** ", more.info, "\n")
  keywords <- scan(what=character(), nlines=n)
  paste(keywords, collapse=",")
}
read.input.number <- function(message = "Indicate a number", more.info=NULL) {
  cat(message, "\n", sep="")
  if(!is.null(more.info)) cat(more.info, "\n", sep="")
  n <- NULL
  while(!is.numeric(n)) {
    n <- scan(n=1)
    if(!is.numeric(n))
      cat("Invalid input number.", "\n")
  }
  n
}
continue.read.inputs <- function(message="Press [y] to continue, [n].") {
  next.it <- ""
  while(next.it != "y" &  next.it != "n") {
    cat(message)
    next.it <- read.input.word(message = "")
  }
  if(next.it == "y") TRUE
  else FALSE
} 
## UTILITY functions
droplevels.annot <- function(a, remove = F) {
  id.col.to.remove <- NULL
  for(j in 1:ncol(a)) {
    if(is.character(a[,j]))
      a[,j] <- as.factor(a[,j])
    if(is.factor(a[,j])) {
      a[,j] <- droplevels(a[,j])
    }
    if(length(levels(a[,j])) > 1)
      print(colnames(a)[j])
    else
      id.col.to.remove <- c(id.col.to.remove, j)
  }
  print(j)
  if(remove) 
    a[,-j]
  else
    a
}
build.model.matrix <- function(pd, intercept = -1, var.int, covariates = NULL, verbose = TRUE) {
  des <- NULL
  if(!is.null(covariates)) {
    f <- formula(paste("~", intercept, "+pd$", var.int, "+pd$", paste(covariates, collapse="+pd$"), sep=""))
    if(verbose) print(f)
    des <- model.matrix(f)
    colnames(des) <- gsub(paste("pd\\$", var.int, sep=""), "", colnames(des))
    for(i in 1:length(covariates))
      colnames(des) <- gsub( paste("pd\\$", covariates[i], sep=""), "", colnames(des))
  }
  else {
    f <- formula(paste("~", intercept, "+pd$", var.int, sep=""))
    if(verbose) print(f)
    des <- model.matrix(f)
    colnames(des) <- levels(as.factor(pd[,var.int]))
  }
  des
}
color.map <- function(labels) {
  levs <- levels(labels)
  colors <- rep(NA, length(labels))
  for(i in 1:length(levs)) {
      colors[which(labels == levs[i])] <- COLORS[i]
  }
  colors
}
hclust2 <- function(x, method="ward.D2", ...) hclust(x, method=method, ...)
dist2 <- function(x, ...) as.dist(1-cor(t(x), method="pearson"))
plot.heatmap <- function(data, labels, top = 500) {
  sidebarcolors <- color.map(labels)
  hvars <- apply(data, 1, var,na.rm=TRUE) 
  hvars <- sort(hvars, decreasing=TRUE) 
  hvars <- hvars[1:top] 
  print(dim(data[names(hvars),]))
  heatmap.2(data[names(hvars),], Rowv=TRUE, 
            scale="column", trace="none", 
            distfun=dist2, hclustfun = hclust2,
            col=redgreen, xlab=NULL, ylab=NULL, 
            labRow = "", dendrogram = "column", cexCol = .7,
            ColSideColors=sidebarcolors)
  par(lend = 1)           # square line ends for the color legend
  legend("bottomleft",      # location of the legend on the heatmap plot
         legend = levels(labels), # category labels
         col = COLORS[1:length(levels(labels))],  # color key
         bty = "n",
         bg = "white",
         border = "white",
         cex = .7,
         lty= 1,            # line style
         lwd = 10           # line width
  )     
}
plot.heatmpa.genes <- function(expr.dat, col, sf = 3, ze = 0.4, pi = 0.4, o = 1, ccolors) {
  require(Biobase)
  require(DFP)
  dfp <- featSelDFP(expr.dat, skipFactor = sf, zeta = ze, piVal = pi, overlapping = o)[[3]]
  dl <- as.matrix(dfp)
  dl <- apply(dl, MARGIN = 2, FUN = function(x) {
    id.na <- which(is.na(x))
    if(length(id.na) > 0) x[id.na] <- "Not assigned"
    x
  })
  feats <- rownames(dfp)
  dat <- expand.grid(var1=1:nrow(dfp), var2=1:ncol(dfp))
  lev <- colnames(attr(dfp,"ifs"))
  dat$var2 <- sapply(dat$var2, function(i) lev[i])
  dat$value <- round(melt(attr(dfp,"ifs"))$value, 2)
  dat$labels <- melt(dl)$value
  print(dat)
  ggplot(dat, aes(x=var1,y=var2))+ 
    geom_point(aes(size = value, colour = labels)) +
    scale_colour_manual(values=ccolors) +
    scale_size(range = c(5, 10), breaks=c(0.25,.5,.75),labels=c("25%","50%","75%")) +
    #scale_size(range = c(5, 10)) + 
    scale_x_continuous(breaks = 1:nrow(dfp), labels = feats) +
    scale_y_discrete(limits = sort(levels(expr.dat$class), decreasing = T)) + 
    labs(size = "Coverage", colour="Gene status") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          panel.background = element_blank(),
          legend.text = element_text(size = 10),
          legend.position="right") +
    xlab("Genes") +
    ylab("Class Labels") 
}
split.top.tables <- function(list.tables, fc = 1.5) { 
  gene.lists <- unlist(lapply(list.tables, function(d) { 
    rownames(d) <- sub("_at",'', rownames(d))
    down.genes <- rownames(d)[which(d$logFC < -log2(fc))]
    if(length(down.genes)==0) down.genes <- "1"
    up.genes <- rownames(d)[which(d$logFC >= log2(fc))]  
    if(length(up.genes)==0) up.genes <- "1"
    list(down.genes, up.genes)}), recursive = FALSE)
  # naming the gene lists 
  names(gene.lists) <- unlist(lapply(names(list.tables), function(x) paste(x,c("down","up"),sep="_")))
  print(length(gene.lists))
  print(lapply(gene.lists, length))
  gene.lists
}
# Goodman and Kruskal's tau measure
GKtau <- function(x,y){
  #  First, compute the IxJ contingency table between x and y
  Nij = table(x,y,useNA="ifany")
  #  Next, convert this table into a joint probability estimate
  PIij = Nij/sum(Nij)
  #  Compute the marginal probability estimates
  PIiPlus = apply(PIij,MARGIN=1,sum)
  PIPlusj = apply(PIij,MARGIN=2,sum)
  #  Compute the marginal variation of y
  Vy = 1 - sum(PIPlusj^2)
  #  Compute the expected conditional variation of y given x
  InnerSum = apply(PIij^2,MARGIN=1,sum)
  VyBarx = 1 - sum(InnerSum/PIiPlus)
  #  Compute and return Goodman and Kruskal's tau measure
  tau = (Vy - VyBarx)/Vy
  tau
}
# DFP-based analysis
featSelDFP <- function(exprData, skipFactor = 3, zeta = 0.5, piVal = 0.9, overlapping = 2) {
  require(DFP)
  mfs <- calculateMembershipFunctions(exprData, skipFactor)
  print("calculateMembershipFunctions <DONE>")
  dvs <- discretizeExpressionValues(exprData, mfs, zeta, overlapping)
  print("discretizeExpressionValues <DONE>")
  fps <- calculateFuzzyPatterns(exprData, dvs, piVal, overlapping)
  print("calculateFuzzyPatterns <DONE>")
  list.fps <- lapply(names(table(exprData$class)), function(c) { 
    fp <- showFuzzyPatterns(fps, c) 
    rownames(exprData)[which((rownames(exprData) %in% names(fp)) == TRUE)]
  })
  #xlistFPs <- lapply(names(table(exprData$class)), FUN= function(c) { fp = showFuzzyPatterns(fps, c); which( (rownames(exprData) %in% names(fp)) == TRUE)})
  # print(xlistFPs)
  dfps <- calculateDiscriminantFuzzyPattern(exprData, fps)
  #plotDiscriminantFuzzyPattern(dfps, overlapping)
  return(list(fps,list.fps,dfps,dvs))
}
merge.combat <- function (esets, covariates = NULL, verbose = T)  {
  raw_merged = merge.array(esets)
  batchInfo = NULL
  for (i in 1:length(esets)) {
    batchInfo = c(batchInfo, rep(i, ncol(esets[[i]])))
  }
  #mod = cbind(as.factor(imp_ann$Nano2), as.factor(imp_ann$TimeExposition),as.factor(imp_ann$cell_type2))
  #print(colnames(pData(raw_merged)))
  c.names <- c("Array name", "Sample name", "Batch")
  saminfo = cbind(rownames(pData(raw_merged)), rownames(pData(raw_merged)), batchInfo)
  if(!is.null(covariates)) {
    for(i in 1:length(covariates)) {
      saminfo <- cbind(saminfo, as.factor(pData(raw_merged)[,covariates[i]]))
      c.names <- c(c.names, covariates[i])
    }
  }
  #colnames(saminfo) = c("Array name", "Sample name", "Batch")
  colnames(saminfo) = c.names
  dat = exprs(raw_merged)
  design <- design.mat(saminfo)
  if(verbose) print(design)
  batches <- list.batch(saminfo)
  n.batch <- length(batches)
  n.batches <- sapply(batches, length)
  n.array <- sum(n.batches)
  B.hat <- solve(t(design) %*% design) %*% t(design) %*% t(as.matrix(dat))
  grand.mean <- t(n.batches/n.array) %*% B.hat[1:n.batch, ]
  var.pooled <- ((dat - t(design %*% B.hat))^2) %*% rep(1/n.array, 
                                                        n.array)
  stand.mean <- t(grand.mean) %*% t(rep(1, n.array))
  if (!is.null(design)) {
    tmp <- design
    tmp[, c(1:n.batch)] <- 0
    stand.mean <- stand.mean + t(tmp %*% B.hat)
  }
  s.data <- (dat - stand.mean)/(sqrt(var.pooled) %*% t(rep(1, n.array)))
  batch.design <- design[, 1:n.batch]
  gamma.hat <- solve(t(batch.design) %*% batch.design) %*% 
    t(batch.design) %*% t(as.matrix(s.data))
  delta.hat <- NULL
  for (i in batches) {
    delta.hat <- rbind(delta.hat, apply(s.data[, i], 1, var, 
                                        na.rm = T))
  }
  gamma.bar <- apply(gamma.hat, 1, mean)
  t2 <- apply(gamma.hat, 1, var)
  a.prior <- apply(delta.hat, 1, aprior)
  b.prior <- apply(delta.hat, 1, bprior)
  gamma.star <- delta.star <- NULL
  for (i in 1:n.batch) {
    temp <- it.sol(s.data[, batches[[i]]], gamma.hat[i, ], 
                   delta.hat[i, ], gamma.bar[i], t2[i], a.prior[i], 
                   b.prior[i])
    gamma.star <- rbind(gamma.star, temp[1, ])
    delta.star <- rbind(delta.star, temp[2, ])
  }
  bayesdata <- s.data
  j <- 1
  for (i in batches) {
    bayesdata[, i] <- (bayesdata[, i] - t(batch.design[i, 
                                                       ] %*% gamma.star))/(sqrt(delta.star[j, ]) %*% t(rep(1, 
                                                                                                           n.batches[j])))
    j <- j + 1
  }
  bayesdata <- (bayesdata * (sqrt(var.pooled) %*% t(rep(1, 
                                                        n.array)))) + stand.mean
  eset = raw_merged
  exprs(eset) = bayesdata
  return(eset)
}
merge.array <- function (esets)  {
  eset1 = esets[[1]]
  annot1 = annotation(eset1)
  for (i in 2:length(esets)) {
    eset2 = esets[[i]]
    d1 = exprs(eset1)
    d2 = exprs(eset2)
    cg = sort(intersect(rownames(d1), rownames(d2)))
    if (length(cg) < (min(dim(d1)[1], dim(d2)[1])/100)) {
      msg(" ! WARNING ! Number of common genes < 1%")
    }
    fData = fData(eset1)[cg, ]
    p1 = pData(eset1)
    p2 = pData(eset2)
    cp = sort(intersect(colnames(p1), colnames(p2)))
    tp = sort(unique(union(colnames(p1), colnames(p2))))
    sp1 = setdiff(colnames(p1), cp)
    sp2 = setdiff(colnames(p2), cp)
    pheno = matrix(NA, ncol = length(tp), nrow = nrow(p1) + 
                     nrow(p2))
    rownames(pheno) = c(rownames(p1), rownames(p2))
    colnames(pheno) = tp
    if (length(cp) != 0) {
      pheno[1:nrow(p1), cp] = as.matrix(p1[, cp])
      pheno[(nrow(p1) + 1):(nrow(p1) + nrow(p2)), cp] = as.matrix(p2[, 
                                                                     cp])
    }
    if (length(sp1) != 0) {
      pheno[1:nrow(p1), sp1] = as.matrix(p1[, sp1])
    }
    if (length(sp2) != 0) {
      pheno[(nrow(p1) + 1):(nrow(p1) + nrow(p2)), sp2] = as.matrix(p2[, 
                                                                      sp2])
    }
    pData = as.data.frame(pheno)
    d1 = d1[cg, , drop = FALSE]
    d2 = d2[cg, , drop = FALSE]
    eset1 = new("ExpressionSet")
    exprs(eset1) = cbind(d1, d2)
    pData(eset1) = pData
    fData(eset1) = fData
    annot1 = c(annot1, annotation(eset2))
  }
  annotation(eset1) = unique(annot1)
  return(eset1)
}
design.mat <- function (saminfo, verbose = TRUE) {
  tmp <- which(colnames(saminfo) == "Batch")
  tmp1 <- as.factor(saminfo[, tmp])
  if(verbose) cat("  => Found", nlevels(tmp1), "batches")
  design <- build.design(tmp1, start = 1)
  ncov <- ncol(as.matrix(saminfo[, -c(1:2, tmp)]))
  if(verbose) cat("  => Found", ncov, "covariate(s)")
  if (ncov > 0) {
    for (j in 1:ncov) {
      tmp1 <- as.factor(as.matrix(saminfo[, -c(1:2, tmp)])[, 
                                                           j])
      design <- build.design(tmp1, des = design)
    }
  }
  design
}
build.design <- function (vec, des = NULL, start = 2) {
  tmp <- matrix(0, length(vec), nlevels(vec) - start + 1)
  for (i in 1:ncol(tmp)) {
    tmp[, i] <- vec == levels(vec)[i + start - 1]
  }
  cbind(des, tmp)
}
list.batch <- function (saminfo) {
  tmp1 <- as.factor(saminfo[, which(colnames(saminfo) == "Batch")])
  batches <- NULL
  for (i in 1:nlevels(tmp1)) {
    batches <- append(batches, list((1:length(tmp1))[tmp1 == levels(tmp1)[i]]))
  }
  batches
}
aprior <- function (gamma.hat){
  m = mean(gamma.hat)
  s2 = var(gamma.hat)
  (2 * s2 + m^2)/s2
}
bprior <- function (gamma.hat) {
  m = mean(gamma.hat)
  s2 = var(gamma.hat)
  (m * s2 + m^3)/s2
}
it.sol <- function (sdat, g.hat, d.hat, g.bar, t2, a, b, conv = 1e-04)  {
  n <- apply(!is.na(sdat), 1, sum)
  g.old <- g.hat
  d.old <- d.hat
  change <- 1
  count <- 0
  while (change > conv) {
    g.new <- postmean(g.hat, g.bar, n, d.old, t2)
    sum2 <- apply((sdat - g.new %*% t(rep(1, ncol(sdat))))^2, 
                  1, sum, na.rm = T)
    d.new <- postvar(sum2, n, a, b)
    change <- max(abs(g.new - g.old)/g.old, abs(d.new - d.old)/d.old)
    g.old <- g.new
    d.old <- d.new
    count <- count + 1
  }
  adjust <- rbind(g.new, d.new)
  rownames(adjust) <- c("g.star", "d.star")
  adjust
}
postmean <- function (g.hat, g.bar, n, d.star, t2)  {
  (t2 * n * g.hat + d.star * g.bar)/(t2 * n + d.star)
}
postvar <- function (sum2, n, a, b)  {
  (0.5 * sum2 + b)/(n/2 + a - 1)
}

## MDS plot
plotMDS <- function (eset, colLabel, symLabel, legend = TRUE, file = NULL,  ctr = FALSE, md = FALSE, ...) {
  if (!is.null(file)) {
    pdf(file, width = 12, height = 7)
  }
  mds = cmdscale(dist(t(exprs(eset))), eig = TRUE)
  colMap = makeColorMap(eset, colLabel)
  colVec = makeColorVec(eset, colLabel, colMap)
  tmp = par()$mar
  if (legend) {
    par(xpd = T, mar = par()$mar + c(0, 0, 0, 4))
  }
  range_x = range(mds$points[, 1])
  range_y = range(mds$points[, 2])
  plot(mds$points, col = colVec, pch = as.numeric(pData(eset)[, 
                                                              symLabel]), panel.first = {
                                                                U = par("usr")
                                                                rect(U[1], U[3], U[2], U[4], col = "azure2", border = "black", 
                                                                     lwd = 3)
                                                              }, lwd = 2, xlab = "", ylab = "", xlim = range_x, ylim = range_y, 
       ...)
  if(ctr) {
    rect(-.5, -.5, .5, .5, border = "black", lty = 3, lwd = 2)
    rect(-1.5, -1.5, 1.5, 1.5, border = "black", lty = 2, lwd = 2)
    rect(-2.5, -2.5, 2.5, 2.5, border = "black", lwd = 2)
  }
  if (legend) {
    x = range_x[2] + (range_x[2] - range_x[1]) * 0.05
    y = range_y[2] - (range_y[2] - range_y[1]) * 0.05
    syms = unique(pData(eset)[, symLabel])
    legend("topleft", legend = syms, pt.lwd = 2, cex = 0.7, pch = as.numeric(syms), 
           box.lwd = 3, bg = "azure2")
    legend("topright", inset=c(-0.3,0), legend = names(colMap), pt.lwd = 2, pch = 19, 
           col = unlist(colMap), box.lwd = 3, bg = "azure2")
    par(xpd = F, mar = tmp)
  }
  if (!is.null(file)) {
    dev.off()
  }
  if(md) return(mds)
}
makeColorMap <- function (eset, label){
  colMap = list()
  vec = unique(as.vector(pData(eset)[, label]))
  for (i in 1:length(vec)) {
    colMap[[vec[i]]] = COLORS[i]
  }
  return(colMap)
}
makeColorVec <- function (eset, label, colMap) {
  labels = as.vector(pData(eset)[, label])
  return(as.vector(unlist(sapply(labels, function(x) {
    colMap[x]
  }))))
}
COLORS <- c("green3","blue","cyan","magenta","yellow","gray","red","orange",
            "darkred", "green","darkblue","darkcyan","darkmagenta","darkorchid1",
            "darkgoldenrod3","aquamarine","darkslategray3","darkolivegreen3","lightcoral",
            "deeppink","gold","Olivedrab1","dimgrey","cornflowerblue","darkgreen",
            "burlywood3","steelblue4","orangered","purple1","khaki1","azure4",
            "blue1","blue2","blue3","blue4","coral1","coral2","coral3","coral4")

#Get Summary of Differential Expression Tables
get_deg_summary <- function(deg_list, names, lfc=0){
        cuts <- c(0.0001, 0.001, 0.01, 0.025, 0.05, 0.10, 1)
        cat("\nStatistical significance summary:\n")
        countsDF <- NULL
        deg_list <- as.list(deg_list)

        for(i in 1:length(deg_list)){
                objectName <- names(deg_list)[i]
                object <- deg_list[[objectName]]
                idx<-which(abs(object$logFC)>lfc); 
                if(length(idx)>0){
                        object <- object[idx,]
                        counts <- sapply(cuts, function(x) c("P.Value"=sum(object$P.Value < x), "adj.P.Val"=sum(object$adj.P.Val < x)))
                }else{
                        counts <- data.frame(sapply(cuts, function(x){c(0,0)}))
                        colnames(counts) <- cuts
                        rownames(counts) <- c("P.Value", "adj.P.Val")
                }
                colnames(counts) <- paste("<", cuts, sep="")
                countName <- names[i]
                #rownames(counts) <- paste0(countName, ";", rownames(counts))
                ##print(counts)
                #countsDF <- rbind(countsDF, counts)
                print(class(counts))
                tmpDF <- data.frame(comparison=countName, score=rownames(counts), counts, stringsAsFactors=FALSE)
                colnames(tmpDF)[3:ncol(tmpDF)] <- paste("<", cuts, sep="")
                rownames(tmpDF) <- paste0(countName, ";", rownames(counts))
                countsDF <- rbind(countsDF, tmpDF)
                print("Summary Table Class:")
                print(class(countsDF))
        }
        return(countsDF)
}

#Get color palette for factor
get_color_palette <- function(iVec, asFactor=FALSE){
        set.seed(1) #Block randomness. Set seed
        print("Getting color palette..")
        if(asFactor){
                print("Input as factor...")
                if(is.factor(iVec)){
                        print("Input is factor...")
                        nm <- levels(iVec)
                }else if(is.vector(iVec)){
                        print("Input is vector...")
                        iVec <- factor(iVec)
                        nm <- levels(iVec)
                }else{
                        print("Neither factor nor vector, returning NULL!")
                        return(NULL)
                }
                print("levels:")
                print(nm)
                #nmPalette <- setNames(randomcoloR::distinctColorPalette(length(nm)), seq(nm))
                nmPalette <- setNames(randomcoloR::randomColor(length(nm), luminosity="dark"), seq(nm))
                print("Palette base:")
                print(nmPalette)
                print("Levels as integer:")
                print(as.integer(iVec))
                colorVec <- nmPalette[as.integer(iVec)]
        }else{
                print("Input as vector...")
                if(is.factor(iVec)){
                        print("Input is factor...")
                        iVec <- as.character(iVec)
                }
                print("iVec:")
                print(iVec)
                #colorVec <- setNames(randomcoloR::distinctColorPalette(length(iVec)), iVec)
                colorVec <- setNames(randomcoloR::randomColor(length(iVec), luminosity="dark"), iVec)
        }
        print("Palette vector:")
        print(colorVec)
        return(colorVec)
}

#Get array columns
get_array_cols <- function(sourceType){
        columns <- switch(sourceType, 
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
        return(columns)
}

#Check the chosen source against raw data file
check_source_type <- function(rawDir, fileName, columns){
        Sys.setlocale('LC_ALL','C')
        print("Checking Source:")
        file <- file.path(rawDir, fileName)
        con <- file(file, "r")
        Found <- FALSE
        i <- 0
        repeat {
                i <- i + 1
                txt <- readLines(con, n = 1)
                #if(i<=35){
                #        print(txt)
                #}
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
        return(Found)
}

