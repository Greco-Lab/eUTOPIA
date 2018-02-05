# eUTOPIA
### A new tool for Inference of NetwOrk Response Modules

A novel computational method and its R and web-based implementations, to perform inference of gene network from transcriptome data and prioritization of key genes with central functional and topological role in the network.

#### Install Dependencies
```R
  #Install CRAN dependencies
  cran_pkgs <- c("V8", "RSQLite", "TopKLists", "doParallel", "foreach", "igraph", "plyr", "shiny", "shinyjs", "shinyBS", "shinydashboard", "colourpicker", "DT", "R.utils", "treemap", "visNetwork", "abind", "radarchart", "randomcoloR", "Rserve", "WriteXLS", "gplots", "ggplot2")
  cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
  if(length(cran_pkgs.inst)>0){
    print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
    for(pkg in cran_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
      print("Installed!!!")
    }
  }
  
  #Install Bioconductor dependencies
  if(!"GOSemSim" %in% rownames(installed.packages())){
    print("Installing GOSemSim from GitHub!")
    devtools::install_github("GuangchuangYu/GOSemSim")
  }
  
  source("http://bioconductor.org/biocLite.R")
  bioc_pkgs <- c("org.Hs.eg.db", "org.Mm.eg.db", "GO.db", "AnnotationDbi", "GSEABase", "minet")
  bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
  if(length(bioc_pkgs.inst)>0){
    source("http://bioconductor.org/biocLite.R")
    print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
    for(pkg in bioc_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      biocLite(pkg, suppressUpdates=TRUE)
      print("Installed!!!")
    }
  }
```

#### How to run INfORM from GitHub
```R
  # Load 'shiny' library
  library(shiny)

  # Using runGitHub
  runGitHub("INfORM", "Greco-Lab", subdir="INfORM-app")

  # Using the archived file
  runUrl("https://github.com/Greco-Lab/INfORM/archive/master.tar.gz", subdir="INfORM-app")
  runUrl("https://github.com/Greco-Lab/INfORM/archive/master.zip", subdir="INfORM-app")
```

#### How to run locally
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/INfORM INfORM_clone

  # Start R session and run by using runApp()
  setwd("./INfORM_clone")
  runApp("INfORM-app/")
```
#### Dependencies and License
Please refer to the 'DESCRIPTION' and 'NAMESPACE' files for information about the license and dependencies required to run INfORM.
