# eUTOPIA
### A solUTion for Omics data PreprocessIng and Analysis

Graphically accessible guided workflow for preprocessing and analysis of omics data. Supports Agilent 2-color, Agilent 1-color, Affymetrix, and Illumina methylation microarray platforms (Ongoing efforts to add support for RNA-Seq data). Discreetly separated steps in analysis designed in R Shiny, incorporates widely used microarray analysis practices and R packages. Reporting is and data interpretation is leveraged from dynamically generated plots. 

#### Install Dependencies
```R
  #Install CRAN dependencies
  cran_pkgs <- c("swamp", "infotheo", "gplots", "RColorBrewer", "shiny", "shinyjs", "shinyBS", "shinydashboard", "shinyFiles",
  "DT", "shinycssloaders", "ggplot2", "ggrepel", "WriteXLS", "rmarkdown", "VennDiagram", "grid", "futile.logger", "reshape2",
  "htmlTable", "devtools", "httr", "randomcoloR")
  cran_pkgs.inst <- cran_pkgs[!(cran_pkgs %in% rownames(installed.packages()))]
  if(length(cran_pkgs.inst)>0){
    print(paste0("Missing ", length(cran_pkgs.inst), " CRAN Packages:"))
    for(pkg in cran_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      install.packages(pkg, repo="http://cran.rstudio.org", dependencies=TRUE)
      print("Installed!!!")
    }
  }
  
  #Install latest version of rhandsontable from GitHub
  print("Installing rhandsontable from GitHub!")
  devtools::install_github("jrowen/rhandsontable")
  
  #Install latest version of UpSetR from GitHub
  print("Installing UpSetR from GitHub!")
  devtools::install_github("hms-dbmi/UpSetR")
  
  #Install Bioconductor dependencies
  #Install latest version of GOSemSim from GitHub
  print("Installing GOSemSim from GitHub!")
  devtools::install_github("GuangchuangYu/GOSemSim")
  
  source("http://bioconductor.org/biocLite.R")
  bioc_pkgs <- c("limma", "sva", "Biobase", "biomaRt", "affy", "affyQCReport", "arrayQualityMetrics", "made4", "vsn", "minfi", 
  "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19", "IlluminaHumanMethylationEPICmanifest", 
  "IlluminaHumanMethylationEPICanno.ilm10b2.hg19", "affyio", "simpleaffy", "yaqcaffy", "GO.db")
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
  runGitHub("eUTOPIA", "Greco-Lab", subdir="eUTOPIA-app")

  # Using the archived file
  runUrl("https://github.com/Greco-Lab/eUTOPIA/archive/master.tar.gz", subdir="eUTOPIA-app")
  runUrl("https://github.com/Greco-Lab/eUTOPIA/archive/master.zip", subdir="eUTOPIA-app")
```

#### How to run locally
```R
  # Clone the git repository
  git clone https://github.com/Greco-Lab/eUTOPIA eUTOPIA_clone

  # Start R session and run by using runApp()
  setwd("./eUTOPIA_clone")
  library(shiny)
  runApp("eUTOPIA-app/")
```
