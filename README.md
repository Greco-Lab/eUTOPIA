# eUTOPIA
### A solUTion for Omics data PreprocessIng and Analysis

Graphically accessible guided workflow for preprocessing and analysis of omics data. Supports Agilent 2-color, Agilent 1-color, Affymetrix, and Illumina methylation microarray platforms (Ongoing efforts to add support for RNA-Seq data). Discreetly separated steps in analysis designed in R Shiny, incorporates widely used microarray analysis practices and R packages. Reporting is and data interpretation is leveraged from dynamically generated plots.

Reference Paper:
> Marwah, V. S., Scala, G., Kinaret, P. A. S., Serra, A., Alenius, H., Fortino, V., & Greco, D. (2019). eUTOPIA: solUTion for Omics data PreprocessIng and Analysis. Source code for biology and medicine, 14(1), 1.

More information at: https://link.springer.com/article/10.1186/s13029-019-0071-7

#### Running the eUTOPIA Docker image (suggested)

If you don't have docker installed on your system you can install it by following the instructions at https://www.docker.com/get-docker.

The eUTOPIA docker image is available at https://hub.docker.com/r/grecolab/eutopia

#### Install OS Dependencies
Pandoc - [https://pandoc.org/installing.html](https://pandoc.org/installing.html)  
Latex Compiling - [http://www.tug.org/texlive/](http://www.tug.org/texlive/)  
Perl - [https://www.perl.org/get.html](https://www.perl.org/get.html)  
Cairo - [https://cairographics.org/download/](https://cairographics.org/download/)  

#### Install R Dependencies
```R

  #Universal Bioconductor package installation function
  install.bioc <- function(pkg){
    vers <- getRversion()
    if (vers >= "3.6"){
      if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
      BiocManager::install(pkg)
    }else{
      if (!requireNamespace("BiocInstaller", quietly = TRUE)){
        source("https://bioconductor.org/biocLite.R")
        biocLite(pkg, suppressUpdates=TRUE)
      }else{
        BiocInstaller::biocLite(pkg, suppressUpdates=TRUE)
      }
    }
  }  

  #Install impute dependency
  install.bioc("impute")

  #Install CRAN dependencies
  cran_pkgs <- c("bibtex", "RMySQL", "progress", "swamp", "infotheo", "gplots", "RColorBrewer",
  "shiny", "shinyjs", "shinyBS", "shinydashboard", "shinyFiles", "DT", "shinycssloaders", "cowplot",
  "ggplot2", "ggrepel", "WriteXLS", "rmarkdown", "VennDiagram", "grid", "futile.logger", "base2grob",
  "reshape2", "htmlTable", "devtools", "httr", "randomcoloR", "doParallel", "foreach", "import")
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
  bioc_pkgs <- c("limma", "sva", "Biobase", "biomaRt", "affy", "affyPLM", 
  "arrayQualityMetrics", "made4", "vsn", "GEOquery", "minfi",
  "IlluminaHumanMethylation450kmanifest", "IlluminaHumanMethylation450kanno.ilmn12.hg19",
  "IlluminaHumanMethylationEPICmanifest", "IlluminaHumanMethylationEPICanno.ilm10b2.hg19",
  "affyio", "simpleaffy", "yaqcaffy", "GO.db", "shinyMethyl")
  bioc_pkgs.inst <- bioc_pkgs[!(bioc_pkgs %in% rownames(installed.packages()))]
  if(length(bioc_pkgs.inst)>0){
    print(paste0("Missing ", length(bioc_pkgs.inst), " Bioconductor Packages:"))
    for(pkg in bioc_pkgs.inst){
      print(paste0("Installing Package:'", pkg, "'..."))
      install.bioc(pkg)
      print("Installed!!!")
    }
  }

  #Install latest version of GOSemSim from GitHub
  print("Installing GOSemSim from GitHub!")
  devtools::install_github("GuangchuangYu/GOSemSim")
```

#### How to run eUTOPIA from GitHub
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
#### Resolve MAX NUMBER DLLs REACHED
```R
  install.packages("usethis")
  usethis::edit_r_environ() #Opens the R enviroenment file in an editor

  #Specify the following line in the file, save and close
  R_MAX_NUM_DLLS=256

  #Restart R session and OS if required
```

### Data for demo:
>>>>>>> 432716d0b6e372d8dd9e613df7daeb05a7271fc3
Sample data can be found at: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE92900

The raw data can be downloaded from:  ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE92nnn/GSE92900/suppl/GSE92900_raw_data_files.tar.gz
** note: this file has to be decompressed

The pheno data matrix can be found at: https://github.com/Greco-Lab/eUTOPIA/blob/master/sample_data/Phenotype_File.tsv

The platform annotation file can be found at: https://github.com/Greco-Lab/eUTOPIA/blob/master/sample_data/GeneList_028005_D_GeneList_20190110.txt
