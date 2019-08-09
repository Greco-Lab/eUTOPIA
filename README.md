Docker for eUTOPIA
=======================

### A solUTion for Omics data PreprocessIng and Analysis

Graphically accessible guided workflow for preprocessing and analysis of omics data. Supports Agilent 2-color, Agilent 1-color, Affymetrix, and Illumina methylation microarray platforms (Ongoing efforts to add support for RNA-Seq data). Discreetly separated steps in analysis designed in R Shiny, incorporates widely used microarray analysis practices and R packages. Reporting is and data interpretation is leveraged from dynamically generated plots.

Reference Paper:
> Marwah, V. S., Scala, G., Kinaret, P. A. S., Serra, A., Alenius, H., Fortino, V., & Greco, D. (2019). eUTOPIA: solUTion for Omics data PreprocessIng and Analysis. Source code for biology and medicine, 14(1), 1.

#### Install Docker application for Windows 10 Home

https://docs.docker.com/toolbox/toolbox_install_windows/

#### Install Docker application for Windows 10 64bit: Pro, Enterprise or Education

https://docs.docker.com/docker-for-windows/install/

#### Install Docker application for Mac

https://docs.docker.com/docker-for-mac/install/

#### Pull the eUTOPIA docker image:

docker pull grecolab/eutopia

#### Run the eUTOPIA docker image, map http port 3838 on the host port 8787:

docker run --rm -p 8787:3838 grecolab/eutopia

#### Using eUTOPIA:

Open your browser and visit (if you choose to map to other port replace 8787 in the url with the correct one): http://localhost:8787/eUTOPIA

