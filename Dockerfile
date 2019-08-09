FROM rocker/shiny:3.6.0

# Install system dependencies
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libv8-3.14-dev \
  libssl-dev \
  libxml2-dev \
  libmariadb-client-lgpl-dev \
  libgit2-dev \
  libssh2-1-dev \
  libmagick++-dev \
  texlive-full

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  --deps TRUE \
  impute \ 
  bibtex \ 
  RMySQL \ 
  progress \ 
  swamp \ 
  infotheo \ 
  gplots \ 
  RColorBrewer \ 
  shiny \ 
  shinyjs \ 
  shinyBS \ 
  shinydashboard \ 
  shinyFiles \ 
  DT \ 
  shinycssloaders \ 
  cowplot \ 
  ggplot2 \ 
  ggrepel \ 
  WriteXLS \ 
  rmarkdown \ 
  VennDiagram \ 
  grid \ 
  futile.logger \ 
  base2grob \ 
  reshape2 \ 
  htmlTable \ 
  devtools \ 
  httr \ 
  randomcoloR \ 
  doParallel \ 
  foreach \ 
  import

# Install GitHub packages
RUN installGithub.r jrowen/rhandsontable \
  hms-dbmi/UpSetR

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  --deps TRUE \
  limma \ 
  sva \ 
  Biobase \ 
  biomaRt \ 
  affy \ 
  affyPLM \ 
  arrayQualityMetrics \ 
  made4 \ 
  vsn \ 
  GEOquery \ 
  minfi \ 
  IlluminaHumanMethylation450kmanifest \ 
  IlluminaHumanMethylation450kanno.ilmn12.hg19 \ 
  IlluminaHumanMethylationEPICmanifest \ 
  IlluminaHumanMethylationEPICanno.ilm10b2.hg19 \ 
  affyio \ 
  simpleaffy \ 
  yaqcaffy \ 
  GO.db \ 
  shinyMethyl

# Install GitHub packages
RUN installGithub.r GuangchuangYu/GOSemSim \
  && rm -rf /tmp/downloaded_packages/

# Create .Renviron with variables
RUN echo -e 'R_LIBS_USER="/home/shiny/Rlibs"\nR_MAX_NUM_DLLS=256' > /root/.Renviron \
  && cp /root/.Renviron /home/shiny/ \
  && chown shiny.shiny /home/shiny/.Renviron \
  && mkdir /home/shiny/Rlibs \
  && chown -R shiny.shiny /home/shiny/Rlibs \
  && chmod -R 777 /home/shiny/Rlibs

# Recursively change bit mode of shiny apps
#RUN chmod -R 777 /srv/shiny-server \

