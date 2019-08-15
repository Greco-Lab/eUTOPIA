FROM rocker/shiny:3.6.0

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  IlluminaHumanMethylation450kmanifest

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  IlluminaHumanMethylation450kanno.ilmn12.hg19

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

# Install system dependencies
RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libv8-3.14-dev \
  libssl-dev \
  libxml2-dev \
  libmariadb-client-lgpl-dev \
  libgit2-dev \
  libssh2-1-dev \
  libmagick++-dev \
  gtk-doc-tools \
  libgeos-dev \
  libgdal-dev \
  libudunits2-dev \
  texlive-full

# Edit install2.r add INSTALL_opts in install.packages() to avoid loading tests
RUN sed -i '10itrace(install.packages, quote({INSTALL_opts<-"--no-test-load";print("INSTALL_opts:");print(INSTALL_opts)}))' /usr/local/bin/install2.r

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  impute \ 
  units \ 
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

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

# Install GitHub packages
RUN installGithub.r jrowen/rhandsontable \
  hms-dbmi/UpSetR

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
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
  minfi

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  IlluminaHumanMethylationEPICmanifest

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  IlluminaHumanMethylationEPICanno.ilm10b2.hg19 

# Clear temporary downloaded R packages
RUN rm -rf /tmp/downloaded_packages/

# Install Bioconductor and CRAN packages
RUN install2.r --error \
  -r "https://cran.rstudio.com" \
  -r "http://www.bioconductor.org/packages/release/bioc" \
  -r "http://www.bioconductor.org/packages/data/annotation" \
  --deps TRUE \
  affyio \ 
  simpleaffy \ 
  yaqcaffy \ 
  GO.db \ 
  shinyMethyl

# Install GitHub packages
RUN installGithub.r GuangchuangYu/GOSemSim \
  && rm -rf /tmp/downloaded_packages/

# Create .Renviron with variables
#RUN echo -e 'R_LIBS_USER="/home/shiny/Rlibs"\nR_MAX_NUM_DLLS=256' > /root/.Renviron \
RUN echo 'R_LIBS_USER="/home/shiny/Rlibs"\nR_MAX_NUM_DLLS=256' > /root/.Renviron \
  && cp /root/.Renviron /home/shiny/ \
  && chown shiny.shiny /home/shiny/.Renviron \
  && mkdir /home/shiny/Rlibs \
  && chown -R shiny.shiny /home/shiny/Rlibs \
  && chmod -R 777 /home/shiny/Rlibs

# Cooy your Shiny apps to the docker image
COPY ./mountpoints/apps/eUTOPIA /srv/shiny-server/eUTOPIA
COPY ./mountpoints/apps/example-app /srv/shiny-server/example-app
COPY ./mountpoints/apps/shinyMethyl_QC /srv/shiny-server/shinyMethyl_QC

# Recursively change bit mode of shiny apps
#RUN chmod -R 777 /srv/shiny-server \

