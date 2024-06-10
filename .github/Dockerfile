# Base image
FROM inseefrlab/onyxia-rstudio:r4.4.0

ENV RENV_VERSION 1.0.7
RUN R -e "install.packages('remotes', repos = c(CRAN = 'https://cloud.r-project.org'))"
RUN R -e "remotes::install_github('rstudio/renv@${RENV_VERSION}')"

WORKDIR ${HOME}/DT-way-to-net-zero
COPY renv.lock renv.lock

ENV RENV_PATHS_LIBRARY renv/library

RUN R -e "renv::restore()"