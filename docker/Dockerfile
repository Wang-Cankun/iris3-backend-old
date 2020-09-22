FROM rocker/rstudio:3.6.3

RUN apt-get update -qq && apt-get -y --no-install-recommends install \
  libxml2-dev \
  libcairo2-dev \
  libsqlite-dev \
  libmariadbd-dev \
  libmariadbclient-dev \
  libpq-dev \
  libssh2-1-dev \
  unixodbc-dev \
  libsasl2-dev \
  && install2.r --error \
  --deps TRUE \
  tidyverse \
  dplyr \
  devtools \
  formatR \
  remotes \
  selectr \
  caTools \
  BiocManager \
  Seurat \
  jsonlite \ 
  hdf5r \ 
  devEMF \
  gam \
  Polychrome \
  plotROC \ 
  NMF \
  igraph \
  clues \
  MLmetrics \
  AUC \ 
  ROSE \
  networkD3 \
  ggrepel \
  xml2 \
  XML \
  seqinr \ 
  grid \
  rlist \ 
  foreach \
  doParallel 


RUN installGithub.r gongx030/DrImpute immunogenomics/harmony

RUN R -e 'BiocManager::install("scater", "GSVA","BioQC"\
  "scran","slingshot","destiny","AUCell",\
  "BSgenome.Hsapiens.UCSC.hg38", "BSgenome.Mmusculus.UCSC.mm10",\
  "TxDb.Hsapiens.UCSC.hg38.knownGene","TxDb.Mmusculus.UCSC.mm10.knownGene",\
  "org.Hs.eg.db","org.Mm.eg.db"\
  )'

