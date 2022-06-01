
## tag: cryanking/temperature_container sha256:39594b0a044d97bdf8114a7b4bf5567e8b2790808a7284fc7234acfd985e2ffb
FROM rocker/tidyverse:4.1.2@sha256:447776015b8ef273fdd96fb1ae804e2e048cf7960f1590486f6dd7a590f972f8

ARG R_VERSION
ARG BUILD_DATE
ARG CRAN
ENV BUILD_DATE ${BUILD_DATE:-2022-05-01}
ENV R_VERSION=${R_VERSION:-4.1.2} \
    CRAN=${CRAN:-https://cran.microsoft.com/snapshot/2022-05-01}} \ 
    TERM=xterm \
    MRAN_BUILD_DATE=2022-05-01


ENV LANG en_US.UTF-8  
ENV LANGUAGE en_US:en  
ENV LC_ALL en_US.UTF-8   


RUN useradd docker \
	&& mkdir /home/docker \
	&& chown docker:docker /home/docker \
	&& addgroup docker staff

RUN  apt-get update \
	&& DEBIAN_FRONTEND="noninteractive" apt-get install -y --no-install-recommends cmake 
	
  RUN echo "options(repos = c(CRAN = 'https://cran.microsoft.com/snapshot/${MRAN_BUILD_DATE}'), download.file.method = 'libcurl')"  >> ./usr/local/lib/R/etc/Rprofile.site

RUN install2.r --error \
    lme4 \ 
    tableone \
    pbkrtest \
    kableExtra
    
    
RUN chmod -R 755 /root


