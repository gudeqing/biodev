ARG BASE_IMAGE=docker-reg.basebit.me:5000/base/enigma2-vnc-app:1.6.1
FROM ${BASE_IMAGE}

ARG VERSION
LABEL software.version="${VERSION}" \
      software.vendor="https://www.basebit.ai/"

RUN apt-get update && apt-get install -y software-properties-common gnupg2 && apt-key adv --keyserver keyserver.ubuntu.com --recv-keys E298A3A825C0D65DFD57CBB651716619E084DAB9 && add-apt-repository 'deb https://cloud.r-project.org/bin/linux/ubuntu bionic-cran35/'

RUN apt-get update \
        && DEBIAN_FRONTEND=noninteractive apt-get install --no-install-recommends -y \
        wget \
        r-base \
        r-base-dev \
        r-base-core \
        r-recommended \
        libnss3 \
        libasound2 \
        make \
        cmake \
        build-essential \
        autoconf \
        automake \
        libtool \
        flex \
        bison \
        libxml2-dev \
        libclang-dev \
        libedit2 \
        libssl-dev \
        lsb-release \
        multiarch-support \
        psmisc \
        procps \
        gfortran \
        vim \
        pandoc \
        pandoc-citeproc \
        libcurl4-gnutls-dev \
        libcairo2-dev \
        libxt-dev \
        libssh2-1-dev \
        libssl1.0.0

RUN wget -P /tmp https://download1.rstudio.org/desktop/bionic/amd64/rstudio-1.3.1073-amd64.deb \
        && apt install --no-install-recommends -y /tmp/rstudio-1.3.1073-amd64.deb \
        && apt-get -y autoremove \
        && apt-get clean && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/* \
        && mkdir -p /enigma

# 安装miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-py37_4.9.2-Linux-x86_64.sh \
    && sh Miniconda3-py37_4.9.2-Linux-x86_64.sh -b \
    # 安装rust，因为有些R包如gifski需要
    && wget https://static.rust-lang.org/rustup/rustup-init.sh \
    && sh rustup-init.sh -y

ENV PATH /root/miniconda3/condabin:$PATH
# 安装R和核心包
RUN conda create -n R3.6.3 \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda install r-essentials 

ENV PATH  /root/miniconda3/envs/R3.6.3/bin:$PATH

# 设置环境变量
ENV RSTUDIO_APP_DIR=/headless/rstudio \
    RSTUDIO_STARTUP_DIR=/headless/rstudio/.rstudio \
    RSTUDIO_WORK_DIR=/headless/rstudio/.rstudio \
    RSTUDIO_CRAN_URL=https://cloud.r-project.org \
    RSTUDIO_CRAN_CERT= \
    VNC_APP_NAME="XDP Sandbox - RStudio"

# 复制文件
COPY ./root /
COPY ./scripts/ $RSTUDIO_APP_DIR/scripts/

# 使用conda安装其他R包
RUN rm /bin/sh && ln -s /bin/bash /bin/sh
RUN source /root/miniconda3/bin/activate R3.6.3 && conda install r-devtools r-openxlsx r-ggthemes r-ggiraph \
    r-ggstance r-GGally r-ggalt r-ggforce r-ggrepel r-ggraph r-ggExtra \
    r-gganimate r-plotROC r-ggnetwork r-car r-Hmisc r-multcomp r-pbkrtest \
    r-mvtnorm r-SparseM r-lme4 r-gbm r-ranger r-ROCR r-multcomp r-hflights \
    r-USAboundaries r-fueleconomy r-nasaweather r-gapminder r-janeaustenr \
    r-RColorBrewer r-gplots r-survival r-ggpubr r-survminer r-gifski \
    r-ggstatsplot r-ClassDiscovery r-glmnet r-pheatmap r-ggfortify r-lars r-randomForest r-wgcna \
    # 使用R安装其他R包
    && Rscript -e "install.packages(c('MatchIt', 'ggplot2movies', 'gunsales', 'geomnet', 'vtreat', 'stationaRy'), repos = 'http://cran.us.r-project.org')" \
    && Rscript -e "devtools::install_github(c('ricardo-bion/ggradar', 'Ather-Energy/ggTimeSeries', 'hafen/housingData','hadley/babynames', 'hadley/neiss', 'hadley/yrbss','hadley/usdanutrients'), dependencies=TRUE)"

# 使用conda 安装常用生信分析包
RUN source /root/miniconda3/bin/activate R3.6.3 && conda install -c bioconda \
    bioconductor-deseq2 bioconductor-limma bioconductor-edgeR bioconductor-airway bioconductor-maftools \
    bioconductor-org.hs.eg.db bioconductor-genefilter bioconductor-clusterprofiler \
    bioconductor-GSEABase bioconductor-GSVA

# 激活R.3.6.3环境
RUN echo "source /root/miniconda3/bin/activate R3.6.3" >> ~/.bashrc

CMD ["bash", "-c", "/headless/rstudio/scripts/cmd.sh"]
