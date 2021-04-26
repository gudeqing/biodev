FROM dceoy/clir:latest
# DEBIAN_FRONTEND这个环境变量，告知操作系统应该从哪儿获得用户输入。如果设置为"noninteractive"，你就可以直接运行命令，而无需向用户请求输入（所有操作都是非交互式的）
ENV DEBIAN_FRONTEND noninteractive

ADD https://raw.githubusercontent.com/dceoy/print-github-tags/master/print-github-tags /usr/local/bin/print-github-tags

# set -e 表示在接下来执行的命令中，如果命令的返回值不为0，那么会使所在的进程或shell退出
RUN set -e \
      && ln -sf bash /bin/sh \
      && ln -s python3 /usr/bin/python \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        ca-certificates curl g++ gcc libjpeg-turbo8-dev libpng-dev libtbb-dev \
        libz-dev make perl perl-doc python3 \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/* \
      && clir update

# 安装bowtie
RUN set -e \
      && cd /usr/local/bin \
      && curl -SL -O https://github.com/BenLangmead/bowtie/releases/download/v1.3.0/bowtie-1.3.0-linux-x86_64.zip \
      && unzip bowtie-1.3.0-linux-x86_64.zip \
      && rm bowtie-1.3.0-linux-x86_64.zip

# 安装bowtie2
RUN set -e \
      && cd /usr/local/bin \
      && curl -SL -O https://github.com/BenLangmead/bowtie2/releases/download/v2.4.2/bowtie2-2.4.2-linux-x86_64.zip \
      && unzip bowtie2-2.4.2-linux-x86_64.zip \
      && rm bowtie2-2.4.2-linux-x86_64.zip

ENV PATH=/usr/local/bin/bowtie2-2.4.2-linux-x86_64/:${PATH}

# 安装star
RUN set -e \
      && cd /usr/local/bin \
      && curl -SL -O https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz \
      && tar -xzf 2.7.8a.tar.gz

ENV PATH=/usr/local/bin/STAR-2.7.8a/bin/Linux_x86_64_static/:${PATH}

# 安装HISAT2: http://daehwankimlab.github.io/hisat2/download/
RUN set -e \
      && curl -SL -O https://cloud.biohpc.swmed.edu/index.php/s/oTtGWbWjaxsQ2Ho/download \
      && unzip download && rm download \
      && mv hisat2-2.2.1 /usr/local/bin

ENV PATH=/usr/local/bin/hisat2-2.2.1/:${PATH}

# 安装rsem
RUN set -eo pipefail \
      && chmod +x /usr/local/bin/print-github-tags \
      && print-github-tags --release --latest --tar deweylab/RSEM \
        | xargs -i curl -SL {} -o /tmp/rsem.tar.gz \
      && tar xvf /tmp/rsem.tar.gz -C /usr/local/src --remove-files \
      && mv /usr/local/src/RSEM-* /usr/local/src/RSEM \
      && cd /usr/local/src/RSEM \
      && make \
      && make ebseq \
      && make pRSEM \
      && make install \
      && cd /usr/local/bin/ \
      && ln -s samtools-1.3/samtools .

ENTRYPOINT ["/usr/local/bin/rsem-calculate-expression"]
CMD ["echo", "usage: docker run -it <image_id> rsem-calculate-expression [options] --paired-end r1.fq,r1.fq r2.fq,r2.fq ref_name sample_name"]
