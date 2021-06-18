FROM dceoy/jupyter:latest AS builder

COPY --from=dceoy/samtools:latest /usr/local/src/samtools /usr/local/src/samtools
COPY --from=dceoy/bcftools:latest /usr/local/src/bcftools /usr/local/src/bcftools

ADD . /tmp/pdbio

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        libbz2-dev libcurl4-gnutls-dev libncurses5-dev libgsl-dev libperl-dev \
        liblzma-dev libssl-dev libz-dev make \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && cd /usr/local/src/samtools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/samtools \
      && make clean \
      && ./configure \
      && make \
      && make install \
      && cd /usr/local/src/bcftools/htslib-* \
      && make clean \
      && ./configure \
      && make \
      && cd /usr/local/src/bcftools \
      && make clean \
      && ./configure --enable-libgsl --enable-perl-filters \
      && make \
      && make install

RUN set -e \
      && pip install -U --no-cache-dir pip /tmp/pdbio \
      && rm -rf /tmp/pdbio


FROM dceoy/jupyter:latest

COPY --from=builder /usr/local /usr/local

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y install --no-install-recommends --no-install-suggests \
        libcurl4-gnutls-dev libgsl-dev libncurses5 \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

ENTRYPOINT ["jupyter"]
