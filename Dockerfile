FROM dceoy/jupyter:latest

ADD . /tmp/pdvcf

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && pip install -U --no-cache-dir pip /tmp/pdvcf \
      && rm -rf /tmp/pdvcf

ENTRYPOINT ["jupyter"]
