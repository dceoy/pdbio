FROM dceoy/jupyter:latest

ADD . /tmp/pdbio

RUN set -e \
      && apt-get -y update \
      && apt-get -y dist-upgrade \
      && apt-get -y autoremove \
      && apt-get clean \
      && rm -rf /var/lib/apt/lists/*

RUN set -e \
      && pip install -U --no-cache-dir pip /tmp/pdbio \
      && rm -rf /tmp/pdbio

ENTRYPOINT ["jupyter"]
