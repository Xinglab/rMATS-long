FROM debian:trixie

RUN apt-get update \
    && apt-get install -y --no-install-recommends \
       ca-certificates \
       curl \
       git \
    && rm -rf /var/lib/apt/lists/* \
    && mkdir /conda \
    && cd /conda \
    && curl -L 'https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh' -O \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /conda/install \
    && /conda/install/bin/conda init \
    && echo '' > /conda/install/.condarc \
    && git clone 'https://github.com/Xinglab/rMATS-long.git' /rMATS-long \
    && cd /rMATS-long \
    && ./install

ENV PATH /rMATS-long/conda_env/bin:${PATH}

# Set defaults for running the image.
# The ENTRYPOINT and CMD are empty to be compatible with
# CWL and WDL implementations that cannot override those values
WORKDIR /rMATS-long
ENTRYPOINT []
CMD []
