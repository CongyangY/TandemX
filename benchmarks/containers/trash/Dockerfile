FROM rocker/r-ver:4.4.3@sha256:3dae5d2eeddf74f10e0a81fb6b7ae350295e288000304f438b844b2c1e00fe2c

# Comparator-only Linux environment. No host R or Python environment is changed.
ENV DEBIAN_FRONTEND=noninteractive MAKEFLAGS=-j2
RUN apt-get update && apt-get install -y --no-install-recommends \
    ca-certificates git mafft hmmer libcurl4-openssl-dev libssl-dev libxml2-dev \
    libbz2-dev liblzma-dev zlib1g-dev && rm -rf /var/lib/apt/lists/*
COPY install.R /opt/install-trash-dependencies.R
RUN mkdir -p /opt/source_archives /opt/provenance && \
    Rscript --vanilla /opt/install-trash-dependencies.R && \
    find /opt/source_archives -type f -exec sha256sum {} + > /opt/provenance/package_archives.sha256 && \
    dpkg-query -W > /opt/provenance/debian_packages.tsv

ARG TRASH_COMMIT=f7d53a013aa4664ec275f5f5544d3a3906a43c00
ARG TRASH2_COMMIT=f290a5efe2341b5525711dcaabfbfc9087edf3cb
RUN git init /opt/TRASH && git -C /opt/TRASH fetch --depth 1 https://github.com/vlothec/TRASH.git ${TRASH_COMMIT} && \
    git -C /opt/TRASH checkout --detach FETCH_HEAD && \
    git init /opt/TRASH_2 && git -C /opt/TRASH_2 fetch --depth 1 https://github.com/vlothec/TRASH_2.git ${TRASH2_COMMIT} && \
    git -C /opt/TRASH_2 checkout --detach FETCH_HEAD && \
    mkdir -p /opt/TRASH/src/mafft-linux64 && \
    ln -s /usr/bin/mafft /opt/TRASH/src/mafft-linux64/mafft.bat && \
    git -C /opt/TRASH rev-parse HEAD > /opt/provenance/trash_commit.txt && \
    git -C /opt/TRASH_2 rev-parse HEAD > /opt/provenance/trash2_commit.txt
WORKDIR /work
RUN apt-get update && apt-get install -y --no-install-recommends time && \
    rm -rf /var/lib/apt/lists/* && dpkg-query -W > /opt/provenance/debian_packages.tsv
CMD ["Rscript", "--version"]
