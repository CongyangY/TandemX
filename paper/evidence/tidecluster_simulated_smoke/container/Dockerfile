FROM teambraker/braker3@sha256:5f8b3c508a9fe1bbc2e9a74dcc013eeed82f91dd5945adca7823514d9c8aecf8 AS tidehunter_builder

ARG TIDEHUNTER_COMMIT=d4a7feb2b9edfefb267a30db5d6f0685fb79ea32
ARG ABPOA_COMMIT=4dc66a346a0f98dd4858e066a1128215100631e6

USER root

RUN sed -i \
      -e 's|http://archive.ubuntu.com/ubuntu|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' \
      -e 's|http://security.ubuntu.com/ubuntu|https://mirrors.tuna.tsinghua.edu.cn/ubuntu|g' \
      /etc/apt/sources.list \
    && apt-get -o Acquire::Retries=5 update \
    && apt-get -o Acquire::Retries=5 install -y --no-install-recommends g++ make \
    && rm -rf /var/lib/apt/lists/* \
    && git clone --branch v1.4.3 --recurse-submodules \
      https://github.com/yangao07/TideHunter.git /tmp/TideHunter \
    && test "$(git -C /tmp/TideHunter rev-parse HEAD)" = "${TIDEHUNTER_COMMIT}" \
    && test "$(git -C /tmp/TideHunter ls-tree HEAD abPOA | awk '{print $3}')" = "${ABPOA_COMMIT}" \
    && cd /tmp/TideHunter \
    && make sse2=1 \
    && test "$(/tmp/TideHunter/bin/TideHunter -v | tail -1)" = "1.4.3"

FROM condaforge/miniforge3@sha256:2052d26e02808699f24f8c67f6900fce705563c04fd4bcc4e5dca434265524f4

ARG TIDECLUSTER_VERSION=1.21.2
ARG MMSEQS_VERSION_COMMIT=747c64cc8db3b4803a0f1194a3f75b3ba9f81bcb
ARG MMSEQS_SSE2_SHA256=2f18cc42a31e64cd33207c5069ef9cf5d54bfd50ec6e550b4d65ee534a4ba402
LABEL org.opencontainers.image.title="TandemX TideCluster comparator"
LABEL org.opencontainers.image.source="https://github.com/kavonrtep/TideCluster"
LABEL org.opencontainers.image.version="${TIDECLUSTER_VERSION}"

RUN mamba create -y -n tidecluster \
      -c conda-forge -c bioconda -c petrnovak \
      python=3.11 zlib=1.2.13 "tidecluster=${TIDECLUSTER_VERSION}" \
    && mamba clean --all --yes \
    && mamba run -n tidecluster TideCluster.py --version \
    && test "$(mamba run -n tidecluster TideHunter -v | tail -1)" = "1.4.3" \
    && mamba run -n tidecluster kitehor --version

ENV PATH=/opt/conda/envs/tidecluster/bin:${PATH}

COPY --from=tidehunter_builder /tmp/TideHunter/bin/TideHunter /opt/conda/envs/tidecluster/bin/TideHunter

RUN test "$(TideHunter -v | tail -1)" = "1.4.3" \
    && curl -L --retry 5 --retry-all-errors \
      https://github.com/soedinglab/MMseqs2/releases/download/16-747c6/mmseqs-linux-sse2.tar.gz \
      -o /tmp/mmseqs-linux-sse2.tar.gz \
    && echo "${MMSEQS_SSE2_SHA256}  /tmp/mmseqs-linux-sse2.tar.gz" | sha256sum -c - \
    && tar -xzf /tmp/mmseqs-linux-sse2.tar.gz -C /opt \
    && ln -sf /opt/mmseqs/bin/mmseqs /opt/conda/envs/tidecluster/bin/mmseqs \
    && test "$(mmseqs version)" = "${MMSEQS_VERSION_COMMIT}" \
    && rm /tmp/mmseqs-linux-sse2.tar.gz \
    && mamba install -y -n tidecluster -c conda-forge time \
    && mamba clean --all --yes \
    && /opt/conda/envs/tidecluster/bin/time --version

ENTRYPOINT ["TideCluster.py"]
CMD ["--help"]
