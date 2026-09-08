FROM condaforge/miniforge3@sha256:2052d26e02808699f24f8c67f6900fce705563c04fd4bcc4e5dca434265524f4

ARG UNITFINDER_COMMIT=e80bff38a1a2718d0b6cf5fd769e6d484fc2ce28
LABEL org.opencontainers.image.title="TandemX unitFinder comparator"
LABEL org.opencontainers.image.source="https://github.com/HuangYicheng-Bio/unitFinder"
LABEL org.opencontainers.image.revision="${UNITFINDER_COMMIT}"

# Isolated comparator environment. The official repository has no dependency
# manifest, so the complete solver result is retained inside the image.
RUN mamba create -y -n unitfinder -c conda-forge -c bioconda \
      python=3.11 biopython numpy matplotlib-base scikit-learn \
      trf=4.09.1 mummer4=4.0.1 clustalo=1.2.4 time git \
    && mamba clean --all --yes

ENV PATH=/opt/conda/envs/unitfinder/bin:${PATH} \
    MPLBACKEND=Agg \
    MPLCONFIGDIR=/tmp/matplotlib

RUN mkdir -p /opt/provenance \
    && git init /opt/unitFinder \
    && git -C /opt/unitFinder fetch --depth 1 \
      https://github.com/HuangYicheng-Bio/unitFinder.git "${UNITFINDER_COMMIT}" \
    && git -C /opt/unitFinder checkout --detach FETCH_HEAD \
    && test "$(git -C /opt/unitFinder rev-parse HEAD)" = "${UNITFINDER_COMMIT}" \
    && test -z "$(git -C /opt/unitFinder status --short)" \
    && (cd /opt/unitFinder \
      && git ls-files -z | sort -z | xargs -0 sha256sum) \
      > /opt/provenance/unitfinder_tracked_files.sha256 \
    && conda list -n unitfinder --explicit > /opt/provenance/conda_explicit.txt \
    && python /opt/unitFinder/bin/unitFinder.py -h > /opt/provenance/help.txt \
    && trf 2>&1 | head -n 1 > /opt/provenance/trf_version.txt \
    && nucmer --version > /opt/provenance/nucmer_version.txt \
    && clustalo --version > /opt/provenance/clustalo_version.txt

WORKDIR /work
ENTRYPOINT ["python", "/opt/unitFinder/bin/unitFinder.py"]
CMD ["-h"]
