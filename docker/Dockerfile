# --- Stage 1: Build Environment ---
FROM ubuntu:22.04 AS builder
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential cmake gfortran liblapack-dev libblas-dev ca-certificates git \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# 1. Build CHARMM
COPY c50b1.tar.gz .
RUN mkdir source && tar -xvzf c50b1.tar.gz -C source --strip-components=1
WORKDIR /build/source
RUN ./configure && \
    sed -i '45,55d' build/cmake/cmake_install.cmake && \
    cd build/cmake && \
    make -j$(nproc) && \
    make install

# 2. Clone and Install MMTSB Toolset
WORKDIR /build
RUN git clone https://github.com/mmtsb/toolset.git
WORKDIR /build/toolset/src
RUN make all

# 3. Fetch specific CPHMD toppar files (Sparse Checkout)
WORKDIR /build
RUN mkdir cphmd_tutorial && cd cphmd_tutorial && \
    git init && \
    git remote add origin https://gitlab.com/shenlab-amber-cphmd/cphmd-tutorial.git && \
    git config core.sparseCheckout true && \
    echo "hphmd_charmm/toppar/" >> .git/info/sparse-checkout && \
    git pull origin main

# --- Stage 2: Final Runtime Image ---
FROM ubuntu:22.04

RUN apt-get update && apt-get install -y \
    libgfortran5 liblapack3 libblas3 libgomp1 perl \
    && rm -rf /var/lib/apt/lists/*

# Create directories
RUN mkdir -p /opt/charmm /opt/mmtsb

# Copy CHARMM from builder
COPY --from=builder /build/source/bin /opt/charmm/bin
COPY --from=builder /build/source/lib /opt/charmm/lib
COPY --from=builder /build/source/toppar /opt/charmm/toppar

# Copy the specific CPHMD file (Keeping Original Name)
COPY --from=builder /build/cphmd_tutorial/hphmd_charmm/toppar/toppar_phmd_c22.str /opt/charmm/toppar/toppar_phmd_c22.str

# Copy MMTSB Toolset from builder
COPY --from=builder /build/toolset /opt/mmtsb

# --- Configuration & Environment Variables ---

# 1. MMTSB Specifics
ENV MMTSBDIR="/opt/mmtsb"
ENV PATH="/opt/mmtsb/perl:/opt/mmtsb/bin:/opt/charmm/bin:${PATH}"

# 2. CHARMM Specifics
ENV CHARMMEXEC="/opt/charmm/bin/charmm"
ENV CHARMMDATA="/opt/charmm/toppar"

# 3. Legacy/Standard CHARMM Variables
ENV CHM_HOME="/opt/charmm"
ENV TOPPAR="/opt/charmm/toppar"
ENV LD_LIBRARY_PATH="/opt/charmm/lib"

RUN ldconfig

WORKDIR /data
CMD ["/bin/bash"]
