# ==============================================================================
# STAGE 1: Build Environment (CHARMM & MMTSB Only)
# ==============================================================================
FROM ubuntu:22.04 AS builder

ENV DEBIAN_FRONTEND=noninteractive

# 1. Install Build Dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    gfortran \
    liblapack-dev \
    libblas-dev \
    git \
    tcsh \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /build

# 2. Build CHARMM
COPY c50b1.tar.gz .
RUN mkdir source && tar -xvzf c50b1.tar.gz -C source --strip-components=1
WORKDIR /build/source
RUN ./configure && \
    sed -i '45,55d' build/cmake/cmake_install.cmake && \
    cd build/cmake && \
    make -j$(nproc) && \
    make install

# 3. Clone and Install MMTSB Toolset
WORKDIR /build
RUN git clone https://github.com/mmtsb/toolset.git
WORKDIR /build/toolset/src
RUN make all

# 4. Fetch specific CPHMD toppar files
WORKDIR /build
RUN mkdir cphmd_tutorial && cd cphmd_tutorial && \
    git init && \
    git remote add origin https://gitlab.com/shenlab-amber-cphmd/cphmd-tutorial.git && \
    git config core.sparseCheckout true && \
    echo "hphmd_charmm/toppar/" >> .git/info/sparse-checkout && \
    git pull origin main


# ==============================================================================
# STAGE 2: Final Runtime Image
# ==============================================================================
FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

# 1. Install Runtime Dependencies
RUN apt-get update && apt-get install -y \
    curl \
    git \
    tcsh \
    libgfortran5 \
    liblapack3 \
    libblas3 \
    libgomp1 \
    libxext6 \
    libx11-6 \
    libbz2-1.0 \
    zlib1g \
    openssh-client \
    && rm -rf /var/lib/apt/lists/*

# 2. Install Miniforge (Bypasses Anaconda ToS)
ENV CONDA_DIR=/opt/conda
RUN curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh" && \
    bash Miniforge3-Linux-x86_64.sh -b -p $CONDA_DIR && \
    rm Miniforge3-Linux-x86_64.sh

# Add Conda to PATH
ENV PATH=$CONDA_DIR/bin:$PATH

# 3. Install AmberTools & Python Libs
RUN conda install -y -n base -c dacase ambertools-dac=25 && \
    conda install -y -n base numpy scipy && \
    pip install MDAnalysis mdtraj && \
    conda clean -afy

# 4. Copy Artifacts from Builder
RUN mkdir -p /opt/charmm /opt/mmtsb

COPY --from=builder /build/source/bin /opt/charmm/bin
COPY --from=builder /build/source/lib /opt/charmm/lib
COPY --from=builder /build/source/toppar /opt/charmm/toppar
COPY --from=builder /build/cphmd_tutorial/hphmd_charmm/toppar/toppar_phmd_c22.str /opt/charmm/toppar/toppar_phmd_c22.str
COPY --from=builder /build/toolset /opt/mmtsb

# 5. Environment Configuration

# --- AMBERHOME ---
ENV AMBERHOME="/opt/conda"

# --- CHARMM & MMTSB ---
ENV MMTSBDIR="/opt/mmtsb"
ENV CHARMMEXEC="/opt/charmm/bin/charmm"
ENV CHARMMDATA="/opt/charmm/toppar"
ENV CHM_HOME="/opt/charmm"
ENV TOPPAR="/opt/charmm/toppar"

# --- PATH & LIBS ---
# FIX: Added /opt/charmm/bin explicitly to PATH (instead of the CHARMMEXEC file)
ENV PATH="${MMTSBDIR}/perl:${MMTSBDIR}/bin:/opt/charmm/bin:${PATH}"
ENV LD_LIBRARY_PATH="/opt/charmm/lib:/opt/conda/lib"

# 6. Finalize
RUN ldconfig

# Source amber.sh for internal python setup
RUN echo "source ${AMBERHOME}/amber.sh" >> /root/.bashrc

WORKDIR /data
CMD ["/bin/bash"]
