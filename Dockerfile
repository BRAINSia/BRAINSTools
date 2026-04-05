# ===========================================================================
# BRAINSTools Docker Build Environment
# Addresses https://github.com/BRAINSia/BRAINSTools/issues/382
#
# Usage:
#   docker build -t brainstools .
#   docker run -it brainstools bash
#
# Multi-stage build:
#   Stage 1 (builder) - installs system deps, builds SuperBuild EPs and
#                        BRAINSTools in a single layer.
#   Stage 2 (runtime) - copies only the install tree into a minimal image.
# ===========================================================================

# ---------------------------------------------------------------------------
# Stage 1: Build
# ---------------------------------------------------------------------------
FROM ubuntu:24.04 AS builder

ARG DEBIAN_FRONTEND=noninteractive
ARG BUILD_TYPE=Release
ARG NPROC=4

# System packages required for the SuperBuild (ITK, VTK, ANTs, TBB, etc.)
RUN apt-get update && apt-get install -y --no-install-recommends \
      build-essential \
      cmake \
      ninja-build \
      git \
      python3 \
      python3-dev \
      # VTK / rendering dependencies
      libgl1-mesa-dev \
      libglu1-mesa-dev \
      libxt-dev \
      libxrender-dev \
      libxext-dev \
      # FFTW (used by several BRAINSTools modules)
      libfftw3-dev \
      # BLAS / LAPACK
      libblas-dev \
      liblapack-dev \
      # TBB (optional system package; SuperBuild can also build its own)
      libtbb-dev \
      # Miscellaneous build requirements
      pkg-config \
      zlib1g-dev \
      libssl-dev \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Verify cmake version meets minimum requirement (3.20.6)
RUN cmake --version

WORKDIR /src/BRAINSTools
COPY . .

# Phase I  – SuperBuild: build external projects (ITK, VTK, ANTs, …)
# Phase II – BRAINSTools proper is built automatically by the SuperBuild
RUN mkdir -p /build && \
    cmake -S /src/BRAINSTools -B /build \
      -G Ninja \
      -DCMAKE_BUILD_TYPE:STRING=${BUILD_TYPE} \
      -DBRAINSTools_SUPERBUILD:BOOL=ON \
      -DBRAINSTools_USE_QT:BOOL=OFF \
      -DBRAINSTools_BUILD_TESTING:BOOL=ON \
      -DUSE_ANTS:BOOL=ON \
      -DCMAKE_INSTALL_PREFIX:PATH=/opt/brainstools \
    && cmake --build /build -j${NPROC} \
    && cmake --install /build

# ---------------------------------------------------------------------------
# Stage 2: Runtime – carry only the install tree
# ---------------------------------------------------------------------------
FROM ubuntu:24.04 AS runtime

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
      libgl1-mesa-glx \
      libglu1-mesa \
      libxt6 \
      libfftw3-3 \
      libblas3 \
      liblapack3 \
      libtbb12 \
      python3 \
      ca-certificates \
    && rm -rf /var/lib/apt/lists/*

COPY --from=builder /opt/brainstools /opt/brainstools

ENV PATH="/opt/brainstools/bin:${PATH}" \
    LD_LIBRARY_PATH="/opt/brainstools/lib:${LD_LIBRARY_PATH}"

# Smoke-test: make sure a core executable is present
RUN ls /opt/brainstools/bin/BRAINSFit || true

WORKDIR /data
CMD ["bash"]
