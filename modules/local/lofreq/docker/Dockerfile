FROM mambaorg/micromamba:latest

# Install apt packages
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    wget \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV PATH /opt/conda/bin:$PATH
WORKDIR /home/mambauser

# Base stage: Common setup for both test and prod
FROM condaforge/miniforge3 AS base

# Install apt packages
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    procps \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Set up conda channels with priority
RUN conda config --add channels defaults \
    && conda config --add channels bioconda \
    && conda config --add channels conda-forge \
    && conda config --set channel_priority strict

# Install conda packages
RUN mamba install -y python=3.11
RUN mamba install -y lofreq=2.1.5
RUN mamba install -y pyvcf
RUN mamba clean --all --yes

# Install pip packages
RUN pip install numpy
RUN pip install matplotlib
RUN pip install scipy

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /home

CMD ["bash"]


