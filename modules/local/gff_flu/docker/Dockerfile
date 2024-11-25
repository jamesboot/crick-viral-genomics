FROM condaforge/miniforge3

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
RUN mamba install -y python=3.10
RUN mamba install -y gcc
RUN mamba install -y miniprot=0.13
RUN mamba install -y blast
RUN mamba clean --all --yes

# Install pip packages
RUN pip install polars==0.19.18
RUN pip install polars-lts-cpu==0.19.18
RUN pip install biopython==1.80
RUN pip install gfflu

ENV PATH=/opt/conda/bin:$PATH
WORKDIR /home

CMD ["bash"]
