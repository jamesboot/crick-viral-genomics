FROM mambaorg/micromamba:latest

# Install apt packages
USER root
RUN apt-get update \
    && apt-get install -y --no-install-recommends \
    wget \
    procps \
    libgomp1 \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install conda packages
USER $MAMBA_USER
COPY --chown=$MAMBA_USER:$MAMBA_USER env.yml /tmp/env.yml
RUN micromamba install -y -n base -f /tmp/env.yml && \
    micromamba clean --all --yes
ENV PATH /opt/conda/bin:$PATH
WORKDIR /home/mambauser

# Download SPAdes
RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz

# Extract the SPAdes tarball
RUN tar -xzf SPAdes-4.0.0-Linux.tar.gz

# Add SPAdes to PATH
ENV PATH="/home/mambauser/SPAdes-4.0.0-Linux/bin:${PATH}"

# Verify the installation by running the test
WORKDIR /home/mambauser/SPAdes-4.0.0-Linux/bin
# RUN spades.py --test
