# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Metadata
LABEL maintainer="Yangyxt yangyxt@gmail.com"

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Update the apt package list
RUN apt update -y && \
    apt upgrade -y && \
    apt install -y wget bzip2

# Install Miniconda3 and add to PATH
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p /root/miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

# Put Miniconda on PATH
ENV PATH="/root/miniconda/bin:${PATH}"

# Create a new conda environment from a yml file
COPY acmg_conda.yml .
RUN conda init bash && . $HOME/.bashrc && \
    conda env create -f acmg_conda.yml

# Clean up
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog
# Activate the environment by default
ENV CONDA_DEFAULT_ENV=acmg
ENV PATH="/root/miniconda/envs/${CONDA_DEFAULT_ENV}/bin:${PATH}"

COPY ./common_bash_utils.sh /scripts/common_bash_utils.sh
COPY ./annotation_per_family.sh /scripts/annotation_per_family.sh
COPY ./data /data

ENTRYPOINT [ "/scripts/annotation_per_family.sh" ]
CMD []
