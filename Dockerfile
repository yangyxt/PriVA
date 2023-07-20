# Use Ubuntu 22.04 as the base image
FROM ubuntu:22.04

# Metadata
LABEL maintainer="Yangyxt yangyxt@gmail.com"

# Avoid warnings by switching to noninteractive
ENV DEBIAN_FRONTEND=noninteractive

# Update the apt package list
RUN apt-get update -y \
    && apt-get upgrade -y \
    && apt-get install -y wget bzip2

# Install Miniconda3
RUN wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh \
    && bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda \
    && rm Miniconda3-latest-Linux-x86_64.sh

# Put Miniconda on PATH
ENV PATH="$HOME/miniconda/bin:${PATH}"

# Run bash init
RUN conda init bash && source $HOME/.bashrc

# Create a new conda environment from a yml file
COPY acmg_conda.yml .
RUN conda env create -f acmg_conda.yml

# Clean up
RUN apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# Switch back to dialog for any ad-hoc use of apt-get
ENV DEBIAN_FRONTEND=dialog
