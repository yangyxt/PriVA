#!/bin/bash

# Build the Docker image
docker build -t my-image .

# Run a container from the image, invoking the bash shell
docker run --name my-container -it my-image /bin/bash << EOF
  # Check Miniconda installation
  conda --version
  if [ $? -eq 0 ]; then
    echo "Miniconda is installed."
  else
    echo "Miniconda is not installed."
    exit 1
  fi

  # Check the conda environment
  conda env list | grep my-environment
  if [ $? -eq 0 ]; then
    echo "Conda environment is created."
  else
    echo "Conda environment is not created."
    exit 1
  fi
EOF
