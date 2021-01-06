FROM continuumio/miniconda

WORKDIR /app


# Create and activate conda environment
ARG conda_env=rxn_reaction_preprocessing
RUN conda create -c rdkit -n ${conda_env} rdkit=2020.09.1.0 python=3.6
ENV PATH /opt/conda/envs/${conda_env}/bin:$PATH
RUN /bin/bash -c "source activate ${conda_env}"

# Copy the files
COPY . ./

# Install the dependencies
ARG GHE_TOKEN_ARG=''
RUN GHE_TOKEN=$GHE_TOKEN_ARG pip install -e .

# Run the tests
RUN pip install pytest
RUN python -m pytest

WORKDIR /