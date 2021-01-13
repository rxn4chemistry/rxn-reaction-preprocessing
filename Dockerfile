# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
ARG IBM_CODE_PATH_DEFAULT=/app

FROM continuumio/miniconda3:4.8.2 as builder
ARG IBM_CODE_PATH_DEFAULT
ENV IBM_CODE_PATH ${IBM_CODE_PATH_DEFAULT}

# Read build arguments
ARG GHE_TOKEN_ARG
ENV GHE_TOKEN=$GHE_TOKEN_ARG

# Setup environment
RUN echo "source activate base" >> ~/.bashrc
ENV PATH /opt/conda/envs/base/bin:$PATH

RUN apt-get -y update \
    && apt-get -y install curl \
    && rm -rf /var/lib/apt/lists/*
# TODO: replace 'add-legal-notice' with 'develop' path after merging
RUN mkdir -p ${IBM_CODE_PATH}/bin \
    && curl -f -o ${IBM_CODE_PATH}/bin/legal_notice https://${GHE_TOKEN}@raw.github.ibm.com/rxn/rxn_utilities/add-legal-notice/docker_utilities/scripts/legal_notice.sh \
    && chmod u+x ${IBM_CODE_PATH}/bin/legal_notice

# Install conda dependencies
COPY requirements_conda.txt /
RUN while read req; do conda install -y -n base $req; done < /requirements_conda.txt

# Install custom Python package
COPY bin/ ${IBM_CODE_PATH}/bin
RUN chmod -R u+x ${IBM_CODE_PATH}/bin
COPY setup.cfg ${IBM_CODE_PATH}/setup.cfg
COPY setup.py ${IBM_CODE_PATH}/setup.py
COPY README.md ${IBM_CODE_PATH}/README.md
COPY pyproject.toml ${IBM_CODE_PATH}/pyproject.toml
COPY LICENSE ${IBM_CODE_PATH}/LICENSE
COPY rxn_reaction_preprocessing ${IBM_CODE_PATH}/rxn_reaction_preprocessing
WORKDIR ${IBM_CODE_PATH}
RUN conda run -n base pip install -e .

FROM continuumio/miniconda3:4.8.2 as internal
ENV FOR_DISTRIBUTION=false
ARG IBM_CODE_PATH_DEFAULT
ENV IBM_CODE_PATH ${IBM_CODE_PATH_DEFAULT}

COPY --from=builder ${IBM_CODE_PATH} ${IBM_CODE_PATH}
COPY --from=builder /opt/conda /opt/conda

RUN mkdir -p /data/output

ENV PATH ${IBM_CODE_PATH}/bin:${PATH}

RUN mkdir -p ${IBM_CODE_PATH}/data
COPY data/standardization-files/pistachio-200302.json ${IBM_CODE_PATH}/data/standardization-files/pistachio-200302.json

WORKDIR /

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

ENTRYPOINT ["legal_notice", "docker_entrypoint"]

# BUILD DISTRIBUTION CONTAINER
FROM internal as distribution
ENV FOR_DISTRIBUTION=true

# Add specification of not-for-distribution dependencies
COPY --from=builder /nodist_requirements_pip.txt ${IBM_CODE_PATH}/nodist_requirements_pip.txt
COPY --from=builder /nodist_requirements_conda.txt ${IBM_CODE_PATH}/nodist_requirements_conda.txt

# Remove not-for-distribution dependencies
COPY nodist_requirements_conda.txt /
RUN while read req; do conda uninstall -n base -y $req; done < /nodist_requirements_conda.txt
COPY nodist_requirements_pip.txt /
RUN conda run -n base pip uninstall -y -r /nodist_requirements_pip.txt
