# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020-2021
# ALL RIGHTS RESERVED
ARG IBM_CODE_PATH_DEFAULT=/app

FROM rxn4chemistry/rxn4chemistry-base:main as builder
ARG IBM_CODE_PATH_DEFAULT
ENV IBM_CODE_PATH ${IBM_CODE_PATH_DEFAULT}

# Read build arguments
ARG GHE_TOKEN_ARG
ENV GHE_TOKEN=$GHE_TOKEN_ARG

# Install conda and PyPI dependencies
RUN conda install -y -n rxn-env tabulate=0.8.7
RUN conda run -n rxn-env pip install xxhash==2.0.0 hydra-core==1.1.0.dev4
COPY requirements_pip_nodeps.txt /
RUN conda run -n rxn-env pip install -r /requirements_pip_nodeps.txt --no-deps

RUN mkdir -p ${IBM_CODE_PATH}/bin \
    && curl -f -o ${IBM_CODE_PATH}/bin/legal_notice https://${GHE_TOKEN}@raw.github.ibm.com/rxn/rxn_utilities/develop/docker_utilities/scripts/legal_notice.sh \
    && chmod u+x ${IBM_CODE_PATH}/bin/legal_notice

# Package up environment
RUN conda install -c conda-forge conda-pack
RUN conda-pack --ignore-missing-files -n rxn-env -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar
RUN /venv/bin/conda-unpack

FROM debian:buster-slim as runtime
ARG IBM_CODE_PATH_DEFAULT
ENV IBM_CODE_PATH ${IBM_CODE_PATH_DEFAULT}

# Install RDKit runtime dependencies
RUN apt-get -y update \
    && apt-get install -yq --no-install-recommends libxrender1 libxext6 \
    && rm -rf /var/lib/apt/lists/*

# Install Python environment
COPY --from=builder /venv /venv
RUN echo "source /venv/bin/activate" > ~/.bashrc
ENV PATH /venv/bin:$PATH

RUN mkdir -p /data/output

# Install custom Python package
COPY bin/ ${IBM_CODE_PATH}/bin
RUN chmod -R u+x ${IBM_CODE_PATH}/bin
COPY setup.cfg ${IBM_CODE_PATH}/setup.cfg
COPY setup.py ${IBM_CODE_PATH}/setup.py
COPY README.md ${IBM_CODE_PATH}/README.md
COPY pyproject.toml ${IBM_CODE_PATH}/pyproject.toml
COPY LICENSE ${IBM_CODE_PATH}/LICENSE
COPY rxn_reaction_preprocessing ${IBM_CODE_PATH}/rxn_reaction_preprocessing
RUN GHE_TOKEN="" /venv/bin/python -m pip install --no-deps ${IBM_CODE_PATH}

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8

# Copy entrypoint scripts
COPY --from=builder ${IBM_CODE_PATH}/bin/legal_notice ${IBM_CODE_PATH}/bin/legal_notice
COPY bin/docker_entrypoint ${IBM_CODE_PATH}/bin/docker_entrypoint

RUN chmod -R u+x ${IBM_CODE_PATH}/bin
ENV PATH ${IBM_CODE_PATH}/bin:$PATH

WORKDIR /

ENTRYPOINT ["legal_notice", "docker_entrypoint"]

FROM runtime as dev
RUN python -m pip install pytest==6.1.2 pytest-cov==2.10.1 mypy==0.790
COPY tests /app/tests

CMD [ "bash", "-c", "cd /app && python -m pytest -sv --cov=rxn_reaction_preprocessing --cov-fail-under=50 && mypy rxn_reaction_preprocessing" ]
