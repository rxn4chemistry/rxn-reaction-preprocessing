# LICENSED INTERNAL CODE. PROPERTY OF IBM.
# IBM Research Zurich Licensed Internal Code
# (C) Copyright IBM Corp. 2020
# ALL RIGHTS RESERVED
FROM continuumio/miniconda3 as builder

RUN mkdir -p /app

ARG GHE_TOKEN_ARG
ENV GHE_TOKEN=$GHE_TOKEN_ARG

RUN apt-get -y update \
 && apt-get -y install curl \
 && rm -rf /var/lib/apt/lists/*

# TODO: replace 'add-legal-notice' with 'develop' path after merging
RUN curl -f -o /app/legal_notice.sh https://${GHE_TOKEN}@raw.github.ibm.com/rxn/rxn_utilities/add-legal-notice/docker_utilities/scripts/legal_notice.sh \
 && chmod u+x /app/legal_notice.sh

FROM continuumio/miniconda3 as main

COPY --from=builder /app/legal_notice.sh /app/legal_notice.sh

ENV PATH /opt/conda/envs/base/bin:$PATH
RUN /bin/bash -c "source activate base"

RUN mkdir -p /data/output

COPY setup.cfg /app/setup.cfg
COPY setup.py /app/setup.py
COPY README.md /app/README.md
COPY pyproject.toml /app/pyproject.toml
COPY LICENSE /app/LICENSE
COPY bin /app/bin
COPY rxn_reaction_preprocessing /app/rxn_reaction_preprocessing
RUN mkdir -p /app/data/
COPY data/standardization-files/pistachio-200302.json /app/data/standardization-files/pistachio-200302.json

WORKDIR /app
RUN conda run -n base pip install -e .

COPY non_distributed_requirements.txt /app/
#RUN conda run -n base pip uninstall -y -r /app/non_distributed_requirements.txt

COPY docker_entrypoint.sh /app/docker_entrypoint.sh
RUN chmod u+x /app/docker_entrypoint.sh

COPY non_distributed_requirements.txt /app/non_distributed_requirements.txt
COPY non_distributed_requirements.sh /app/non_distributed_requirements.sh
COPY download_non_distributed_requirements.sh /app/download_non_distributed_requirements.sh

WORKDIR /

ENV LC_ALL C.UTF-8
ENV LANG C.UTF-8
ENV IBM_CODE_PATH /app
ENTRYPOINT ["/app/legal_notice.sh", "/app/docker_entrypoint.sh"]
