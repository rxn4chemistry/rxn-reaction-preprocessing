FROM continuumio/miniconda

ENV RUNNING_IN_DOCKER 1
ARG conda_env=data-preprocessor
RUN mkdir -p /data/output

WORKDIR /app


COPY . ./

RUN conda env create -f environment.yml
ENV PATH /opt/conda/envs/${conda_env}/bin:$PATH
RUN /bin/bash -c "source activate ${conda_env}"

ENTRYPOINT ["python", "./scripts/preprocess.py"]