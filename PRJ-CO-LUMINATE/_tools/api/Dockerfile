FROM rocker/r-ver:4.3.3

RUN apt-get update && apt-get install -y --no-install-recommends \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  gfortran \
  && rm -rf /var/lib/apt/lists/*

RUN R -q -e "install.packages(c('plumber','jsonlite','dplyr','mirt'), repos='https://cloud.r-project.org')"

WORKDIR /app
COPY . /app

ENV COLUMINATE_API_HOST=0.0.0.0
ENV COLUMINATE_API_PORT=8080
ENV COLUMINATE_SCORING_SCRIPT=/app/PRJ-CO-LUMINATE/_tools/api/scoring_runtime.R

EXPOSE 8080

WORKDIR /app/PRJ-CO-LUMINATE
CMD ["Rscript", "_tools/api/run_api.R"]
