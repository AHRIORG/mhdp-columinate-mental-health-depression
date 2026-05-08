FROM rocker/r-ver:4.4.1

RUN apt-get update && apt-get install -y --no-install-recommends \
    curl \
    gfortran \
    libcurl4-openssl-dev \
    libopenblas-dev \
    libssl-dev \
    libxml2-dev \
    && rm -rf /var/lib/apt/lists/*

WORKDIR /app

COPY . /app

RUN R -q -e "install.packages('pak', repos = 'https://cloud.r-project.org')" \
    && R -q -e "pak::pkg_install(c('shiny', 'mirt'))"

ENV PORT=8080
ENV PORTABLE_IRT_HOST=0.0.0.0
ENV PORTABLE_IRT_PORT=8080

EXPOSE 8080

HEALTHCHECK --interval=30s --timeout=5s --start-period=30s --retries=3 \
  CMD curl -fsS "http://127.0.0.1:${PORT}/" >/dev/null || exit 1

CMD ["Rscript", "deploy/run-shiny.R"]
