FROM debian:bullseye

# Install dependencies and TeX Live
RUN apt-get update && apt-get install -y \
    texlive-full \
    latexmk \
    curl \
    git \
    wget \
    build-essential \
    && apt-get clean && rm -rf /var/lib/apt/lists/*
