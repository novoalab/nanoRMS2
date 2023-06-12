# Set the base image to guppy image
FROM genomicpariscentre/guppy-gpu:3.6.1

# File Author / Maintainer
MAINTAINER Leszek Pryszcz <lpryszcz@crg.es>

ARG PACKAGE_VERSION=1.0a
ARG BUILD_PACKAGES="wget apt-transport-https"
ARG DEBIAN_FRONTEND=interactive

ARG PYTHON_VERSION=3.7.12
ARG MEME_VERSION=5.3.0

ENV TZ=Europe/Madrid
RUN ln -snf /usr/share/zoneinfo/$TZ /etc/localtime && echo $TZ > /etc/timezone

# system-dependencies
RUN apt update || :
RUN apt install -y \
    build-essential \
    ghostscript \
    libbz2-dev \
    liblzma-dev \
    libhdf5-dev \
    libexpat1-dev \
#    libxml2 libxslt1.1 \
    make \
    mpich \
    wget \
    zlib1g-dev \
# python3: tk-dev checkinstall 
    libreadline-gplv2-dev libncursesw5-dev libssl-dev libsqlite3-dev libgdbm-dev libc6-dev libbz2-dev libffi-dev \
    && rm -rf /var/lib/apt/lists/*

# python
RUN cd /usr/src && \
    wget -q https://www.python.org/ftp/python/$PYTHON_VERSION/Python-$PYTHON_VERSION.tgz \
    && tar xzf Python-$PYTHON_VERSION.tgz \
    && cd Python-$PYTHON_VERSION \
    && ./configure --enable-optimizations \
    && make install \
    && rm -r /usr/src/Python-$PYTHON_VERSION* \    
    && cd /usr/local/bin \
    && ln -s python3 python \
    && ln -s pip3 pip

# meme
RUN wget -q https://meme-suite.org/meme/meme-software/$MEME_VERSION/meme-$MEME_VERSION.tar.gz \
    && tar xpfz meme-$MEME_VERSION.tar.gz \
    && cd meme-$MEME_VERSION \
    && ./configure --prefix=/opt/meme --enable-build-libxml2 --enable-build-libxslt \
    && make install 
#    && rm -r meme-$MEME_VERSION*
ENV PATH /opt/meme/bin:/opt/meme/libexec/meme:$PATH

# app
#RUN wget -q https://github.com/novoalab/$tool/archive/refs/heads/main.tar.gz && \
#    tar xpfz main.tar.gz && rm main.tar.gz && mv $tool-main /opt/app
RUN wget -q https://public-docs.crg.es/enovoa/public/lpryszcz/src/nanoRMS2.tar.gz \
    && tar xpfz nanoRMS2.tar.gz && rm nanoRMS2.tar.gz \
    && mv nanoRMS2 /opt/app
WORKDIR /opt/app
# tomto requires pip<20 and numpy<1.20 installed before
RUN pip install "pip<20" "numpy<1.20" && pip install -r requirements.txt

# link guppy models
RUN for ext in cfg jsn; do for f in /opt/app/data/*.$ext; do echo $f; ln -s $f /opt/ont/guppy/data/; done; done

## need to comment out os.setpgrp() lines in /run
RUN sed -i '/os.setpgrp/d' src/predict_mods.py