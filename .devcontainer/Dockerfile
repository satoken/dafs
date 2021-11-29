# See here for image contents: https://github.com/microsoft/vscode-dev-containers/tree/v0.195.0/containers/cpp/.devcontainer/base.Dockerfile
# [Choice] Debian / Ubuntu version (use Debian 11/9, Ubuntu 18.04/21.04 on local arm64/Apple Silicon): debian-11, debian-10, debian-9, ubuntu-21.04, ubuntu-20.04, ubuntu-18.04
ARG VARIANT=debian-11
FROM mcr.microsoft.com/vscode/devcontainers/cpp:0-${VARIANT}

# [Optional] Uncomment this section to install additional packages.
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends pkg-config libglpk-dev 

# build from the source
ARG VIENNA_VER=2.4.18
RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-${VIENNA_VER}.tar.gz \
    && tar zxvf ViennaRNA-${VIENNA_VER}.tar.gz \
    && cd ViennaRNA-${VIENNA_VER} \
    && ./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin \
    && make && make install \
    && cd .. && rm -rf ViennaRNA-${VIENNA_VER} ViennaRNA-${VIENNA_VER}.tar.gz 
