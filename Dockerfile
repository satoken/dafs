# From ubuntu:18.04
From debian:11

ENV DEBIAN_FRONTEND=noninteractive

WORKDIR /workspaces

# use the official package
# ADD https://www.tbi.univie.ac.at/RNA/download/debian/debian_10/viennarna_2.4.15-1_amd64.deb .
# ADD https://www.tbi.univie.ac.at/RNA/download/debian/debian_10/viennarna-dev_2.4.15-1_amd64.deb .
# ADD https://www.tbi.univie.ac.at/RNA/download/debian/debian_10/python3-rna_2.4.15-1_amd64.deb .

RUN apt-get update \
    && apt-get -y install build-essential wget cmake \
            libglpk-dev libgsl-dev libgmp-dev pkg-config \
    && apt-get autoremove -y \
    && apt-get clean -y \
    && rm -rf /var/lib/apt/lists/*

# build from the source
ARG VIENNA_VER=2.4.18
RUN wget -q https://www.tbi.univie.ac.at/RNA/download/sourcecode/2_4_x/ViennaRNA-${VIENNA_VER}.tar.gz \
    && tar zxvf ViennaRNA-${VIENNA_VER}.tar.gz \
    && cd ViennaRNA-${VIENNA_VER} \
    && ./configure --without-perl --without-python --without-python3 --without-forester --without-rnalocmin \
    && make && make install \
    && cd .. && rm -rf ViennaRNA-${VIENNA_VER} ViennaRNA-${VIENNA_VER}.tar.gz 

COPY . .

RUN rm -rf build && mkdir build \
    && cd build \
    && cmake -DCMAKE_BUILD_TYPE=Release .. \
    #&& cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_EXE_LINKER_FLAGS='-static' -DCMAKE_FIND_LIBRARY_SUFFIXES='.a' -DBUILD_SHARED_LIBRARIES=OFF .. \ # static link
    && make && make install
