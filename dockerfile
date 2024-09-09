########## Ubuntu-based

FROM ubuntu:latest

########## Install utilized tools

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y python3 python3-pip tree ccache git make vim nano && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

########## Install prerequisites

WORKDIR /usr/local/bin

RUN apt-get update && \
    DEBIAN_FRONTEND=noninteractive && \
    apt-get install -y cutadapt ncbi-blast+ clustalo iqtree && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN git clone https://github.com/rcedgar/usearch12.git && \
    cd usearch12/src && \
    make && \
    cd ../.. && \
    mv ./usearch12/bin/usearch12 ./usearch12/bin/usearch && \
    chmod +x ./usearch12/bin/usearch

########## Install required python packages & edna_bp (eDNA bioinformatics pipeline)

WORKDIR /usr/src/app

COPY . .

RUN pip3 install -r requirements.txt --break-system-packages && pip3 install -e . --break-system-packages

########## Build workplace

WORKDIR /

RUN mkdir workplace

WORKDIR /workplace

CMD ["bash"]