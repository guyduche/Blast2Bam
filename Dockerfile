FROM ubuntu:20.04

LABEL author="Aurelien Guy-Duche"

RUN apt update && apt install -y make gcc libxml2-dev xsltproc zlib1g-dev

ADD src /blast2bam

RUN cd blast2bam && make

WORKDIR /data

ENTRYPOINT ["blast2bam"]

CMD ["--help"]
