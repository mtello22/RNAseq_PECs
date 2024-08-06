FROM debian:bullseye-slim

LABEL image.author.name "Marco Tello"
LABEL image.authot.email "mtello@student.ubc.ca"

RUN apt-get update && apt-get install -y curl

RUN curl -sSL https://github.com/COMBINE-lab/salmon/releases/download/v1.2.0/salmon-1.2.0_linux_x86_64.tar.gz | tar xz \
&& mv /salmon-*/bin/* /usr/bin/ \
&& mv /salmon-*/lib/* /usr/lib/

RUN mkdir -p /home/mtello/GitHub/RNAseq_PECs/data/transcriptome
