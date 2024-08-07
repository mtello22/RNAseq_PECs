FROM debian:bullseye-slim

LABEL image.author.name "Marco Tello"
LABEL image.authot.email "mtello@student.ubc.ca"

RUN apt-get update && apt-get install -y curl

RUN wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.9.0/salmon-1.9.0_linux_x86_64.tar.gz
RUN tar zxf salmon-1.9.0_linux_x86_64.tar.gz
RUN export PATH=$PATH:$PWD/salmon-1.9.0_linux_x86_64/bin
RUN mkdir -p /home/mtello/GitHub/RNAseq_PECs/data/transcriptome
