FROM ubuntu:20.04

RUN apt-get update -y
RUN apt-get install -y gcc 
RUN apt-get install -y libblas-dev liblapacke-dev 
RUN apt-get install -y make 
RUN apt-get install -y clang-format clang-tidy
RUN apt-get install -y valgrind
RUN apt-get install -y doxygen graphviz
RUN apt-get install -y python3 python3-pip
RUN pip3 install matplotlib


WORKDIR /app
COPY . .

ENV HOSTNAME=docker
RUN mkdir -p bin && \
  make all && \
  make run