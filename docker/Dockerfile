FROM ubuntu:20.04

RUN apt-get update -y && \
  apt-get install -y gcc libblas-dev liblapacke-dev make vim clang-format clang-tidy

WORKDIR /app
COPY . .

ENV HOSTNAME=docker
RUN mkdir -p bin && \
  make all && \
  make run