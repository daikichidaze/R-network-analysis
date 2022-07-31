FROM python:3
USER root

RUN apt-get update

WORKDIR /workdir


COPY requirements.txt requirements.txt 

RUN pip install --upgrade pip
RUN pip install -r requirements.txt