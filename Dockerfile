FROM python

#PREAMBLE

WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

RUN apt-get --assume-yes update \
	&& apt-get --assume-yes upgrade

#MAIN

RUN python setup.py instal
	&& rm -rf *.tgz *.tar *.zip \
	&& rm -rf /var/cache/apk/* \
	&& rm -rf /tmp/*

#ENVIRONMENT

ENC LC_ALL C
ENV LANG C
