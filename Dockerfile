FROM alpine
# if you would like to include genomescope2, use "FROM dmolik/genomescope2" Although @molikd commited a pull request to genomescope2, so this may change

LABEL container="merfin" \ 
  about.summary="Improved variant filtering and polishing via k-mer validation" \
  about.home="https://github.com/arangrhie/merfin"

#PREAMBLE
WORKDIR /home/genomics
COPY . /home/genomics
RUN cd /home/genomics

# This is as small as possible alpine container for merfin (and meryl) 
# Stnadard update 
RUN apk update \
	&& apk upgrade

#MAIN
# minimum required packages
RUN apk add --no-cache --upgrade make bash gcc g++ git coreutils linux-headers perl libexecinfo-dev
# if using dmolik/genomescope2, you'd only need to add git and libexecinfo-dev 

#described install step
RUN git clone https://github.com/arangrhie/merfin.git \
  && cd merfin/src \
  && make -j 12

# put merfin in path
ENV PATH=${PATH}:/home/genomics/merfin/build/bin

#CLEANUP
RUN apk del git g++
# if using dmolik/genomescope2, only delete git

RUN rm -rf *.tgz *.tar *.zip \
	&& rm -rf /var/cache/apk/* \
	&& rm -rf /tmp/*
