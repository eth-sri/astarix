# BUILD:
# $ sudo docker build -t astarix .
# RUN
# $ sudo docker run -it --rm --name astarix-container astarix

FROM ubuntu:18.04

# install prerequisites
# - build-essential: for sparsehash
# - automake1.11: for sparsehash
# - libc-dev: contains argp
# - zlib1g: zlib package
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
		git \
		build-essential \
		libc-dev \
		automake1.11 \
		zlib1g-dev \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# install sparsehash
RUN git clone https://github.com/sparsehash/sparsehash.git /sparsehash && \
	cd /sparsehash && \
	./configure && \
	make && \
	make install

# copy 
COPY . /astarix

# set working directory
WORKDIR /astarix

# compile and test
RUN make && \
	make test
