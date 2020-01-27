# BUILD and TEST:
# $ sudo docker build -t astarix .
# RUN
# $ sudo docker run -it --rm --name astarix-container astarix

FROM ubuntu:18.04

# install prerequisites
# - build-essential: for make
# - libc-dev: contains argp
# - zlib1g: zlib package
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
		build-essential \
		libc-dev \
		#zlib1g-dev \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# copy 
COPY . /astarix

# set working directory
WORKDIR /astarix

# compile and test
RUN make && \
	make test
