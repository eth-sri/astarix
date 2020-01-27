# BUILD and TEST:
# $ sudo docker build -t astarix .
# RUN
# $ sudo docker run -it --rm --name astarix-container astarix

FROM ubuntu:18.04

# install prerequisites
# - build-essential: for make
# - libc-dev: contains argp
# - python3, python3-pip, pandas: for testing
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y \
		build-essential \
		libc-dev \
		python3 \
		python3-pip \
	&& apt-get clean && rm -rf /var/lib/apt/lists/*

# copy 
COPY . /astarix

# set working directory
WORKDIR /astarix

RUN pip3 install pandas

# compile and test
RUN make && \
	make test

