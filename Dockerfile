FROM zwn97/theiasfm-docker-base

ADD . /theiasfm

RUN mkdir -p /theiasfm/build

WORKDIR /theiasfm/build

RUN cmake ..
RUN make
CMD make test

