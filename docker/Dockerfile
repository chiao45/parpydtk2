FROM unifem/cht-coupler:mapper-dev
LABEL maintainer "Qiao Chen <benechiao@gmail.com>"

USER $DOCKER_USER
WORKDIR /tmp

RUN git clone --depth 1 -b parallel https://github.com/chiao45/parpydtk2.git && \
  cd parpydtk2 && \
  pip3 install . --user && \
  rm -rf /tmp/*

USER $DOCKER_USER
WORKDIR $DOCKER_HOME
