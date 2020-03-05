FROM ubuntu:latest


############################ Setup: OS #########################################
# Set up locales properly
RUN apt-get update && \
    apt-get install --yes --no-install-recommends locales && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

RUN echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
    locale-gen

ENV LC_ALL en_US.UTF-8
ENV LANG en_US.UTF-8
ENV LANGUAGE en_US.UTF-8

# Use bash as default shell, rather than sh
ENV SHELL /bin/bash

# run critical system updates ... one could also use the artful-curl pack
RUN apt-get update && \
    apt-get install --yes --no-install-recommends \
        wget bzip2 ca-certificates curl git zip unzip && \
    apt-get purge && \
    apt-get clean && \
    rm -rf /var/lib/apt/lists/*

# set up a user, to avoid running containers as root
ENV NB_USER container_user
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password \
    --gecos "Default user" \
    ${NB_USER}

# set up the user directory structure
RUN mkdir ${HOME}/tmp && \
    mkdir ${HOME}/bin
ENV PATH ${HOME}/bin:$PATH
################################################################################


############################ Setup: non-conda tools ############################
# install local scripts
RUN mkdir -p ${HOME}/scripts
COPY cromwell_monitor.sh /bin
################################################################################


############################ Deployment ########################################
# set wd to user home
WORKDIR ${HOME}

# set script as the default entry point
ENTRYPOINT ["/bin/bash", "/bin/cromwell_monitor.sh"]

# set the user
USER ${NB_USER}
################################################################################
