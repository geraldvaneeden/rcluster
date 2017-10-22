FROM rocker/r-base
WORKDIR /usr/local/src/rscripts

# Setup SSH
RUN apt-get update && apt-get install -y openssh-server
RUN mkdir /var/run/sshd
COPY ssh/id_rsa /root/.ssh/id_rsa
COPY ssh/id_rsa.pub /root/.ssh/id_rsa.pub
COPY ssh/id_rsa.pub /root/.ssh/authorized_keys
COPY ssh/config /root/.ssh/config
RUN chmod 600 /root/.ssh/id_rsa && chmod 600 /root/.ssh/id_rsa.pub && chmod 600 /root/.ssh/authorized_keys

### Install packages
COPY src/packages.R /usr/local/src/rscripts/packages.R
RUN Rscript packages.R

### Inject src
COPY src /usr/local/src/rscripts
RUN R CMD INSTALL snowfall_1.84-6.1.tar.gz

CMD ["/user/sbin/sshd", "-D"]