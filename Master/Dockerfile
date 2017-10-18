FROM rocker/r-base

### Install packages
RUN R --save -e "install.packages('Metrics', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('reshape', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('VIM', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('jomo', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('mice', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('missForest', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('ggplot2', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('gridExtra', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('doRNG', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('snow', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('doSNOW', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('foreach', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('MASS', repos = 'http://cran.mirror.ac.za')" \
	&& R --save -e "install.packages('sn', repos = 'http://cran.mirror.ac.za')"

RUN apt-get update && apt-get install -y openssh-server
RUN mkdir /var/run/sshd

#RUN R --save -e "system(\"echo 'PermitEmptyPasswords yes' >> /etc/ssh/sshd_config\")"\
 #&& R --save -e "system(\"echo 'PermitRootLogin yes' >> /etc/ssh/sshd_config\")"\
 #&& R --save -e "system(\"sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd\")"\
 #&& R --save -e "system(\"ssh-keygen -t rsa -N '' -f ~/.ssh/id_rsa\")"\
 #&& R --save -e "system(\"echo 'root:rclust' | chpasswd\")"\
 #&& R --save -e "system(\"echo 'export VISIBLE=now' >> /etc/profile\")"

#RUN sed -i 's/PermitRootLogin without-password/PermitRootLogin yes/' /etc/ssh/sshd_config
#RUN echo 'PermitRootLogin yes' >> /etc/ssh/sshd_config

# SSH login fix. Otherwise user is kicked off after login
#RUN sed 's@session\s*required\s*pam_loginuid.so@session optional pam_loginuid.so@g' -i /etc/pam.d/sshd

#RUN ssh-keygen -t rsa -N "" -f ~/.ssh/id_rsa
COPY ssh/id_rsa /root/.ssh/id_rsa
COPY ssh/id_rsa.pub /root/.ssh/id_rsa.pub
COPY ssh/id_rsa.pub /root/.ssh/authorized_keys
COPY ssh/config /root/.ssh/config

COPY src /usr/local/src/rscripts

WORKDIR /usr/local/src/rscripts

#CMD ["/bin/bash"]

CMD ["Rscript", "aither-knn.R"]
