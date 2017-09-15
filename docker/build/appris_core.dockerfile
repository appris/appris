# APPRIS-CORE ----
# our base image
FROM ubuntu

# Update apt package
RUN apt-get update -y

# Install packages
RUN apt-get install -y acl git vim perl-doc nodejs npm 

# Clone APPRIS code
RUN cd /opt && git clone https://github.com/appris/appris.git
# delete packages because they produces some erros
RUN rm -rf /opt/appris/modules/lib/perl5/XML /opt/appris/modules/lib/perl5/x86_64-linux-gnu-thread-multi/auto/XML* /opt/appris/modules/lib/perl5/Module

# Install mysql-server
ENV DEBIAN_FRONTEND noninteractive
RUN { \
        echo mysql-community-server mysql-community-server/data-dir select ''; \
        echo mysql-community-server mysql-community-server/root-pass password ''; \
        echo mysql-community-server mysql-community-server/re-root-pass password ''; \
        echo mysql-community-server mysql-community-server/remove-test-db select false; \
    } | debconf-set-selections \
&& apt-get install -y mysql-server
COPY build/entrypoint.sh /sbin/entrypoint.sh
RUN chmod 755 /sbin/entrypoint.sh

# Install CPAN and Perl modules
RUN apt-get install -y libcam-pdf-perl libdbi-perl libdbd-mysql-perl libexpat1 libexpat1-dev libxml-sax-expat-perl gawk
RUN cpan -f -i XML::DOM
#COPY build/cpanConfig.pm /tmp/.
RUN cpan -f  -i Moose
RUN cpan -f  -i Module::Build
RUN cpan -f  -i Module::Runtime
RUN cpan -f  -i Module::Implementation
RUN cpan -f  -i IPC::Run

# Setting up the enviroment of 'root' and 'appris' users
COPY build/setup.root.sh /tmp/.
RUN cat "/tmp/setup.root.sh" >> /root/.bashrc
RUN chmod -R og+w /opt/appris/cache /opt/appris/tmp

RUN addgroup appris
RUN adduser --ingroup appris --disabled-password --gecos '' appris
#RUN apt-get install -y sudo
#RUN echo 'appris.appris' | adduser --ingroup appris --gecos '' appris
#RUN usermod -aG sudo appris
COPY build/setup.appris.sh /tmp/.
RUN cat "/tmp/setup.appris.sh" >> /home/appris/.bashrc

# Setting up the environment variables
#USER appris
#ENV HOME /home/appris
WORKDIR /home/appris
ENV APPRIS_HOME "/opt/appris"
ENV APPRIS_WORKSPACE "/home/appris/ws"
RUN mkdir "${APPRIS_WORKSPACE}"

# Setting up firestar DB
COPY build/FireDB_22Aug2013.sql.gz /tmp/.
RUN gzip -d /tmp/FireDB_22Aug2013.sql.gz
RUN /etc/init.d/mysql start && \
	mysql -h localhost -e 'create database FireDB; create user firedb@localhost; grant all privileges on FireDB.* to firedb@localhost' && \
	mysql FireDB -h localhost -u root < /tmp/FireDB_22Aug2013.sql && \
	rm /tmp/FireDB_22Aug2013.sql


# APPRIS-DATABASE ----

# APPRIS-SERVER ----
# Setting up the environment variables
#RUN /bin/bash -c "source ${APPRIS_HOME}/conf/apprisrc.WS"

# Install npm packages
#RUN cd /opt/appris/server && npm install express

# tell the port number the container should expose
#EXPOSE 3000

# Mount local directory to container 
#ADD db /mnt/

# Mysql user
#RUN /etc/init.d/mysql start && mysql -h localhost -e 'create user appris@localhost; grant all privileges on *.* to appris@localhost'

# COPY DATABASES for the code
#COPY myfirstapp/ /opt/appris/db/.

# run the application
#CMD service mysql start && tail -F /var/log/mysql/error.log
#CMD ["nodejs", "/opt/appris/server/scripts/web-server.js"]

