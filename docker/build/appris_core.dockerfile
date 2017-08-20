# our base image
FROM ubuntu

# Update apt package
RUN apt-get update -y

# Install packages
RUN apt-get install git vim perl-doc nodejs npm -y

# Clone APPRIS code
RUN cd /opt && git clone https://github.com/appris/appris.git

# Install mysql-server
ENV HOME /root
ENV DEBIAN_FRONTEND noninteractive
RUN { \
        echo mysql-community-server mysql-community-server/data-dir select ''; \
        echo mysql-community-server mysql-community-server/root-pass password ''; \
        echo mysql-community-server mysql-community-server/re-root-pass password ''; \
        echo mysql-community-server mysql-community-server/remove-test-db select false; \
    } | debconf-set-selections \
&& apt-get install -y mysql-server

# APPRIS-CORE

# Setting up the environment variables
COPY build/appris.bashrc /tmp
RUN cat /tmp/appris.bashrc >> /root/.bashrc

# APPRIS-SERVER

# Install npm packages
RUN cd /opt/appris/server && npm install express

# tell the port number the container should expose
EXPOSE 5000

# Mount local directory to container 
#ADD mnt /mnt/
ADD db /mnt/

# run the application
#CMD service mysql start && tail -F /var/log/mysql/error.log
#CMD ["nodejs", "/opt/appris/server/scripts/web-server.js"]
