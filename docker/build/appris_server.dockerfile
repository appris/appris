# APPRIS-DATABASE ----

# our base image
FROM appris/core

# Copy confing file
COPY build/config_2017_08.v24.reduced.json /opt/appris/conf/.

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

