Installation from the Docker repository
=======================================

> __Note:__ For more detail, read the related section in the [INSTALL readme file](http://apprisws.bioinfo.cnio.es/pub/docs/INSTALL.md)

Build the image
===============

Build the image of APPRIS/core
------------------------------
Now that you have your Dockerfile, you can build your image. The docker build command does the heavy-lifting of creating a docker image from a Dockerfile.

```
$ docker build -t appris/core -f build/appris_core.dockerfile .
```
The following files are required:

+ build/appris_core.dockerfile
+ build/setup.sh
+ build/entrypoint.sh
+ build/FireDB_{data_version}.sql.gz


Running APPRIS-Docker container
===============================

Run the APPRIS/core container
-------------------------------
The next step is to run the *__appris/core__* image mounting the database directory and the _working directory_

```
$ docker run -itd \
    -v ${database_dir}:/opt/appris/db \
    -v ${working_dir}:/home/appris/ws \
    appris/code
```

> __Note:__ For more information about the *__{database_dir}__* and *__{working_dir}__*, read the related sections in the [INSTALL readme file](INSTALL.md#installation-from-the-docker-repository)

Now, we run APPRIS commands in a running container. We have to run the container with _appris_ user:

```
$ docker exec --user appris -it {CONTAINER_ID} \
    bash -c "source /opt/appris/conf/apprisrc.docker && \
            appris_run_appris \
                -c /home/appris/ws/test.ws.env \
                -m firestar,matador3d,matador3d2,spade,corsair,thump,crash,proteo,appris \
                -l info"
```
