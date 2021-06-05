Installation from the Docker repository
---------------------------------------

> __Note:__ For more detail, read the section *__Installation from the Docker repository__* in the text file [INSTALL
.md](https://apprisws.bioinfo.cnio.es/pub/docs/INSTALL.md).

Run the APPRIS/core container
-------------------------------
The next step is to run the *__appris/core__* image mounting the database directory and the _working directory_

```
$ docker run -itd \
    -v ${database_dir}:/opt/appris/db \
    -v ${working_dir}:/home/appris/ws \
    appris/code
```

> __Note:__ For more detail of the *__{database_dir}__* and *__{working_dir}__*, read the section in the text file [INSTALL.md](https://apprisws.bioinfo.cnio.es/pub/docs/INSTALL.md)

Now, we will be ablle to run APPRIS commands in a running container with _appris_ user.
```
$ docker exec --user appris -it {CONTAINER_ID} \
    bash -c "source /opt/appris/conf/apprisrc.docker && \
            appris_run_appris \
                -c /home/appris/ws/test.ws.env \
                -m firestar,matador3d,matador3d2,spade,corsair,corsair_alt,thump,crash,proteo,appris \
                -l info"
```

### Samples

For example, using the files in the "docker" samples. We run the following APPRIS commands with _appris_ user:

```
$ docker run -itd \
    -v /local/appris/db:/opt/appris/db \
    -v /local/appris/docker/ws:/home/appris/ws \
    appris/code
```

+ Run simple test (cached):

```
$ docker exec --user appris -it {CONTAINER_ID} \
    bash -c "source /opt/appris/conf/apprisrc.docker && \
            appris_run_appris \
                -c /home/appris/ws/test.ws.env \
                -m firestar,matador3d,matador3d2,spade,corsair,corsair_alt,thump,crash,proteo,appris \
                -l info"
```


