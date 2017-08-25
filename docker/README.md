Build the image
===============

Build the image of APPRIS-CORE
------------------------------
Now that you have your Dockerfile, you can build your image. The docker build command does the heavy-lifting of creating a docker image from a Dockerfile.

```
$ docker build -t appris/core -f build/appris_core.dockerfile .
```

Running APPRIS-Docker container
-------------------------------

```
$ docker run appris/appris
```

Running APPRIS-Docker container
-------------------------------

docker run -it -v /local/jmrodriguez/appris/db/:/mnt/ 278d76e17d86


