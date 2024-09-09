# Docker Version

**Step 1.** Download [Docker Desktop](https://docs.docker.com/get-docker/) on your own machine.

After installation, run `docker version` in terminal. If the version is displayed, it means the installation was successful.

**Step 2.** Clone the GitHub repo:
```sh
git clone https://github.com/ComputationalAgronomy/eDNA_bioinformatics_pipeline.git
```

**Step 3.** Build an image according to the *Dockerfile*.
```sh
docker build -t [ImageName] .
```

After the Dockerfile is successfully exported to an image, the [ImageName] should appear in the repo list if you use `docker image ls` to check.

**Step 4.** launch a new container from the Docker image that was just built:
```sh
docker run -it [ImageName]
```

If the launch is successful, your terminal should display something like `(base) root@93f4d3cf355f:/#`.

`(base)` indicates the conda environment you are using, where all dependent Python packages are installed (don't deactivate this!). `93f4d3cf355f` indicates the ID of this container.

**Step fin.** The Container is Ready to Work! Let's try [the first example](#simplist-example)!

#### Useful Docker Commands

(base) root@93f4d3cf355f:/# `exit` or `Ctrl+Z`: Exit the container.

`docker cp [container ID]:/path/to/file /host/destination/folder`: Copying files from Docker container to host.

`docker cp /path/to/file [container ID]:/container/destination/folder`: Copying files from host to Docker container.

`docker exec -it [container ID or Name] bash`: Enter a running container's shell.

`docker container ls -a`: List all containers

`docker stop [container ID or Name]`: Stop a running container.

`docker start [container ID or Name]`: Start a stopped container.

`docker rmi [ImageName]`: Remove a Docker image.

`docker container rm [container ID or Name]`: Remove a container.

`docker system prune (--force)`: Remove \<none> TAG images (be careful when using this command).