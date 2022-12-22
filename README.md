
Assuming in all cases we are in the project folder

## Run targets pipeline

``` r
# See documentation in targets package
targets::tar_make()
```

## Run in docker

  - Install Docker

  - Build docker image

<!-- end list -->

``` bash
docker build -t uoc:final .
```

  - Execute `run_in_docker.sh`

Notes:

  - `run_in_docker.sh` script assumes that the name of the image is
    `uoc:final` if the name of the image is different please modify the
    script accordingly.

  - The `_target` directory will be owned by root when the pipeline is
    executed using docker, this happens because Docker runs under root,
    as far as documentation indicates this cannot be changed. Remember
    to use `chown` to modify the ownership.