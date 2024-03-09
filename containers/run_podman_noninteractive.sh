podman-hpc run -i -t --rm \
--name gitr_container \
--gpu \
-v $(pwd):/host \
docker.io/stonecoldhughes/gitr:latest
