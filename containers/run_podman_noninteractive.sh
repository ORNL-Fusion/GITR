podman-hpc run -i -t --rm \
--name gitr_container \
--gpu \
-v /pscratch/sd/h/hayes/scaling-test:/host \
docker.io/stonecoldhughes/gitr:latest
