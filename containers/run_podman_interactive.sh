podman-hpc run -i -t --rm \
--name gitr_interactive \
--gpu \
-v $(pwd):/host \
docker.io/stonecoldhughes/gitr:gpu_gitr_interactive
