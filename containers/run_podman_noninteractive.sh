podman-hpc run -i -t --rm \
--name gitr_interactive \
--gpu \
-v $(pwd):/host \
--volume="/etc/group:/etc/group:ro" \
--volume="/etc/passwd:/etc/passwd:ro" \
--volume="/etc/shadow:/etc/shadow:ro" \
docker.io/stonecoldhughes/gitr:gpu_gitr_noninteractive
