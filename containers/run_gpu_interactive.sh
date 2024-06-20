docker run -i -t --rm \
--name gitr_interactive \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
--user $(id -u):$(id -g) \
--volume="/etc/group:/etc/group:ro" \
--volume="/etc/passwd:/etc/passwd:ro" \
--volume="/etc/shadow:/etc/shadow:ro" \
stonecoldhughes/gitr:gpu_gitr_interactive
