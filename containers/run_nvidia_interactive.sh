docker run -i -t --rm \
--name gitr_container \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
--user $(id -u):$(id -g) \
--volume="/etc/group:/etc/group:ro" \
--volume="/etc/passwd:/etc/passwd:ro" \
--volume="/etc/shadow:/etc/shadow:ro" \
gitr_interactive
