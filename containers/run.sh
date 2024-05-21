docker run -i -t --rm \
--name cuda_12_container \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
--user $(id -u):$(id -g) \
--volume="/etc/group:/etc/group:ro" \
--volume="/etc/passwd:/etc/passwd:ro" \
--volume="/etc/shadow:/etc/shadow:ro" \
cuda_12_interactive
