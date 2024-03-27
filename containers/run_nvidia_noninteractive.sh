# -t option might not be needed...
docker run --rm \
--name gitr_container \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
gitr_noninteractive
