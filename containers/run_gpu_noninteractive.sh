# -t option might not be needed...
docker run --rm \
--name gitr_noninteractive \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
stonecoldhughes/gitr:gpu_gitr_noninteractive
