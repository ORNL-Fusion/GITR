docker run -i -t --rm \
--name gitr_container \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
gitr
