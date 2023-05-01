docker run -i -t --rm \
--name nvidia_interactive_container \
--runtime=nvidia \
--gpus all \
-v $(pwd):/host \
nvidia_interactive_image
