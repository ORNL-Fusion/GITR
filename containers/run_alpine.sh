docker run -i -t --rm \
--name alpine_interactive_container \
-v $(pwd):/host \
alpine_interactive_image
