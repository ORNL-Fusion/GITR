docker run -i -t --rm \
--name alpine_interactive_container \
-v $(pwd):/host \
stonecoldhughes/gitr:cpu_gitr_alpine_interactive
