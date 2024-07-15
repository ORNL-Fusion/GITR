#DOCKER_BUILDKIT=1 \
docker build \
--file containers/Dockerfile_alpine \
--tag alpine_interactive_image .
