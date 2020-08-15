FROM debian:9.5-slim

LABEL "com.github.actions.name"="Dependency"
LABEL "com.github.actions.description"="Build a code dependency"
LABEL "com.github.actions.icon"="mic"
LABEL "com.github.actions.color"="purple"

LABEL "repository"="http://github.com/tyounkin/GITR-1"
LABEL "homepage"="http://github.com/actions"
LABEL "maintainer"="Tim <tyounkin19@gmail.com>"

ADD entrypoint.sh /entrypoint.sh
RUN chmod +x /entrypoint.sh
ENTRYPOINT ["/entrypoint.sh"]
