#!/bin/bash

if ! docker login --username ${GH_USERNAME} --password ${GH_TOKEN} ghcr.io; then
    echo "Error: Failed to login to ghcr.io. Please check your GH_USERNAME and GH_TOKEN env vars"
    exit 1
fi

if ! docker build -t ghcr.io/taddyb33/t-route-dev:${TAG} -f Dockerfile.troute .; then
    echo "Error: Failed to build Docker image. Please check your Dockerfile and build context."
    exit 1
fi

if ! docker push ghcr.io/taddyb33/t-route-dev:${TAG}; then
    echo "Error: Failed to push Docker image. Please check your TAG env var"
    exit 1
fi

echo "Successfully built and pushed Docker image with tag: ${TAG}"
