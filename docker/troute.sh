#!/bin/bash
###################################################################################
# Helper script that allows you to execute the docker container like
# it was a local install of t-route e.g. 
# ~ python -m nwm_routing -V3 -f test/LowerColorado_TX/test_AnA.yaml
# Can be executed as
# ~ docker/troute.sh -V3 -f test/LowerColorado_TX/test_AnA.yaml
# Although the former would fail as you need to be in the same dir as the .yaml
# This script mounts the .yaml directory and makes it the workdir of the container
################################################################################### 

image_name="troute"
container_yaml_dir="/config"  # The directory in the container where the YAML file will be mounted
yaml_file_path=""
yaml_filename=""
parent_folder=""
# Loop through all arguments to find the .yaml file
for arg in "$@"; do
    if [[ $arg == *.yaml ]]; then
        # Get absolute path and just the filename
        yaml_file_path="$(realpath "$arg")"
        yaml_filename="$(basename "$arg")"
	parent_folder="$(dirname "$yaml_file_path")"
        break
    fi
done

# Check if yaml_file_path was found
if [ -z "$yaml_file_path" ]; then
    echo "No .yaml file provided in arguments."
    exit 1
fi

# Uncomment to print command
#echo "docker run -v ${parent_folder}:${container_yaml_dir}/ -w ${container_yaml_dir} ${image_name} ${@/$arg/$yaml_filename}"

# Mount the YAML file, set the working directory, and run the Docker container with all arguments (replacing the full YAML path with just its filename)
docker run -v "${parent_folder}:${container_yaml_dir}/" -w "${container_yaml_dir}" ${image_name} "${@/$arg/$yaml_filename}"

