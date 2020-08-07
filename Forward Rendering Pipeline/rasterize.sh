#!/usr/bin/env bash

# SYNOPSIS
# ./render [-h] [-s sceneid]

# DESCRIPTION
# A script to rasterize the scene xml files in the directory with the given 
# directory id. If there is no directory id given, rasterizes all the input 
# scenes one by one.

# OPTIONS
# -h                    Print the synopsis 
# -l                    List the input directories with their number
# -d input_directory    Run the rasterize binary with all the input xml files 
#                       inside the input_directory

# EXAMPLES
# ./rasterize -d 1      Run the rasterize binary with the xml files inside the 
#                       directory number 1
# ./rasterize           Run the rasterize binary with all the xml files inside 
#                       all the directories

OPTIONS=d:hl
LONGOPTS=input_directory:,list_directories,help

declare -A INPUT_DIRECTORIES=(
    ["1"]="clipping_example"
    ["2"]="culling_disabled_inputs"
    ["3"]="culling_enabled_inputs"
    ["4"]="different_projection_type_example"
)

usage() {
    echo "Usage: ${0} [-d input_directory] [-h] [-l]"
}

rasterize_all() {
    for directory in "${!INPUT_DIRECTORIES[@]}"; do
        rasterize "${INPUT_DIRECTORIES[${directory}]}"
    done
}

rasterize() {
    for input_file in "${1}"/*; do
        [ -e "${input_file}" ] || continue;
        if [ -d "${input_file}" ]; then
            rasterize "${input_file}"
        elif [[ ${input_file} == *.xml ]]; then
            echo "Running rasterize with ${input_file}"
            ./rasterizer "${input_file}"
        fi
    done
}

case "${1}" in
    -d|--input_directory)
        if [[ ${2} == +([0-9]) ]]; then
            rasterize "${INPUT_DIRECTORIES[${2}]}"
            exit 0

        else
            usage
            exit 1
        fi
    ;;

    -h|--help)
        usage
        exit 0
    ;;

    -l|--list_directories)
        for directory in "${!INPUT_DIRECTORIES[@]}"; do
            echo "${directory}. ${INPUT_DIRECTORIES[${directory}]}"
        done
        exit 0
    ;;

    *)
        rasterize_all
        exit 0
    ;;
esac
