#!/usr/bin/env bash

# SYNOPSIS
# ./render [-h] [-s sceneid]

# DESCRIPTION
# A script to render the given scene id using the raytracer binary. If there is
# no scene id, renders all the input scenes one by one.

# OPTIONS
# -h                Print the synopsis 
# -s sceneid        Render the input named "input<sceneid>.xml" in the input 
#                   directory

# EXAMPLES
# ./render -s 02    Renders the input02.xml
# ./render          Renders all the input scenes

INPUT_DIRECTORY="inputs"
OPTSTRING=":hs:"

usage() {
    echo "Usage: ${0} [-h] [-s sceneid]"
}

render() {
    ./raytracer "${INPUT_DIRECTORY}/input${1}.xml"
}

render_all() {
    for input_scene in "${INPUT_DIRECTORY}"/*.xml; do
        [ -e "${input_scene}" ] || continue;

        ./raytracer "${input_scene}"
    done
}

while getopts "${OPTSTRING}" option; do
    case ${option} in
        h)
            usage
            exit 0
        ;;

        s)
            if [[ ${OPTARG} == +([0-9]) ]]; then
                render ${OPTARG}
                exit 0

            else
                usage
                exit 1
            fi
        ;;

        *)
            usage
            exit 1
        ;;
    esac
done

render_all
exit 0
