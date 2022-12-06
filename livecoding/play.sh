#!/bin/env sh
DIRREL="$( dirname -- "${BASH_SOURCE[0]}"; )";   # Get the directory name
BASEDIR="$( realpath -e -- "$DIRREL"; )";    # Resolve its full path if need be
first_arg="$1"
DIR="$( realpath -e -- "$BASEDIR"; )/$first_arg";    # Resolve its full path if need be
shift
asciinema play "$DIR/main.cast" $@
echo "DONE"
