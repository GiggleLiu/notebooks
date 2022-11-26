#!/bin/env sh
first_arg="$1"
shift
asciinema play "$first_arg/main.cast" $@
echo "DONE"
