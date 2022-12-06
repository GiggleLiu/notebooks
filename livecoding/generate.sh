#!/bin/bash
DIRREL="$( dirname -- "${BASH_SOURCE[0]}"; )";   # Get the directory name
BASEDIR="$( realpath -e -- "$DIRREL"; )";    # Resolve its full path if need be
DIR="$( realpath -e -- "$BASEDIR"; )/$1";    # Resolve its full path if need be
julia --project=$DIR -e "using Pkg; Pkg.instantiate()"
julia --project=$DIR -e "using AsciinemaGenerator, InteractiveUtils; cast_file(\"$DIR/main.jl\"; output_file=\"$DIR/main.cast\", mod=Main, tada=true, height=30)"
read -p "Run? [y/n] " -n 1 -r
echo    # (optional) move to a new line
if [[ ! $REPLY =~ ^[Yy]$ ]]
    source $BASEDIR/play.sh $1
then
    exit 1
fi