#!/bin/bash
DIRREL="$( dirname -- "${BASH_SOURCE[0]}"; )";   # Get the directory name
BASEDIR="$( realpath -e -- "$DIRREL"; )";    # Resolve its full path if need be
for FILE in 1.basic  4.multipledispatch 2.array  5.performance  3.types 6.metaprogramming
do
    echo "Generating $FILE..."
    DIR="$( realpath -e -- "$BASEDIR"; )/$FILE";    # Resolve its full path if need be
    julia --project=$DIR -e "using Pkg; Pkg.instantiate()"
    julia --project=$DIR -e "using AsciinemaGenerator, InteractiveUtils; cast_file(\"$DIR/main.jl\"; output_file=\"$DIR/main.cast\", mod=Main, tada=true, height=30)"
done
