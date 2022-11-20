#!/bin/env sh
julia -e "using Pkg; Pkg.activate(\"$1\"); Pkg.instantiate()"
julia -e "using AsciinemaGenerator, InteractiveUtils; cast_file(\"$1/main.jl\"; output_file=\"$1/main.cast\", mod=Main, tada=true)"
echo "DONE"
