#!/bin/env sh
julia --project=$1 -e "using Pkg; Pkg.instantiate()"
julia --project=$1 -e "using AsciinemaGenerator, InteractiveUtils; cast_file(\"$1/main.jl\"; output_file=\"$1/main.cast\", mod=Main, tada=true, height=30)"
echo "DONE"
