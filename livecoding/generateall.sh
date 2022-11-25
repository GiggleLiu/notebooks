#!/bin/env sh

for FILE in 1.basic  4.multipledispatch 2.array  5.performance  3.types 6.metaprogramming
do
    sh generate.sh $FILE
done
