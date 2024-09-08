#!/bin/bash

SOURCE_FILES="\
main.cpp \
DataStructs.cpp \
rk4.cpp \
FluxFunctions.cpp \
RHSoperator.cpp"

mpic++ -g -O3 -pg -D_DOUBLE_ $SOURCE_FILES -Iincludes -o double.p

mpic++ -g -O3 -pg $SOURCE_FILES -Iincludes -o single.p
