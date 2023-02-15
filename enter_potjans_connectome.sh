#!/bin/bash
#
# Use test.jld2 to decide whether or not to purge data pre-existing data files, pertaining to a stale simulation run.
#

#mkdir -p ~/data_ben


#filename='~/data_ben/test.jld2'
# Check the file is exists or not
#if [ -f $filename ]; then
#   rm test.txt
#   echo "$filename is removed"
#fi

./tsn.sh init ${PWD}/src
#./tsn.sh train cpu ${PWD}/src
#./tsn.sh test cpu ${PWD}/src
#./tsn.sh plot ${PWD}/src/test.jld2
