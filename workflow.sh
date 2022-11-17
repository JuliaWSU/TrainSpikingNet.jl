#!/bin/bash
#
# Use test.jld2 to decide whether or not to purge data pre-existing data files, pertaining to a stale simulation run.
#
filename='test.jld2'
# Check the file is exists or not
if [ -f $filename ]; then
   rm test.txt
   echo "$filename is removed"
fi

./tsn.sh init ~/data_ben
./tsn.sh train gpu ~/data_ben
./tsn.sh test gpu ~/data_ben
./tsn.sh plot ~/data_ben/test.jld2
