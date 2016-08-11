#!/bin/bash

if [ $# -ne 1 ]; then
	echo "usage: batchRun.sh filelist"
	exit 1
fi

while read file; do

    ./Conudl -c $file

done < $1
