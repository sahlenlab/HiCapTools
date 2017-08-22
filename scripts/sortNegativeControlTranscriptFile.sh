#!/bin/sh

#Sort Transcript file with header

if [ -z "$1" ]; then
	echo "Script to sort the transcript file downloaded from UCSC to format required for HiCapTools"
	echo "Usage: sh sortTranscriptFile.sh /path/to/file_name"
	exit 1
fi

filename=$1
namewoext="${filename%.*}"
ext="${filename##*.}"

sortedfilename=$namewoext
sortedfilename="$sortedfilename.sorted."
sortedfilename=$sortedfilename$ext


(head -n 1 $filename && tail -n +2 $filename | sort -k13,13) > $sortedfilename
