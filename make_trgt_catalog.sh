#!/bin/bash

# Check if the input file is provided as an argument
if [ -z "$1" ]; then
  echo "Usage: $0 <input_file>"
  exit 1
fi

# Input file from the first argument
input_file="$1"

# Use a temporary file to store modifications
temp_file=$(mktemp)
#echo "Processing input file: ${input_file}" >&2
# Process each line in the input file
while IFS=$(printf '\t') read -r chr start end col4 motifs; do
  # Create ID and formatted STRUC fields
  #echo "Creating ID and STRUC for motifs: ${motifs}">&2
  #echo "Processing line: ${chr}\t${start}\t${end}\t${col4}\t${motifs}">&2
  id="${chr}_${start}_${end}"
  struc="(${motifs})n"
  # Print the reformatted line to the temporary file
  echo -e "${chr}\t${start}\t${end}\tID=${id};MOTIFS=${motifs};STRUC=${struc}" >> "${temp_file}"
done < "${input_file}"

# Overwrite the original file with the modified content
mv "${temp_file}" "${input_file}"