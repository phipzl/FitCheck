#!/bin/bash


for filename in *.raw; do
    output_filename="${filename%.raw}.mnc"
    rawtominc -transverse -input "$filename" "$output_filename" 4 4 3
done
