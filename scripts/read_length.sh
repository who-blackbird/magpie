#!/usr/bin/env bash
#Author: Alba Sanchis-Juan (as2635@cam.ac.uk)
#Compute reads length distribution from a fastq file

zcat $1 | awk 'NR%4 == 2 {lengths[length($0)]++} END {for (l in lengths) {print l, lengths[l]}}'