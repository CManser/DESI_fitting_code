#!/bin/bash
d=`pwd` 

for file in /Users/christophermanser/Storage/PhD_files/DESI/WD_ascii/*b-arm*
  do
    python3 DA_fitter.py $file # >> results_da2014.out
  done
