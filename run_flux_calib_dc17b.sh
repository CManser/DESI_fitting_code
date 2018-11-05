#!/bin/bash
d=`pwd` 

for file in /Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/uncalibrated/uncalibrated_frame*-b0-*.dat
  do
    python3 flux_calib_dc17b.py $file # >> results_da2014.out
  done
