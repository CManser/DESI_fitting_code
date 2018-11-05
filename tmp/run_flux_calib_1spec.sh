#!/bin/bash
d=`pwd` 

for file in /Users/christophermanser/Storage/PhD_files/DESI/BOSS_testing/spPlate_flux_text/spPlate*.dat
  do
    python3 flux_calib_1spec.py $file # >> results_da2014.out
  done
