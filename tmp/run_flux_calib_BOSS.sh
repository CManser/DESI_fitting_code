#!/bin/bash
d=`pwd` 

for file in /Users/christophermanser/Storage/PhD_files/DESI/BOSS_testing/wd_obs/uncalib_spec*-b?-*.dat
  do
    python3 flux_calib_BOSS.py $file # >> results_da2014.out
  done
