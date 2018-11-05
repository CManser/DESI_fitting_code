# RUN IN PYTHON 3333333333333!!!!!!!!!!!!!!!!!!
# Outputs:
# 0 the object ID for desi data
# 1 T1
# 2 T1_err
# 3 logg1
# 4 logg1_err
# 5 T2
# 6 T2_err
# 7 logg2
# 8 logg2_err
# 9 rv
# 10 StN
#

import numpy as np
import astropy.io.fits as fits
ptl = 'petal9'
tmp = '/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/uncalibrated/'
f = open('dc17b_fluxcalib_{}.out'.format(ptl),'r')
f_out = open('/Users/christophermanser/Storage/PhD_files/DESI/flux_calibration_testing/fluxcalib_dc17b_{}.lst'.format(ptl), 'w')
lines = f.readlines()
for i in range(len(lines)):
  if lines[i][:2] == "/U":
    l = ""
    print(lines[i][:-1]) 
    l += lines[i][:-1][len(tmp):]
  if lines[i][:-1] == "First solution":
    tmp1 = lines[i+1][:-1].split(' ')
    tmp2 = lines[i+2][:-1].split(' ')
    tmp3 = lines[i+3][:-1].split(' ')
    l += ' ' + tmp1[3] + ' ' + tmp1[4] + ' ' + tmp2[3] + ' ' + tmp2[4]
  if lines[i][:-1] == "Second solution":
    tmp1 = lines[i+1][:-1].split(' ')
    tmp2 = lines[i+2][:-1].split(' ')
    tmp4 = lines[i+3][:-1].split(' ')
    l += ' ' + tmp1[3] + ' ' + tmp1[4] + ' ' + tmp2[3] + ' ' + tmp2[4] + ' ' + tmp3[2] + ' ' + tmp4[2] + '\n'
    print(l)
    f_out.write(l)
print("RUN IN PYTHON 3")