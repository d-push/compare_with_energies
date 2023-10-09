# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:56:08 2018

@author: Дмитрий
"""
import sys, os
import numpy as np
#import matplotlib.pyplot as plt
#import matplotlib as mpl
#import scipy.optimize as sp_opt
import data_proc as dp
import shutil as sh

folder_ac = sys.argv[1]
ext = sys.argv[2]
filter_substring = sys.argv[3]

foldernames_ac = [f for f in next(os.walk(folder_ac))[1] if filter_substring in f]
print(foldernames_ac)
for foldername_ac in foldernames_ac:
	print(foldername_ac)
	foldername_ac = os.path.join(folder_ac, foldername_ac)

	filenames_ac = [f for f in sorted(os.listdir(foldername_ac)) if (f.endswith(ext) and '_' in f.split(os.sep)[-1])]
	filenames_ac_info = [f.split("_") for f in filenames_ac]

	for i in range(0, len(filenames_ac_info)):
		filename_ac = os.path.join(foldername_ac, "_".join(filenames_ac_info[i]))
		filename_ac_new = os.path.join(foldername_ac, "_".join(filenames_ac_info[i][1:]))
		os.rename(filename_ac, filename_ac_new)

print("Done!")
