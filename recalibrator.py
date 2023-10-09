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
import argparse

def parse_cmd_line():
	parser = argparse.ArgumentParser(description='Print acoustic data in a loop.')
	parser.add_argument('folder_ac', type=str, help="path to the folder containing subfolders with data to be matched")
	parser.add_argument('ext', type=str, help="Recalibrated files extension.")
	parser.add_argument('--folder_calibr', '-fold_c', type=str, help="Folder needed to recalibrate by average.")
	parser.add_argument('--filter_substring', '-filt_s', type=str, help="Only folders containing substring in their names will be recalibrated.", default='')

	args = parser.parse_args()
		
	return(args)

def main(args):
	'''
	Recalibrator.
	'''
	
	#%%Ask user for calibration parameters.
	A,B = [float(element) for element in input("Введите А и B\n").split()]
	
	#%% Константы
	exceed_coeff = 1.25 #Если энергия превышает exceed*en_max, то папка перемещается в "подозрительные", калибровка пересчитывается по среднему.
	#deficit_coeff
	exceed_num = 3
	ext_ac = '.bin' #Расширение для калибровочных файлов с акустикой.
	correction_eps = 0.05 #Коридор (в процентах от нового значения), при выходе за пределы которого будет осуществлена рекалибровка.
	threshold = 0.4 #Граница "содержательных" (не пустых) кадров в мДж.
	
	calibr_col = 0 #Номер индекса в файле с калибровкой, в котором содержится значение в попугаях.
	
	#%% Расчёт калибровочных коэффициентов по данным файла с калибровкой. Поиск максимальной и минимальной энергии в калибровочном файле.
	if args.folder_calibr:
		calibr_filenames = [os.path.join(args.folder_calibr, f) for f in os.listdir(args.folder_calibr) if f.endswith(ext_ac)]
		parrots = np.array([f.split(os.sep)[-1].split('_')[calibr_col] for f in calibr_filenames], dtype = float)
		parrots_max = np.amax(parrots)
		en_max = A*parrots_max+B
		
		parrots_signif_kalibr = parrots[parrots > threshold] #Нужно для расчёта en_av.
		en_av = np.mean(A*parrots_signif_kalibr+B) #Нужно для расчёта нового коэффициента А, А = en_av/parrots_signif.
		
		exceed = exceed_coeff*en_max
	
	#A,B = [float(element) for element in input("Введите А и B\n").split()] #Read calibration coefficients.
	
	foldernames_ac = [f for f in next(os.walk(args.folder_ac))[1] if args.filter_substring in f]
	for foldername_ac in foldernames_ac:
		suspicious = False
		print(foldername_ac)
		foldername_ac = os.path.join(args.folder_ac, foldername_ac)
	
		filenames_ac = [f for f in sorted(os.listdir(foldername_ac)) if (f.endswith(args.ext) and '_' in f.split(os.sep)[-1])]
		if len(filenames_ac) == 0:
			print("Папка с акустическими файлами пуста.")
			continue 
		filenames_ac_info = sorted([f.split("_") for f in filenames_ac], key = lambda x: float(x[0]))
	
		parrots = np.array([filename_ac_info[0] for filename_ac_info in filenames_ac_info], dtype = float)
		energies = A*parrots+B
		parrots_signif = parrots[energies > threshold]
		
		print(f'Average_energy (energy > threshold, threshold = {threshold} mJ):')
		print(f'Energy = {np.mean(energies[energies>threshold])} +- {np.std(energies[energies > threshold])}')
	
		#en_av_new = np.mean(A*parrots_signif + B)
		#print("en_av = {}, en_av_new = {}".format(en_av, en_av_new))
	
		if args.folder_calibr:
			susp_en = np.zeros(exceed_num)
			exceed_count = 0
			for energy in energies:
				if energy > exceed:
					susp_en[exceed_count] = energy
					exceed_count += 1
					if exceed_count >= exceed_num:
						suspicious = True
						break
			print("susp_en (current maximum energy) = {}, exceed (energy threshold) = {}".format(susp_en, exceed))
	
		'''
		susp_en = np.amax(energies)
		if susp_en > exceed and recal_by_av or recal_by_av_all:
			suspicious = True
		'''
	
		if suspicious:
			print("WARNING: Suspicious calibration. Data will be recalibated by average.")
			print("Suspicious energies = {}, en_max = {}".format(susp_en, en_max))
			A_new = en_av/np.mean(parrots_signif)
			print("en_av = {}, np.mean(parrots_signif) = {}".format(en_av, np.mean(parrots_signif)))
			energies_new = A_new*parrots
			en_max_new = np.amax(energies_new)
			if abs(en_max_new - np.amax(susp_en)) < correction_eps*en_max_new:
				for i in range(0,len(energies)):
					filename_ac = os.path.join(foldername_ac, "_".join(filenames_ac_info[i]))
					filename_ac_new = os.path.join(foldername_ac, "{:.2f}_".format(energies[i])+"_".join(filenames_ac_info[i]))
					os.rename(filename_ac, filename_ac_new)
			else:
				for i in range(0, len(parrots)):
					filename_ac = os.path.join(foldername_ac, "_".join(filenames_ac_info[i]))
					filename_ac_new = os.path.join(foldername_ac, "{:.2f}_".format(energies_new[i])+"_".join(filenames_ac_info[i]))
					os.rename(filename_ac, filename_ac_new)
				foldername_ac_new = os.path.join(args.folder_ac, "Corrected_by_average")
				foldername_ac_new = os.path.join(foldername_ac_new, foldername_ac.split(os.sep)[-1]+"_Anew={:.2f}e-3".format(A_new*1e3))
				print(foldername_ac_new)
				#Удаление папки, если она существует.
				if os.path.exists(foldername_ac_new):
					sh.rmtree(foldername_ac_new)
				sh.move(foldername_ac, foldername_ac_new)
			suspicious = False
		else:
			for i in range(0,len(energies)):
				filename_ac = os.path.join(foldername_ac, "_".join(filenames_ac_info[i]))
				filename_ac_new = os.path.join(foldername_ac, "{:.2f}_".format(energies[i])+"_".join(filenames_ac_info[i]))
				os.rename(filename_ac, filename_ac_new)
	
	print("Done!")

if __name__=='__main__':
	args = parse_cmd_line() #Command line argument parsing.
	main(args)