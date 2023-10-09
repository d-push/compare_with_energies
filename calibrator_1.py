# -*- coding: utf-8 -*-
"""
Created on Tue Apr 24 16:56:08 2018

@author: Дмитрий
"""
import sys, glob, os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from scipy.signal import find_peaks

import argparse

import data_proc as dp

mpl.rcParams['agg.path.chunksize'] = 10000 #Enables large file plotting.
mpl.rc('text', usetex=True)
mpl.rcParams.update({'font.size': 24})
mpl.rcParams['text.latex.preamble']=r"\boldmath"

def run_av(array, av_size):
	average = np.zeros(len(array)-av_size)
	for i in range(0, len(array)-av_size):
		average[i] = np.average(array[i:i+av_size])
	return(average)

def parse_cmd_line():
	parser = argparse.ArgumentParser(description='Print acoustic data in a loop.')
	parser.add_argument('folder_ac', type=str, help="path to the folder containing subfolders with data to be matched")
	parser.add_argument('--en_chan', type=int, help="Number of the channel that contains signal from the energy detector.")
	parser.add_argument('--use_two_channels', '-utc', help="Use both channels from Rudnev in data processing: one for signals from diode and the second for signals from calorimeter.")
	parser.add_argument('--xlims', type=float, nargs=2, help="X-axis limits.")
	parser.add_argument('--ylims', type=float, nargs=2, help="Y-axis limits.")
	parser.add_argument('--old_osc', help="read data as for old version of oscilloscope (ms count). Otherwise, read as for new version (use time stamp with days and hours). This parameter is not compatible with --same.", action="store_true")
	parser.add_argument('--new_Rud', '-nR', action='store_true', help="data was recorded using new Rudnev module with new LV program.")
	parser.add_argument('--rem_out', type=float, help="Points distanced more than <rem_out> sigma values from the approximation, will not be taking into account when approximation coefficients will be calculated.")
	parser.add_argument('--debug', action='store_true', help="Print debug info.")
	parser.add_argument('--threshold', '-trs', type=float, help="Energies higher than this level are not taken into account for calibration. Default is 0.45.", default = 0.45)
	parser.add_argument('--plot_file_type', '-pft', type=str, help="File type for plot graph. Default is .png.", default = '.png')
	
	args = parser.parse_args()
	
	if args.en_chan is None:
		args.en_chan == 0
		
	return(args)

def check_main_maxima_number(filenames_ac, args):
	'''
	Определяет, где нужно обрезать waveform'у, чтобы в неё не попал 3-й большой максимум.
	'''
	#%%Constants.
	t_start_num = 5 # Отсуп на графике в начале waveform (в связи с наличием провала в начале).
	threshold_main_max = 0.5
	DIST_BETW_MAX = 26.74e-6 #Сдвиг между двумя самыми высокими максимумами.
	diff_threshold = 0.5
	
	#%%Initial parameters
	main_maxima = [] #Array with the maxima positions values.
	main_maxima_values = [] #Array with the values of the main maxima.
	
	for filename in filenames_ac:
		#%% Read data from file.
		wf_data = dp.read_bin_basing_on_args(filename, args)
		print(filename)
		if wf_data is None:
			continue
		else:
			two_channels = wf_data[0]
			dt, wfs = wf_data[1:]
			
		if len(wfs) > 1:
			wf = wfs[args.en_chan]
		else:
			wf = wfs[0]
			
		max_pos = np.argmax(wf[t_start_num:])
		main_maxima.append(max_pos) #Position (coordinate) of the main maximum
		main_maxima_values.append(wf[t_start_num+max_pos])
		
	main_maxima = np.array(main_maxima)
	main_maxima_values = np.array(main_maxima_values)
	
	#Отбрасываем максимумы ниже, чем threshold_main_max*<самый большой максимум>.
	main_maxima_values[main_maxima_values < (main_maxima_values[-1]*threshold_main_max)] = 0
	non_empty_files_indices = np.argwhere(main_maxima_values).flatten()
	main_maxima_values = main_maxima_values[non_empty_files_indices]
	main_maxima = main_maxima[non_empty_files_indices]
	main_maxima = np.sort(main_maxima)

	diffs = np.sort(np.diff(main_maxima))
	max_diff = diffs[-1]; second_max_diff = diffs[-2]
	
	args.xlims = [0]*2
	if (second_max_diff/max_diff < diff_threshold) or (max_diff < int(round(0.5*DIST_BETW_MAX/dt))):
		args.xlims[0] = 0; args.xlims[1] = main_maxima[-1]*dt + DIST_BETW_MAX/2.0
	else:
		args.xlims[0] = 0; args.xlims[1] = main_maxima[-1]*dt - DIST_BETW_MAX/2.0
		
	return(args)

def main(args):
	#%% Константы
	t_start_num = 5 # Отсуп на графике в начале waveform (в связи с наличием провала в начале).
	dt_between_max = 10e-6 #Нижняя граница расстояния между двумя максимумами на wf.
	SHIFT_LEFT = 37.75e-6 #Сдвиг левой границы области усреднения (1-го максимума) относительно главного максимума.
	SHIFT_RIGHT = 30.75e-6 #Сдвиг правой границы области усреднения (1-го максимума) относительно главного максимума.
	DIST_BETW_MAX = 26.74e-6 #Сдвиг между двумя самыми высокими максимумами.
	DIST_BETW_MAX_EPS = 0.05 #Ширина корридора между максимумами в долях расстояния.
	MAX_DIFF_BETWEEN_THE_MAXIMA = 3 #Найденные "большие" максимумы могут отличаться не более чем в это количество раз, иначе один из них - паразитный.
	dt_av = 1.0e-6 #Диапазон усреднения для бегущего среднего.
	dt_av_for_max_pick = 2e-6 #Диапазон усреднения для бегущего среднего для поиска главных максимумов.
	FON_SHIFT = 10e-6 #Сдвиг левой границы относительно shift left для подсчёта фона.
	SHIFT_LEFT_DIODE = 4e-6 #Сдвиг левой границы относительно максимума сигнала с диода для подсчёта фона.
	#fon_plateau_dur = 1e-6 #Длина части выборки, по которой подсчитывается ноль для сигнала с диода.
	V_to_J = 100.0 #Коэффициент пересчёта из вольт в Джоули.
	
	threshold = args.threshold #Энергии ниже этого уровня не учитываются при калибровке. Default is 0.45.
	
	#Recalculated
	
	old_path = os.getcwd()
	os.chdir(args.folder_ac)
	
	filenames_ac = glob.glob('*_*.bin')
	if len(filenames_ac) == 0:
		filenames_ac = glob.glob('*.bin')
	filenames_ac_info = [f.split("_") for f in filenames_ac]
	if len(filenames_ac_info[0]) > 1:
		filenames_ac_info = sorted(filenames_ac_info, key=lambda x: float(x[0]))
		filenames_ac = ["_".join(f) for f in filenames_ac_info]
	
	calibration = [] #Empty list for calibration data
	filenames_below_threshold = []
	#%%Main loop.
	two_high_maxima_files = []; one_high_maximum_files = []
	
	if args.xlims is None:
		print("Limits to find main maxima will be obtained automatically.")
		args = check_main_maxima_number(filenames_ac, args)
		print("Found limits are:")
		print(args.xlims)
		
	for filename in filenames_ac:
		#%% Read data from file.
		wf_data = dp.read_bin_basing_on_args(filename, args)
		if wf_data is None:
			continue
		else:
			two_channels = wf_data[0]
			dt, wfs = wf_data[1:]
			
		if len(wfs) > 1:
			wf = wfs[args.en_chan]
		
		#%%Cut according to args.xlims.
		wf = wf[int(round(args.xlims[0]/dt)):int(round(args.xlims[1]/dt))]
		
		#%% Поиск максимумов
		#en_wf_plot(t_array, wf[left_border:right_border+1], filename_to_save) #Построение графика.
		#wf_max = np.amax(wf[t_start_num:]) #Величина максимума
		#maxima = np.argwhere(wf == wf_max).flatten() # Координаты (отсч.) всех точек, значения в которых равны максмиальному.
		wf_averaged = dp.run_av(wf[t_start_num:], int(round(dt_av_for_max_pick/dt)) )
		init_zero = np.mean(wf_averaged[:int(round(FON_SHIFT/dt))]) #Нулевое приближение фона. Нужно для корректного определения амплитуд максимумов (и расчёта их отношения).
		wf_averaged = wf_averaged - init_zero
		
		maxima, _ = find_peaks(wf_averaged, distance=int(round(dt_between_max/dt)) )
		maxima = sorted(maxima, key=lambda x: wf_averaged[x])[-2:]
		
		dist_betw_max = abs(maxima[1] - maxima[0])*dt
		if args.debug:
			print(f'dist_betw_max = {dist_betw_max}, wf_averaged[maxima[0]] = {wf_averaged[maxima[0]]}, wf_averaged[maxima[1]] = {wf_averaged[maxima[1]]}')
		if (dist_betw_max < DIST_BETW_MAX*(1+DIST_BETW_MAX_EPS)) and (dist_betw_max > DIST_BETW_MAX*(1-DIST_BETW_MAX_EPS)) and (wf_averaged[maxima[0]]*MAX_DIFF_BETWEEN_THE_MAXIMA > wf_averaged[maxima[1]]):
			maxima = np.sort(maxima)
			max_coord = maxima[1] #Координата последнего из двух самых высоких максимумов.
			shift_left = SHIFT_LEFT; shift_right = SHIFT_RIGHT
			two_high_maxima_files.append(filename)
		else:
			max_coord = maxima[1] #Координата самого высокого максимума, считаем, что второй не поместился.
			shift_left = SHIFT_LEFT - DIST_BETW_MAX
			shift_right = SHIFT_RIGHT - DIST_BETW_MAX
			one_high_maximum_files.append(filename)
		
		if max_coord <= int(round(shift_left/dt)): #Проверка того, что максимум не слишком близко к началу выборки.
			continue
			
		#print(maxima[0]*dt*10**6, maxima[1]*dt*10**6, wf[maxima[0]], wf[maxima[1]], wf[maxima[0]]*MAX_DIFF_BETWEEN_THE_MAXIMA, wf[maxima[0]]*MAX_DIFF_BETWEEN_THE_MAXIMA > wf[maxima[1]])
		
		#%% Пересчёт сдвигов из секунд в отсчёты.
		shift_left_num = int(round(shift_left/dt))
		shift_right_num = int(round(shift_right/dt))
		fon_shift_num = int(round(FON_SHIFT/dt))
		
		#Установка левой границы, правой границы, и границы области для вычисления нулевого значения.
		left_border = max_coord - shift_left_num
		right_border = max_coord - shift_right_num
		if left_border - fon_shift_num >=0:
			fon_left_border = left_border - fon_shift_num
		else:
			fon_left_border = 0
	
		dnum_av = int(round(dt_av/dt)) #Ширина полосы для бегущего среднего в отсчётах.
		average = run_av(wf[left_border:right_border+1], dnum_av) #Вычисление бегущего среднего.
		en_V = np.amax(average) # Энергия в вольтах без учёта фона.
		fon = np.mean(wf[fon_left_border:left_border])
		energy = (en_V-fon)*V_to_J
	
		if energy >= threshold:
			if (len(wfs) > 1) and args.use_two_channels:
				wfd = wfs[1-args.en_chan]
				wfd = dp.run_av(wfd, window = int(round(dt_av/dt)) )
				max_coord_diode = np.argmax(wfd)
				shift_left_diode_num = int(round(SHIFT_LEFT_DIODE/dt))
				diode_fon_right_coord = max_coord_diode - shift_left_diode_num
				diode_fon_left_coord = diode_fon_right_coord - fon_shift_num
				fon_diode = np.mean(wfd[diode_fon_left_coord:diode_fon_right_coord])
				parrots = wfd[max_coord_diode] - fon_diode
			else:
				parrots = filename.split('_')[0]
			calibration.append([parrots, energy])
		else:
			filenames_below_threshold.append((filename, energy))
		
		###TEMP###
		#Graph plotting - points where averaging area start.
		'''
		plt.figure(figsize=(10.5, 9.0), dpi=300)
		plt.rcParams['text.latex.preamble'] = [r'\boldmath']
		plt.xlabel(r'\textbf{Energy, arb. un.}')
		plt.ylabel(r'\textbf{Energy, mJ}')
		t_array = np.linspace(200e-6, 250e-6, int(round((50e-6)/dt)))*10**6
		plt.plot(t_array, dp.run_av(wf, window = dnum_av)[int(round(200e-6/dt)):int(round(250e-6/dt))], '-', color='k', lw=1.5)
		plt.grid()
		plt.minorticks_on()
		plt.scatter([fon_left_border*dt*10**6, left_border*dt*10**6], wf[[fon_left_border, left_border]], color='g')
		plt.scatter((right_border+1)*dt*10**6, wf[right_border+1], color='b')
		plt.savefig(filename.split('.bin')[0]+'_temp.png', bbox_inches='tight')
		plt.close()
		'''
		###
	
	#%% Processing the obtained calibration array.
	parrots = np.zeros(len(calibration))
	energies = np.zeros(len(calibration))
	
	for i in range(0, len(calibration)):
		parrots[i] = calibration[i][0]
		energies[i] = calibration[i][1]
	
	#Linear approximation of the calibration data.
	print(parrots)
	print(energies)
	regr_params, cov = np.polyfit(parrots, energies, deg=1, cov=True)
	A, B = regr_params
	sigma_A, sigma_B = np.sqrt(np.diag(cov))
	print(A,B)
	print(sigma_A, sigma_B)
	
	#Remove the points that deviate from mean for more than 5 sigma.
	A_old = 0; B_old = 0
	if args.rem_out:
		while A != A_old or B != B_old:
			dev = np.abs(energies - (A*parrots+B)) #Модули отклонения от аппроксимирующей прямой.
			mean_dev = np.sqrt(np.mean(dev**2))
			#print(np.amin(dev))
			vpsk = np.argmax(dev)
			print(vpsk, np.amax(dev), parrots[vpsk], A*parrots[vpsk]+B, energies[vpsk])
			#print(np.mean(dev))
			print(f'mean_dev = {mean_dev}, max_dev = {np.amax(dev)}')
			#print(energies)
			#print(A*parrots + B)
			#print("dev:")
			#print(dev)
			print(len(parrots))
			parrots = parrots[dev < args.rem_out*mean_dev]
			energies = energies[dev < args.rem_out*mean_dev]
			print(len(parrots))
			print(dev[dev >= args.rem_out*mean_dev])
		
			A_old = A; B_old = B
			#Linear approximation of the new calibration data.
			regr_params, cov = np.polyfit(parrots, energies, deg=1, cov=True)
			A, B = regr_params
			sigma_A, sigma_B = np.sqrt(np.diag(cov))
			print(A,B)
			print(sigma_A, sigma_B)
	
	#Parameters for approximation plot.
	appr_x = np.array([np.amin(parrots), np.amax(parrots)])
	
	#Graph plotting
	plt.figure(figsize=(10.5, 9.0), dpi=300)
	plt.xlabel(r'\textbf{Signal from diode, V}')
	plt.ylabel(r'\textbf{Energy, mJ}')
	plt.plot(parrots, energies, '.', color='k', lw=1.5)
	plt.plot(appr_x, A*appr_x+B, '-', color='r', lw=1.5)
	plt.grid()
	plt.minorticks_on()
	plt.tick_params(axis='x', pad=7)
	plt.tick_params(axis='y', pad=5)
	s = r'$A = ({:.2f}\pm{:.2f})$'.format(A*1e3, sigma_A*1e3) + '$\cdot 10^{-3}$' + r'\textbf{mJ/un.}\,' + '\n$B = {:.3f}\pm{:.3f}$'.format(B, sigma_B) + r'\,\textbf{mJ}'
	plt.text(0.15, 0.78, s, bbox=dict(facecolor='none', edgecolor='black'), transform=plt.gcf().transFigure)
	filename_to_plot = os.path.basename(os.path.normpath(args.folder_ac)) + args.plot_file_type
	plt.savefig(filename_to_plot, bbox_inches='tight')
	plt.close()
	
	if args.debug:
		print('\nFilenames_below_threshold:')
		print(f'{filenames_below_threshold}\n')
		
		print(f"Обнаружено {len(two_high_maxima_files)} с двумя высокими максимумами и {len(one_high_maximum_files)} - с одним.")
		print("Файлы с двумя высокими максимумами:")
		print(two_high_maxima_files)
		print("Файлы с одним высоким максимумом:")
		print(one_high_maximum_files)
	
	print(f'Average_energy (energy > threshold, threshold = {threshold} mJ):')
	print(f'Energy = ({np.mean(energies[energies>threshold])} +- {np.std(energies[energies > threshold])}) mJ\n')
	
	RMS = np.sqrt( np.sum( (energies - A*parrots - B)**2 )/len(energies) )
	print(f'RMS deviation from the approximation: {RMS} mJ')
	
	os.chdir(old_path)

if __name__=='__main__':
	args = parse_cmd_line() #Command line argument parsing.
	print(args)
	main(args)