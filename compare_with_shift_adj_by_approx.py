import sys, os
import numpy as np
from scipy.optimize import curve_fit
import scipy.stats as sps
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import shutil as sh
import argparse

from scipy import ndimage


import data_proc as dp

def fit_function(x, a, b):
	return a * x + b

def calc_av_squared_deviation_from_trend(fit_function, a, b, X, Y):
	'''
	Calculate averaged squared deviation of a data set from linear trend (sigma). Only y deviation is considered.

	fit_function --- fitting function (linear).
	a, b --- slope and intercept
	X, Y --- arrays of X and Y data values.
	'''

	dev = np.std(Y - fit_function(X, a, b))

	return(dev)

def f_rem_out(fit_function, foldername_ac, a, b, a1, b1, energies_parrots, maxima_new, is_zero, N_out=5.0):
	av_dev_lim = N_out*calc_av_squared_deviation_from_trend(fit_function, a, b, energies_parrots, maxima_new)
	energies_parrots_rem_out = np.zeros_like(energies_parrots)
	maxima_new_rem_out = np.zeros_like(maxima_new)
	for i in range(len(energies_parrots)):
		if abs(maxima_new[i] - fit_function(energies_parrots[i], a, b)) < av_dev_lim:
			energies_parrots_rem_out[i] = energies_parrots[i]
			maxima_new_rem_out[i] = maxima_new[i]
	energies_parrots_rem_out = energies_parrots_rem_out[np.abs(energies_parrots_rem_out) > is_zero]
	maxima_new_rem_out = maxima_new_rem_out[np.abs(maxima_new_rem_out) > is_zero]

	#Построение нового графика (без выпадающих точек).
	calibration_file = foldername_ac+"/Калибровка_same_rem_out" + ".png"
	dp.en_wf_plot(energies_parrots_rem_out, maxima_new_rem_out, calibration_file, style = '.')
	##########NEW##############
	plt.figure(figsize=(10.5, 9.0), dpi=200)
	plt.plot(energies_parrots_rem_out, maxima_new_rem_out, '.', color='k', lw=1.5)
	plt.plot([np.amin(energies_parrots_rem_out),np.amax(energies_parrots_rem_out)], [a*np.amin(energies_parrots_rem_out)+b, a*np.amax(energies_parrots_rem_out)+b], '-', lw=1.5, label = '1')
	plt.plot([np.amin(energies_parrots_rem_out),np.amax(energies_parrots_rem_out)], [a1*np.amin(energies_parrots_rem_out)+b1, a1*np.amax(energies_parrots_rem_out)+b1],  '-', lw=1.5, label = '2')
	plt.grid()
	plt.legend()
	plt.savefig(calibration_file, bbox_inches='tight')
	plt.close()

def parse_cmd_line():
	parser = argparse.ArgumentParser(description='Match the data with corresponding pulse energies.')
	parser.add_argument('folder_en', type=str, help="path to the folder containing energy files or folders (in case of energies reading from oscillograms, aquired using Rudnev).")
	parser.add_argument('folder_ac', type=str, help="path to the folder containing subfolders with data to be matched")
	parser.add_argument('ft', type=str, help="type of the data to be matched. Possible values: ac, mode, spectr, interf, lum")
	parser.add_argument('ext', type=str, help="file extension. Possible values: .dat, .png, .tif, .bin")
	parser.add_argument('--old_osc', help="read data as for old version of oscilloscope (ms count). Otherwise, read as for new version (use time stamp with days and hours). This parameter is not compatible with --same.", action="store_true")
	parser.add_argument('--en_from_images', '-efIm', help="Extract energies from images, obtained by CCD after the telescope.", action="store_true")
	parser.add_argument('--en_from_Rudnev', '-efR', help="Extract energy from curves, aquired using new Rudnev USB oscilloscope.", action="store_true")
	parser.add_argument('--en_channel', type=int, help="Number of the channel, on what the energy was recorded. Parameter is used only for the new oscilloscope (when --en_from_rudnev is true). 0 or 1 is possible. Default is 0.")
	parser.add_argument('--try', dest='try_no_losses', help="check if there were no lost files, and if so, compare them one by one (without trying to correct lost)", action="store_true")
	parser.add_argument('--no_losses', dest='no_losses', help="assume that there were no lost files (can be used to match files without information about time.", action="store_true")
	parser.add_argument('--use_map', help="use breakdown map (for modes, luminescence and spectra_measurements)", action="store_true")
	parser.add_argument('--subtract_plane', help="subtract_plane (for modes, luminescence and spectra measurements)")
	parser.add_argument('--same', dest='same_computer', help="for data obtained at the same computer. Comparation directly by time. This parameter is not compatible with --old_osc. Automatically set --subtract_plane as False.", action="store_true")
	parser.add_argument('--inv', help="for inverted acoustics data", action="store_true")
	parser.add_argument('-ro','--rem_out', dest='rem_out', help="remove points that outlines trend more than N_out sigma values.", action="store_true")
	parser.add_argument('--ac_lims', type=float, nargs=2, help="time range whithin which the acoustic maxima should lie. Two values (minimum and maximum) should be provided in seconds.")
	parser.add_argument('--no_trig', help="Don't use strob column for energies files. Assume that the energy file recording started by turning on the strob (first acoustic waveform correspond to one of the first rows in the file with energies (or was obtained just before the first energy value was recorded)). Otherwise, the col_trig data will be used to define, when the recording of the acoustic started", action="store_true")
	parser.add_argument('--run_av', help="Use running average for acoustic data processing. Default averaging is by 11 points then by 20 points.", action="store_true")
	parser.add_argument('--shift_lims', type=int, nargs=2, help="shift (between first data file and first energy entry) maximum and minimum values (in shots). If not specified, limits are set to -3, 5")
	parser.add_argument('--mshift', type=int, help = "match the data using manually provided shift")
	parser.add_argument('--make_dict', help="Create dictionary using data names. Do NOT copy files.", action = "store_true")

	args = parser.parse_args()
	if args.old_osc and args.same_computer:
		print("\n --old_osc and --same parameters are not compatible with each other as old_osc doesn't allow to define time in hours unambiguously.\n")
	if args.ac_lims is not None:
		args.ac_lims = np.asarray(args.ac_lims)
	if args.shift_lims is None:
		args.shift_lims = [-3, 5]
	if args.same_computer:
		args.subtract_plane = False

	return(args)

def max_from_diode_wf(wf, dt, SHIFT_LEFT_DIODE = 3e-6, SHIFT_LEFT = 10e-6, dt_av = 1e-6):
	'''
	Calculate energy in parrots from the diode waveform taking zero level into account.
	'''
	
	wf = dp.run_av(wf, window = int(round(dt_av/dt)))
	max_coord_diode = np.argmax(wf)
	shift_left_diode_num = int(round(SHIFT_LEFT_DIODE/dt)); fon_shift_num = int(round(SHIFT_LEFT/dt))
	diode_fon_right_coord = max_coord_diode - shift_left_diode_num
	diode_fon_left_coord = diode_fon_right_coord - fon_shift_num
	fon_diode = np.mean(wf[diode_fon_left_coord:diode_fon_right_coord])
	energy = wf[max_coord_diode] - fon_diode
	
	return(energy)

def main():
	#%% Константы
	filetype_en = ".png" #Если энергия считывается с Руднева (args.en_from_Rudnev), автоматически задаётся filetype_en = ".bin"
	#shift_min = -100
	#shift_max = 100
	set_limit = 20 # Минимальная длина выборки, начиная с которой данные считаются успешно сопоставленными.
	use_fon = False # Использовать ли сигнал "фона" (земли) с диода. По умолчанию True.
	col_fon = 6 #"старая" версия - 8, "новая" - 6.
	col_trig = 8 #"старая" версия - 6, "новая" - 8.
	col_en = 9 # номер столбца с энергиями. По умолчанию - 9.
	use_trig = True # Использовать ли колонку с сигналом строба, или считать начало записи выборки в начале файла с энергиями.
	r2_coeff = 0.9 # Если в режиме без коррекции пропусков (--try) r^2 > r2_coeff*r^2 в обычном режиме сопоставления, то используется режим без коррекции пропусков.
	low_r2_threshold = 0.15 # Если r^2 ниже этого значения, сопоставление считается проведённым неуспешно, выполняется попытка автоматического поиска shift (по начальным временам сопоставляемых данных), и сопоставление при этом значении shift.
	bounds_inv = (np.array([-np.inf, -np.inf]), np.array([-1e-5, np.inf])) #Границы подбора параметров аппроксимации для "перевёрнутой" акустики.
	bounds_not_inv = (np.array([1e-5, -np.inf]), np.array([np.inf, np.inf])) #Границы подбора параметров аппроксимации для "не перевёрнутой" акустики.
	DELTA = 50 #[мс] Максимально допустимый сдвиг текущего кадра по времени относительно рассчётного времени. По умолчанию DELTA = 15 мс.
	is_zero = 1e-8 #Условие "обращения в ноль" значений float при проверках.
	is_zero_int = 1e-1 #Условие "обращения в ноль" значений float при проверках, когда сравниваются целые числа.
	N_out = 5.0 #Если данные отклоняются от линейной зависимости более, чем на N_out*sigma, то они выбрасываются из завершающего сопоставления.
	en_str_length_default = 17 #Длина строки в файле с энергиями (значение по умолчанию).
	bg_from_empty = None
	bd_map_path = "for_breakdown_map/bd_map.txt" #Путь к файлу с картой пробоев.
	bd_map_path_modes = "for_breakdown_map/bd_map_chinese.txt"
	empty_en_folders = []
	info_delim = '__'
	
	SHIFT_LEFT_DIODE = 3e-6; SHIFT_LEFT = 10e-6 #Параметры для поиска и вычитания нулевого уровня на сигнале с диода.
	DT_AV = 1e-6 #Длина окна усреднения бегущего среднего при определении максимума сигнала с диода.

	#%% Начальные параметры
	#no_gaps = False #Используется обычный метод сопоставления (с коррекцией пропусков) - False, или без коррекции пропусков (в режиме --try). Флаг, используемый в программе для вывода в файл режима сопоставления, начальное значение менять не нужно.

	#%% Parse command line arguments.
	args = parse_cmd_line()
	print(args)

	folder_en = args.folder_en
	folder_ac = args.folder_ac
	ft = args.ft
	ext = args.ext

	use_trig = not args.no_trig # Использовать ли колонку с сигналом строба, или считать начало записи выборки в начале файла с энергиями.

	if args.use_map:
		if args.ft=='spectr':
			bd_mult, bd_single = dp.read_bd_map(bd_map_path)
		elif args.ft == 'mode':
			bd_mult, bd_single = dp.read_bd_map(bd_map_path_modes)
		else:
			print("ERROR: could not recognise file type (ft)! Breakdown map will not be taken into account.")
			bd_mult = []
			bd_single = []
	else:
		bd_mult = []
		bd_single = []
	
	if args.en_from_Rudnev:
		args.use_fon = False
		filetype_en = ".bin"
		print("Fon (backgrond) usage will be disabled for energy calaculations because --en_from_Rudnev parameter has been passed. filetype_en was assigned to '.bin'.")

	#%% Print some reminders for user.
	print("\n!!!")
	print("use_fon = {}".format(use_fon))
	print("use_trig = {}".format(use_trig))
	print(f'run_av = {args.run_av}\n')
	print(f'DELTA = {DELTA}\n')

	print(folder_en)

	print("Use try_no_losses={}, no_losses={}, same_computer={}, old_osc={}, inv={}".format(args.try_no_losses, args.no_losses, args.same_computer, args.old_osc, args.inv))

	#%% Main.

	f_params = open(os.path.join(folder_ac, "Calibration_parameters.txt"), 'w')
	f_params.write("foldername_ac \t best_shift \t max_r^2 \t Number of the matched files\n")

	if ft=='mode' or ft=='lum': #Задание области, в которой будет производиться интегрирование (для люминесценции)
		area = list(int(element) for element in input("Введите границы области интегрирования, в которой наблюдается сигнал. По умолчанию: 0 1920 0 1200\n").split())
	else:
		area=(0,1920,0,1200) #По умолчанию.

	if len(area) == 0:
		area = (0,1920,0,1200)
		print(f"Использую значение по умолчанию. area = {area}")
	else:
		print(f"area = {area}")

	#Путь к папке для копирования сопоставленных мод.
	folder_to_write = os.path.abspath(os.path.join(folder_ac, os.pardir))+os.sep+folder_ac.split(os.sep)[-2] + "_matched" + os.sep

	#Удаление папки, если она существует, и создание её заново.
	if os.path.exists(folder_to_write):
		sh.rmtree(folder_to_write)
	os.makedirs(folder_to_write)

	#Список файлов с энергиями и файлов с акустикой, сортировано в алфавитном порядке.
	if not args.en_from_Rudnev and not args.en_from_images:
		filenames_en = [os.path.join(folder_en, f) for f in sorted(os.listdir(folder_en)) if f.endswith(filetype_en)]
	else:
		filenames_en = sorted(next(os.walk(folder_en))[1])
	foldernames_ac = sorted(next(os.walk(folder_ac))[1])

	#Вычисление длины строки в файле с энергиями.
	if not args.en_from_Rudnev and not args.en_from_images:
		en_str_length = dp.autodetect_en_line_length(filenames_en[0])
		if en_str_length == "FAIL":
			en_str_length = en_str_length_default #Присваиваем значение по умолчанию.
	else:
		en_str_length = None

	for filename_en in filenames_en:
		#Название файла (или папки - для нового Руднева) с энергиями (без расширения и пути).
		if not args.en_from_Rudnev and not args.en_from_images:
			filename_en_splitted = filename_en.split(os.sep)[-1]
			filename_en_splitted = filename_en_splitted.split(filetype_en)[0]
		else:
			filename_en_splitted = filename_en.split(os.sep)[-1]
			filename_en = os.path.join(folder_en, filename_en)

		for foldername_ac in foldernames_ac:
			foldername_ac_splitted = foldername_ac.split(os.sep)[-1]
			foldername_ac = os.path.join(folder_ac, foldername_ac)

			if filename_en_splitted == foldername_ac_splitted:
				if args.en_from_images:
					filenames_en_inner = [os.path.join(filename_en, f) for f in sorted(os.listdir(filename_en)) if f.endswith(filetype_en)]
					time_en = np.zeros(len(filenames_en_inner))
					energies = np.zeros(len(filenames_en_inner))
					for i, filename_en_inner in enumerate(filenames_en_inner):
						time_cur = filename_en_inner.split(os.sep)[-1].split(filetype_en)[0]
						time_cur = time_cur.split('__')[-1]
						time_en[i] = dp.calc_lum_time(time_cur)
						data_en, width, height = dp.read_modes_basing_on_ext(filename_en_inner, filetype_en, mult=False, dims=(1280,960))
						data_en = ndimage.median_filter(data_en, size = 2)
						data_en = dp.subtract_min_corner(data_en)
						energies[i] = np.sum(data_en)
					time_en = time_en*1000
				elif not args.en_from_Rudnev:
					try:
						time_en, energies, i_start, lc, signal, fon = dp.read_en_all_data(filename_en, line_length = en_str_length, col_en=col_en, col_fon=col_fon, col_trig=col_trig)
					except EmptyEnFile: #Если файл с энергиями оказался пустой.
						continue
					graph_file = foldername_ac+"/Фон.png"
					fon_i = np.linspace(0, len(fon), len(fon))
					dp.en_wf_plot(fon_i, fon, graph_file, style = '-')
					en_param_list=None
				else:
					filenames_en_inner = [os.path.join(filename_en, f) for f in sorted(os.listdir(filename_en)) if f.endswith(filetype_en)]
					time_en = np.zeros(len(filenames_en_inner))
					energies = np.zeros(len(filenames_en_inner))
					for i, filename_en_inner in enumerate(filenames_en_inner):
						wf_data = dp.read_bin_new_Rudnev(filename_en_inner)
						if wf_data is None:
							energies[i] = None
							continue
						two_channels = wf_data[0]
						if two_channels:
							T, dt, wf0, wf1 = wf_data[1:]
							wf0 = wf0[2:]
							wf1 = wf1[2:]
						else:
							T, dt, wf = wf_data[1:]
							wf = wf[2:]
						time_en_list = T.replace(',','.').split("-")
						time_en[i] = float(time_en_list[-1]) + 60*float(time_en_list[-2]) + 3600.0*float(time_en_list[-3])
						if two_channels:
							if args.en_channel == 1:
								energies[i] = max_from_diode_wf(wf1, dt, SHIFT_LEFT_DIODE = SHIFT_LEFT_DIODE, SHIFT_LEFT = SHIFT_LEFT, dt_av = DT_AV)
							else:
								energies[i] = max_from_diode_wf(wf0, dt, SHIFT_LEFT_DIODE = SHIFT_LEFT_DIODE, SHIFT_LEFT = SHIFT_LEFT, dt_av = DT_AV)
						else:
							energies[i] = max_from_diode_wf(wf, dt, SHIFT_LEFT_DIODE = SHIFT_LEFT_DIODE, SHIFT_LEFT = SHIFT_LEFT, dt_av = DT_AV)
					time_en = time_en*1000
					p = time_en.argsort()
					time_en = time_en[p]
					energies = energies[p]
					i_start = 0; lc = len(filenames_en_inner)
					en_param_list = [time_en, energies, i_start, lc]
					if lc <= 2:
						print("Folder with energies is (almost) empty")
						empty_en_folders.append(filename_en)
						continue

				#i_start = 0

				if args.old_osc:
					filenames_ac_times, filenames_ac_info = dp.make_file_list_to_compare(foldername_ac, ext)
				else:
					filenames_ac_times, filenames_ac_info = dp.make_file_list_to_compare_new_program(foldername_ac, ext)
					
				if len(filenames_ac_times) == 0:
					print("В папке с сопоставляемыми данными нет файлов искомого типа.")
					break

				#%% Calculate bd value from empty frames (for modes and spectra).
				if (ft == 'mode' or ft == 'spectr'):
					if ext == '.RAW':
						filenames_ac = [os.path.join(folder_ac, foldername_ac, '-'.join(f)+ext) for f in filenames_ac_info]
					else:
						filenames_ac = [os.path.join(folder_ac, foldername_ac, '__'.join(f)+ext) for f in filenames_ac_info]
					#bg_from_empty = dp.calc_bg_from_empty_frames(filenames_ac, ext, subtr_plane=True, bd_mult=bd_mult, bd_single=bd_single)
					#if bg_from_empty:
					#	bg_from_empty = bg_from_empty[0]
					#else:
					#	bg_from_empty=None

				###TEMP####
				print("\n!!!")
				print(foldername_ac_splitted)
				###TEMP####
				if 'kalibr' in foldername_ac_splitted or 'calibr' in foldername_ac_splitted:
					print("HERE1!!!\n")
					print(f'args.ac_lims = {args.ac_lims}')
					maxima = dp.read_maxima(foldername_ac, filenames_ac_times, filenames_ac_info, ext, area=area, fon_coeff=1.0, old_osc=args.old_osc, limit_max = True, inv=False, ac_lims=None, use_run_av=args.run_av, bg_from_empty=None, subtr_plane = args.subtract_plane, bd_mult=bd_mult, bd_single=bd_single)
				else:
					print("HERE2!!!\n")
					print(f'args.ac_lims = {args.ac_lims}')
					maxima = dp.read_maxima(foldername_ac, filenames_ac_times, filenames_ac_info, ext, area=area, fon_coeff=1.0, old_osc=args.old_osc, limit_max = False, inv=args.inv, ac_lims=args.ac_lims, use_run_av=args.run_av, bg_from_empty=None, subtr_plane=args.subtract_plane, bd_mult=bd_mult, bd_single=bd_single)

				if args.same_computer:
					#print((time_en*1000).astype(int))
					#print(filenames_ac_times)
					energies_parrots = []; maxima_new = []
					folder_ac_to_write = folder_to_write + filename_en_splitted + os.sep
					os.makedirs(folder_ac_to_write)
					for i, time_en_i in enumerate(time_en):
						for j, filename_ac_time in enumerate(filenames_ac_times):
							#print(f'time_en_i = {time_en_i}, filename_ac_time = {filename_ac_time}')
							#print(f'energies[i] = {energies[i]}, maxima[j] = {maxima[j]}')
							if abs(time_en_i - filename_ac_time) < DELTA and abs(energies[i]) > is_zero and abs(maxima[j]) > is_zero:
								if filename_ac_time < time_en_i:
									j += 1
									if j >= len(filenames_ac_times):
										continue
								energies_parrots.append(energies[i])
								maxima_new.append(maxima[j])
								
								if ext=='.RAW':
									filename_ac_short = "-".join(filenames_ac_info[j]) + ext
								else:
									filename_ac_short = info_delim.join(filenames_ac_info[j]) + ext
									
								filename_ac = os.path.join(foldername_ac, filename_ac_short)
								
								if not args.make_dict:
									filename_ac_new = os.path.join(folder_ac_to_write, str(energies[i]) + "_" + str(i) + "_" + filename_ac_short)
									try:
										sh.copyfile(filename_ac, filename_ac_new)
									except FileNotFoundError:
										if info_delim == '__':
											info_delim = '_'
										else:
											info_delim = '__'
									filename_ac_short = info_delim.join(filenames_ac_info[j]) + ext
									filename_ac = os.path.join(foldername_ac, filename_ac_short)
									sh.copyfile(filename_ac, filename_ac_new)
								else:
									small_dict = {'number': i, 'filename': filename_ac}
									en_dict[energies[i]].append(small_dict)
									
								break
					'''
					calibration = dp.compare_not_save_same_computer(filename_en, foldername_ac, ext=ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon=use_fon, DELTA=DELTA, line_length = en_str_length, en_param_list=en_param_list)
					### Подготовка данных для графика.
					j_start = 0
					energies_parrots = []
					maxima_new = []
					for i in range(0, len(filenames_ac_times)):
						for j in range(j_start, calibration.shape[0]):
							if abs(filenames_ac_times[i] - calibration[j, 1]) < is_zero_int and abs(calibration[j,0]) > is_zero and abs(maxima[i]) > is_zero:
								energies_parrots.append(calibration[j,0])
								maxima_new.append(maxima[i])
								j_start = j
								break
					'''
					energies_parrots = np.array(energies_parrots)
					maxima_new = np.array(maxima_new)

					### Построение графика (для случая "same_computer").
					calibration_file = foldername_ac+"/Калибровка_same" + ".png"
					dp.en_wf_plot(energies_parrots, maxima_new, calibration_file, style = '.')

					#########NEW#############
					a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new)[0:3]
					cur_r_2 = cur_r_2**2
					a1, b1 = sps.linregress(energies_parrots, maxima_new)[0:2]
					f_params.write(foldername_ac+" \t --- \t {:.4f} \t {:.0f} \t 'Same' method.\n".format(cur_r_2, len(energies_parrots)))

					if args.rem_out:
						f_rem_out(fit_function, foldername_ac, a, b, a1, b1, energies_parrots, maxima_new, is_zero, N_out=N_out)
					###########END_NEW#########
					#!!! Пока только строится новый график! Реально данные не исправляются.
					#########################
				else:
					best_shift = -1
					best_i = 0
					max_r_2 = -1
					num_matched = 0
					a = 0; b = 0
					best_shift_old = -1
					best_i_old = 0
					max_r_2_old = -1
					num_matched_old = 0
					a_old = 0; b_old = 0

					for shift in range(args.shift_lims[0], args.shift_lims[1]+1):
						energies_parrots = []
						maxima_new = []
						print("shift = {}".format(shift))

						calibrations = []; calibrations_old = [] #Инициализация. Нужно для проверки, будет одна или несколько калибровок для данного shift.
						if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc, line_length = en_str_length, en_param_list = en_param_list):
							energies_parrots_old = []
							maxima_new_old = []
							try:
								calibration = dp.compare_not_save_no_lost(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA, line_length = en_str_length, en_param_list=en_param_list) # No gaps
							# old -> general method (with gaps, old)
								calibrations_old = dp.compare_not_save_new_method(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA, line_length = en_str_length, en_param_list=en_param_list) #Для согласованности надо не забыть поменять метод ниже.
							except (IndexError):
								print("INDEX_ERROR in compare_not_save_new_method or compare_not_save_no_lost. Move to the next shift.")
								continue
							#For general (old) method.
							energies_parrots_list_old = []; maxima_new_list_old = []
							for calibration_old in calibrations_old:
								energies_parrots_old = []; maxima_new_old = []
								for j,filename_ac_info_calibr in enumerate(calibration_old[1]):
									for i,filename_ac_info_init in enumerate(filenames_ac_info):
										if filename_ac_info_calibr == filename_ac_info_init:
											energies_parrots_old.append(calibration_old[0][j])
											maxima_new_old.append(maxima[i])
											break

								energies_parrots_old = np.array(energies_parrots_old)
								maxima_new_old = np.array(maxima_new_old)
								energies_parrots_list_old.append(energies_parrots_old[:])
								maxima_new_list_old.append(maxima_new_old[:])
							
							if energies_parrots_list_old[0].size <= 1:
								print('Size of the matched energies array <= 1 for shift {}'.format(shift))
								continue

							for i, energies_parrots_old in enumerate(energies_parrots_list_old):
								a_old, b_old, cur_r_2_old = sps.linregress(energies_parrots_old, maxima_new_list_old[i])[0:3]
								cur_r_2_old = cur_r_2_old**2
								if cur_r_2_old > max_r_2_old and len(energies_parrots_old) >= set_limit:
									max_r_2_old = cur_r_2_old
									best_shift_old = shift
									best_i_old = i
									num_matched_old = len(energies_parrots_old)
						elif args.no_losses:
							try:
								calibration = dp.compare_not_save_no_lost(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA, line_length = en_str_length, en_param_list=en_param_list) # No gaps
							except(IndexError):
								print("INDEX_ERROR in compare_not_save_no_lost. Move to the next shift.")
								continue
						else:
							try:
								calibrations = dp.compare_not_save_new_method(filename_en, foldername_ac, ext, shift = shift, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, use_trig=use_trig, old_osc=args.old_osc, DELTA=DELTA, line_length = en_str_length, en_param_list=en_param_list) #Для согласованности надо не забыть поменять метод ниже.
							except (IndexError):
								print("INDEX_ERROR in compare_not_save_new_method. Move to the next shift.")
								continue
						if len(calibrations) == 0:
							calibrations.append(calibration) #Если была одна калибровка - пишем её в список calibrations для единообразия.
						# Калибровка для нового метода (если используется флаг --try и не были пропущены стробы.
						# Либо калибровка для старого (general) метода, если используется только он.
						energies_parrots_list = []; maxima_new_list = []
						for calibration in calibrations:
							energies_parrots = []; maxima_new = []
							for j,filename_ac_info_calibr in enumerate(calibration[1]):
								for i,filename_ac_info_init in enumerate(filenames_ac_info):
									if filename_ac_info_calibr == filename_ac_info_init:
										energies_parrots.append(calibration[0][j])
										maxima_new.append(maxima[i])
										break

							energies_parrots = np.array(energies_parrots)
							maxima_new = np.array(maxima_new)
							energies_parrots_list.append(energies_parrots[:])
							maxima_new_list.append(maxima_new[:])

						#%% Построение графиков.
						if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc, line_length = en_str_length, en_param_list = en_param_list):
							# No gaps and general method.
							for i, energies_parrots in enumerate(energies_parrots_list):
								calibration_file = foldername_ac+"/Калибровка_" + str(shift) + "_" + str(i) + "_no_gaps.png"
								dp.en_wf_plot(energies_parrots, maxima_new_list[i], calibration_file, style = '.', dpi=200)
								
							for i, energies_parrots_old in enumerate(energies_parrots_list_old):
								calibration_file = foldername_ac+"/Калибровка_" + str(shift) + "_" + str(i) + ".png"
								dp.en_wf_plot(energies_parrots_old, maxima_new_list_old[i], calibration_file, style = '.', dpi=200)
						else:
							# General method.
							for i, energies_parrots in enumerate(energies_parrots_list):
								calibration_file = foldername_ac+"/Калибровка_" + str(shift) + "_" + str(i) + ".png"
								dp.en_wf_plot(energies_parrots, maxima_new_list[i], calibration_file, style = '.', dpi=200)

						if energies_parrots_list[0].size <= 1:
							print('Size of the matched energies array <= 1 for shift {}'.format(shift))
							continue

						for i, energies_parrots in enumerate(energies_parrots_list):
							try:
								a, b, cur_r_2 = sps.linregress(energies_parrots, maxima_new_list[i])[0:3]
							except(ValueError):
								print("\n ERROR during the approximation. Probably, empty input arrays.\n")
								continue
							cur_r_2 = cur_r_2**2
							if cur_r_2 > max_r_2 and len(energies_parrots) >= set_limit:
								max_r_2 = cur_r_2
								best_shift = shift
								best_i = i
								num_matched = len(energies_parrots)
								
						#best_shift, max_r_2, a, b, num_matched = fit_and_check_shift(fit_function, energies_parrots, maxima_new, best_shift, shift, max_r_2, num_matched_old, set_limit, bounds=bounds_inv)

					#%% Если максимальный r2 меньше порога, запускаем автоматический поиск "правильного" сдвига.
					if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc, line_length = en_str_length, en_param_list = en_param_list):
						if max_r_2 < low_r2_threshold and max_r_2_old < low_r2_threshold:
							best_shift = dp.shift_search(filenames_ac_times, dt = 0.1, shift_border_min = -3, shift_border_max = 5)
							shift_corr = True
							print("Best_shift was recalculated automatically.")
							print("New best_shift = {}".format(best_shift))
						else:
							shift_corr = False
					else:
						if max_r_2 < low_r2_threshold:
							if args.no_losses:
								best_shift = 0
								shift_corr = True
								print("Best_shift was set to 0.")
								print("New best_shift = {}".format(best_shift))
							else:
								#best_shift = dp.shift_search(filenames_ac_times, dt = 0.1, shift_border_min = -3, shift_border_max = 5)
								best_shift = 1 #TEMP_for_ac. Normal case -> upper string.
								shift_corr = True
								print("Best_shift was recalculated automatically.")
								print("New best_shift = {}".format(best_shift))
						else:
							shift_corr = False

					# осталось запустить разметку интерферограмм с нужным best_shift
					folder_ac_to_write = os.path.join(folder_to_write, filename_en_splitted)
					os.makedirs(folder_ac_to_write)
					if args.try_no_losses == True and dp.check_if_there_were_lost(filename_en, foldername_ac, ext, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_fon = use_fon, old_osc=args.old_osc, line_length = en_str_length, en_param_list = en_param_list):
						print('max r_2 = {}, max r_2 old = {}, optimal_shift = {}, optimal_shift_old = {}, optimal_i_old = {}'.format(max_r_2, max_r_2_old, best_shift, best_shift_old, best_i_old))
						if max_r_2 > max_r_2_old * r2_coeff:
							print("Используется метод без пропусков.")
							dp.compare_no_lost(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA, line_length = en_str_length, en_param_list=en_param_list, make_dict=args.make_dict)
							no_gaps = True
						else:
							print("Используется обычный метод сопоставления.")
							dp.compare_new_method(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift_old, best_i = best_i_old, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA, line_length = en_str_length, en_param_list=en_param_list, make_dict=args.make_dict) #Для согласованности надо не забыть поменять метод выше.
							no_gaps = False
					elif args.no_losses:
						print('max r_2: {}, optimal_shift = {}'.format(max_r_2, best_shift))
						dp.compare_no_lost(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA, line_length = en_str_length, en_param_list=en_param_list, make_dict=args.make_dict)
						no_gaps = True
					else:
						print('max r_2: {}, optimal_shift = {}, best_i = {}'.format(max_r_2, best_shift, best_i))
						dp.compare_new_method(filename_en, foldername_ac, folder_ac_to_write, ext=ext, shift = best_shift, best_i = best_i_old, use_fon = use_fon, col_en=col_en, col_fon=col_fon, col_trig=col_trig, use_trig=use_trig, old_osc=args.old_osc, DELTA = DELTA, line_length = en_str_length, en_param_list=en_param_list, make_dict=args.make_dict) #Для согласованности надо не забыть поменять метод выше.
						no_gaps=False

					if shift_corr:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t Low r^2. Shift was found automatically\n".format(best_shift, max_r_2, num_matched))
					elif no_gaps:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t no gaps method\n".format(best_shift, max_r_2, num_matched))
					else:
						f_params.write(foldername_ac+" \t {:.0f} \t {:.4f} \t {:.0f} \t general method\n".format(best_shift, max_r_2, num_matched))

				#calibr_file = os.path.join(folder_ac, foldername_ac, "Calibr_" + str(shift) + ".png")
				#dp.en_wf_plot(energies_parrots, maxima_new, calibr_file, style = '.')

	f_params.close()
	
	if len(empty_en_folders) != 0:
		print("The following folders with energy are (almost empty):")
		for folder in empty_en_folders:
			print(folder)
		print("No matching has been performed for these folders.")

if __name__=='__main__':
	main()
