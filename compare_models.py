# Compare Models 
#
# Author: Joseph A'Hearn
# Created 01/05/2021
#
# This program compares the Uranus and Neptune models
#   that are ready to be input into GYRE
#   to Chris Mankovich's Saturn model

import numpy as np 
import matplotlib.pyplot as plt 
import plot_assistant as pa 
from scipy.optimize import curve_fit
#import convert as cv 
import matplotlib.ticker as ticker

def read_model(filename):
	data = np.loadtxt(filename, skiprows=1)
	r = data[:,1]
	m = data[:,2]
	#L = data[:,3]
	P = data[:,4]
	T = data[:,5]
	rho = data[:,6]
	nabla = data[:,7]
	N_sq = data[:,8]
	Gamma_1 = data[:,9]
	nabla_ad = data[:,10]
	delta = data[:,11]

	return r, m, P, T, rho, nabla, N_sq, Gamma_1, nabla_ad, delta 


def plot_models():
	#files = ['model.gyre']
	#files = ['Saturn_CM_model.gyre', 'Uranus_CM_model.gyre', 'U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'U3model.gyre', 'Uimodel.gyre', 'N1model.gyre', 'N2model.gyre', 'N3model.gyre']
	#files = ['Uranus_CM_model.gyre', 'U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'U3model.gyre', 'Uimodel.gyre']
	#files = ['model.gyre', 'U1model.gyre', 'U2model.gyre', 'Uimodel.gyre', 'N1model.gyre', 'N2model.gyre']
	files = ['Uranus_CM_model.gyre', 'U3model.gyre']
	#files = ['Ur_alt_model.gyre', 'U3_alt_model.gyre']
	
	fig, axs = plt.subplots(3, 3)
	axs[0, 0].set_title("interior mass")
	axs[0, 0].set_ylabel(r'$M_r$' + ' [g]')

	axs[1, 0].set_title("pressure")
	axs[1, 0].set_ylabel(r'$P$' + ' [dyn cm' + r'$^{-2}$' + ']')
	axs[1, 0].set_yscale('log')

	axs[2, 0].set_title("temperature")
	axs[2, 0].set_ylabel(r'$T$' + ' [K]')

	axs[0, 1].set_title("density")
	axs[0, 1].set_ylabel(r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']')
	#axs[0, 1].set_yscale('log')

	axs[1, 1].set_title('d ln' + r'$T$' + ' / d ln' + r'$P$')
	axs[1, 1].set_ylabel(r'$\nabla$')
	axs[1, 1].set_ylim(0.025, 0.4)

	axs[2, 1].set_title("Brunt-Vaisala frequency")
	axs[2, 1].set_ylabel(r'$N^2$' + ' [s' + r'$^{-2}$' + ']')
	axs[2, 1].set_ylim(0, 1.3)

	axs[0, 2].set_title('(d ln' + r'$P$' + ' / d ln' + r'$\rho$' + ')' + r'$_{ad}$')
	axs[0, 2].set_ylabel(r'$\Gamma_1$')
	axs[0, 2].set_ylim(1.5, 3)

	axs[1, 2].set_title('(d ln' + r'$T$' + ' / d ln' + r'$P$' + ')' + r'$_{ad}$')
	axs[1, 2].set_ylabel(r'$\nabla_{ad}$')
	axs[1, 2].set_ylim(0.025, 0.4)

	axs[2, 2].set_title('-(d ln' + r'$\rho$' + ' / d ln' + r'$T$' + ')' + r'$_{P}$')
	axs[2, 2].set_ylabel(r'$\delta$')
	axs[2, 2].set_ylim(0, 1)

	for file in files:
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(file)
		s=1
		
		axs[0, 0].scatter(r, m, s=s)
		axs[1, 0].scatter(r, P, s=s)
		axs[2, 0].scatter(r, T, s=s)
		axs[0, 1].scatter(r, rho, s=s, label=file[0:2])
		axs[1, 1].scatter(r, Delta, s=s)
		axs[2, 1].scatter(r, N_sq, s=s)
		axs[0, 2].scatter(r, Gamma_1, s=s)
		axs[1, 2].scatter(r, Delta_ad, s=s)
		axs[2, 2].scatter(r, delta, s=s)

		#if(file == 'U3model.gyre'):
		#	# polynomial fit 
		#	polynomial_coeff = np.polyfit(r, Gamma_1, 4)
		#	r_fit = np.linspace(r[0], r[-1], len(r))
		#	predict = np.poly1d(polynomial_coeff)
		#	#gamma_fit = (coeff1 * (r_fit**2)) + (coeff2 * r_fit) + coeff3
		#	gamma_fit = predict(r_fit)
		#	axs[0, 2].scatter(r_fit, gamma_fit, s=s, c='g')

	axs[0, 1].legend(loc=2, fontsize=6, bbox_to_anchor=(1,1))
	fig.tight_layout()
	fig.savefig("model_comparison.png")
	fig.clf() 

def plot_models2(black=False, dpi=200):
	if(black == True):
		c = 'w'
		#fig, axs = plt.subplots(2, 2, facecolor='black')
		fig, axs = plt.subplots(1, 3, facecolor='black', figsize=(10,5))
	else:
		c = 'k'
		#fig, axs = plt.subplots(2, 2)
		fig, axs = plt.subplots(1, 3)

	#files = ['Uranus_CM_model.gyre', 'Uthin4096model.gyre', 'Umedium4096model.gyre', 'Uthick4096model.gyre', 'U3model_with_new_derivative_equation.gyre', 'U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'Uimodel.gyre']
	files = ['Uranus_CM_model.gyre', 'Uthin4096model.gyre', 'Umedium4096model.gyre', 'Uthick4096model.gyre', 'U3model_with_new_derivative_equation.gyre']
	#files = ['Uranus_CM_model.gyre', 'U3model_with_new_derivative_equation.gyre']

	#files = ['Ur_alt_model.gyre', 'U3_alt_model.gyre']

	#fig, axs = plt.subplots(2, 2)

	#facecolor='black'
	#axs[0, 0].set_title("interior mass")
	#axs[0, 0].set_ylabel(r'$M_r/M$')
	#axs[0, 0].set_yscale('log')

	axs[0].set_title("density")
	axs[0].set_xlabel(r'$r/R$')
	axs[0].set_ylabel(r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']')
	#axs[0, 1].set_yscale('log')

	axs[1].set_title("pressure")
	axs[1].set_xlabel(r'$r/R$')
	axs[1].set_ylabel(r'$P$' + ' [TPa]')
	#axs[1, 0].set_yscale('log')

	axs[2].set_title("temperature")
	axs[2].set_xlabel(r'$r/R$')
	axs[2].set_ylabel(r'$T$' + ' [K]')
	#axs[1, 1].set_yscale('log')

	



	for file in files:
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(file)
		s=1
		
		if(file == 'Uranus_CM_model.gyre'):
			label = 'Mankovich'
		elif(file == 'U3model.gyre'):
			label = 'Scheibe et al. 2019'
		elif(file == 'U3model_with_new_derivative_equation.gyre'):
			label = 'Scheibe et al. 2019'
		elif(file == 'U_model.gyre'):
			label = 'Nettelmann et al. 2016'
		elif(file == 'U1model.gyre'):
			label = 'Nettelmann et al. 2013a'
		elif(file == 'U2model.gyre'):
			label = 'Nettelmann et al. 2013b'
		elif(file == 'Uimodel.gyre'):
			label = 'Bethkenhagen et al. 2017'
		elif(file == 'Uthin4096model.gyre'):
			label = 'Aramona & Mankovich - thin'
		elif(file == 'Umedium4096model.gyre'):
			label = 'Aramona & Mankovich - medium'
		elif(file == 'Uthick4096model.gyre'):
			label = 'Aramona & Mankovich - thick'
		M = np.amax(m)
		R = np.amax(r)
		#axs[0, 0].scatter(r / R, m / M, s=s)
		axs[0].scatter(r / R, rho, s=s, label=label)
		axs[1].scatter(r / R, P / 1.0E+11, s=s)
		axs[2].scatter(r / R, T, s=s)
		

	#axs[0].legend(loc=2, fontsize=7, bbox_to_anchor=(1,1))
	axs[0].legend(loc=0, fontsize=6)
	#fig.tight_layout()
	if(black == True):
		for j in range(3):
			axs[j].spines['bottom'].set_color('white')
			axs[j].spines['left'].set_color('white')
			axs[j].title.set_color('white')
			axs[j].yaxis.label.set_color('white')
			axs[j].xaxis.label.set_color('white')
			axs[j].tick_params(axis='x', colors='white')
			axs[j].tick_params(axis='y', colors='white')
		fig.savefig("model_comparison2.png", facecolor=fig.get_facecolor(), dpi=dpi, transparent=True)
	else:
		fig.savefig("model_comparison2.png", dpi=dpi)
	fig.clf() 

def plot_boundaries(planet='Uranus', black=False, rowspan=10, colspan=10):
	r_min = 0.9
	r_max = 1.0
	files = ['Uranus_CM_model.gyre', 'U3model.gyre']
	sources = ['Mankovich', 'Scheibe et al. 2019']

	# Density
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	fig, ax = pa.title_and_axes(fig, ax, 'Density', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', r_min, r_max, 0, 0.4)
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		ax.scatter(r / R, rho, label=sources[i])
	ax.legend(loc=0)
	pa.save_and_clear_plot(fig, ax, filename=planet + '_density_from_' + str(r_min) + '_to_' + str(r_max), black=black)

	# Pressure
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	fig, ax = pa.title_and_axes(fig, ax, 'Pressure', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$P$' + ' [TPa]', r_min, r_max, 0, 0.4)
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		ax.scatter(r / R, P / 1.0E+11, label=sources[i])
	ax.legend(loc=0)
	pa.save_and_clear_plot(fig, ax, filename=planet + '_pressure_from_' + str(r_min) + '_to_' + str(r_max), black=black)

	# Temperature
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	fig, ax = pa.title_and_axes(fig, ax, 'Temperature', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$T$' + ' [K]', r_min, r_max, 0)
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		ax.scatter(r / R, T, label=sources[i])
	ax.legend(loc=0)
	pa.save_and_clear_plot(fig, ax, filename=planet + '_temperature_from_' + str(r_min) + '_to_' + str(r_max), black=black)

def plot_derivatives():
	#files = ['model.gyre', 'U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'U3model.gyre', 'Uimodel.gyre', 'N1model.gyre', 'N2model.gyre', 'N3model.gyre']
	#files = ['Uranus_CM_model.gyre', 'U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'U3model.gyre', 'Uimodel.gyre']
	#files = ['Uranus_CM_model.gyre', 'Uimodel.gyre']
	files = ['Uranus_CM_model.gyre', 'U1model.gyre']
	fig, (ax1, ax2, ax3) = plt.subplots(3, 1, sharex=True)
	#ax1.set_title('(d ln' + r'$P$' + ' / d ln' + r'$\rho$' + ')' + r'$_{ad}$')
	ax1.set_ylabel(r'$\Gamma_1 = $' + '(d ln' + r'$P$' + ' / d ln' + r'$\rho$' + ')' + r'$_{ad}$')

	#ax2.set_title('(d ln' + r'$T$' + ' / d ln' + r'$P$' + ')' + r'$_{ad}$')
	ax2.set_ylabel(r'$\nabla_{ad} = $' + '(d ln' + r'$T$' + ' / d ln' + r'$P$' + ')' + r'$_{ad}$')

	#ax3.set_title('-(d ln' + r'$\rho$' + ' / d ln' + r'$T$' + ')' + r'$_{P}$')
	ax3.set_ylabel(r'$\delta = $' + '-(d ln' + r'$\rho$' + ' / d ln' + r'$T$' + ')' + r'$_{P}$')

	ax1.set_xlim([0,1])





	for file in files:
		r, m, P, T, rho, nabla, N_sq, Gamma_1, nabla_ad, delta = read_model(file)
		R = np.amax(r)
		

		if(file == 'Uranus_CM_model.gyre'):
			delta = -delta
		#if(file == 'model.gyre'):
		#	linewidth = 5
		#elif(file == 'U3model.gyre'):
		#	linewidth = 5
		#else:
		#	linewidth = 1	
		#ax1.plot(r / R, Gamma_1, linewidth=linewidth)
		#ax2.plot(r / R, Delta_ad, linewidth=linewidth)
		#ax3.plot(r / R, delta, label=file[0:2], linewidth=linewidth)
		ax1.scatter(r / R, Gamma_1)
		ax2.scatter(r / R, nabla_ad)
		ax3.scatter(r / R, delta, label=file[0:2])
		#if(file == 'model.gyre'):
		#	ax1.plot(r / R, Gamma_1, linewidth=linewidth)
		#	ax2.plot(r / R, Delta_ad, linewidth=linewidth)
		#	ax3.plot(r / R, delta, label=file[0:2], linewidth=linewidth)
		#elif(file == 'U3model.gyre'):
		#	ax1.plot(r / R, Gamma_1, linewidth=linewidth)
		#	ax2.plot(r / R, Delta_ad, linewidth=linewidth)
		#	ax3.plot(r / R, delta, label=file[0:2], linewidth=linewidth)

	# dashed lines with text
	ax1.plot([0, 1], [5 / 3, 5 / 3], linestyle='dashed', color='k')
	#ax1.text(0.2, 5 / 3, 'non-relativistic / monatomic ideal gas', verticalalignment='bottom', fontsize=10, color='k')
	ax1.plot([0, 1], [4 / 3, 4 / 3], linestyle='dashed', color='k')
	#ax1.text(0.2, 4 / 3, 'relativisitic / pure radiation gas', verticalalignment='bottom', fontsize=10, color='k')

	ax2.plot([0, 1], [ 0.4,  0.4], linestyle='dashed', color='k')
	ax3.plot([0, 1], [-1.5, -1.5], linestyle='dashed', color='k')

	ax3.legend(loc=3, fontsize=6)
	#fig.tight_layout()
	fig.savefig("derivative_comparison.png")
	fig.clf() 

def plot_density_and_BV(planet='Uranus', black=False, rowspan=10, colspan=10):
	if(planet == 'Uranus'):
		files = ['Uranus_CM_model.gyre', 'Uthin4096model.gyre', 'Umedium4096model.gyre', 'Uthick4096model.gyre', 'U3model_with_new_derivative_equation.gyre']
		#files = ['Uranus_CM_model.gyre', 'U1model.gyre']
		#files = ['U_model.gyre', 'U1model.gyre', 'U2model.gyre', 'U3model.gyre', 'Uimodel.gyre', 'Uranus_CM_model.gyre']
		#sources = ['Nettelmann et al. 2016', 'Nettelmann et al. 2013 slow', 'Nettelmann et al. 2013 fast', 'Scheibe et al. 2019', 'Bethkenhagen et al. 2018', 'Mankovich']
		#sources = ['Mankovich', 'Aramona & Mankovich - thin', 'Aramona & Mankovich - medium', 'Aramona & Mankovich - thick', 'Scheibe et al. 2019']
		sources = ['original', 'thin', 'medium', 'thick', 'Scheibe et al. (2019)']
		#sources = ['Mankovich', 'Nettelmann et al. 2013 slow']
	elif(planet == 'Neptune'):
		#files = ['N1model.gyre', 'N2model.gyre', 'N3model.gyre']
		#sources = ['Nettelmann et al. 2013 slow', 'Nettelmann et al. 2013 fast', 'Scheibe et al. 2019']
		files = ['neptune1637260138327702.gyre']
		sources = ['Mankovich']
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	#fig, ax = pa.title_and_axes(fig, ax, 'Density', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', 0, 1, 0)
	fig, ax = pa.title_and_axes(fig, ax, '', r'$r/R$', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', 0, 1, 0, 8.5)
	if(black == True):
		ax.set_title('density (solid bold) and Brunt-Vaisala frequency (dotted)', color='w', fontsize=20)
	else:
		ax.set_title('density (solid bold) and Brunt-Vaisala frequency (dotted)', color='k', fontsize=20)
	ax2 = ax.twinx()
	ax2.set_ylabel(r'$N$' + ' [s' + r'$^{-1}$' + ']', fontsize=18, color='c')
	ax2.set_ylim([0.0,0.005])
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		ax.scatter(r / R, rho,  label=sources[i])
		ax2.plot(  r / R, np.sqrt(N_sq), label=sources[i], linestyle='dotted')
	#ax.legend(loc=9, fontsize=20)
	#pa.save_and_clear_plot(fig, ax, filename=planet + '_density_and_BV_profiles', black=black)
	dpi=300
	if(black == True):
		ax.spines['bottom'].set_color('white')
		ax.spines['left'].set_color('white')
		ax.spines['top'].set_color('white')
		ax.spines['right'].set_color('white')
		ax.yaxis.label.set_color('white')
		ax.xaxis.label.set_color('white')
		ax.tick_params(axis='x', colors='white')
		ax.tick_params(axis='y', colors='white')
		ax.tick_params(axis='both', direction='in')
		ax.get_xaxis().tick_bottom()
		ax.get_yaxis().tick_left()

		ax2.spines['bottom'].set_color('white')
		ax2.spines['left'].set_color('white')
		ax2.spines['top'].set_color('white')
		ax2.spines['right'].set_color('white')
		ax2.yaxis.label.set_color('white')
		ax2.xaxis.label.set_color('white')
		ax2.tick_params(axis='x', colors='white')
		ax2.tick_params(axis='y', colors='white')
		ax2.tick_params(axis='both', direction='in')
		ax.figure.savefig(planet + '_density_and_BV_profiles.png', dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
	else:
		#ax.figure.savefig(planet + '_density_and_BV_profiles.png', dpi=dpi)
		ax.figure.savefig(planet + '_density_and_BV_profiles.pdf')
	plt.clf()

def plot_density_and_BV_for_paper(rowspan=12, colspan=12, black=False):
	files = ['Uranus_CM_model.gyre', 'Uthin4096model.gyre', 'Umedium4096model.gyre', 'Uthick4096model.gyre', 'U3model_with_new_derivative_equation.gyre', 'neptune1637260138327702.gyre']
	#sources = ['original', 'thin', 'medium', 'thick', 'Scheibe et al. (2019)', 'Mankovich']
	sources = ['shallow', 'thin', 'medium', 'thick', 'adiabatic (Scheibe et al. 2019)', 'Mankovich']
	#fig, (ax1, ax2) = plt.subplots(2, sharex=True)
	if(black):
		fig = plt.figure(figsize=(16,16), facecolor='black')
	else:
		fig = plt.figure(figsize=(16,16))
	ax1 = plt.subplot2grid((12,12),(0, 0), rowspan=int(rowspan / 2) - 1, colspan=colspan, fig=fig)
	ax2 = plt.subplot2grid((12,12),(6, 0), rowspan=int(rowspan / 2) - 1, colspan=colspan, fig=fig)
	ax1.set_ylabel(r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', fontsize=20)
	ax1.set_xlim([0,1])
	y_max = 8.5
	ax1.set_ylim([0,y_max])
	fig, ax2 = pa.title_and_axes(fig, ax2, '', r'$r/R$', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', 0, 1, 0, y_max)
	plt.rc('xtick', labelsize=20)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=20)    # fontsize of the tick labels
	#ax1.set_title('density (solid bold) and Brunt-Vaisala frequency (dotted)', fontsize=20)
	ax3 = ax1.twinx()
	ax3.set_ylabel(r'$N$' + ' [s' + r'$^{-1}$' + ']', fontsize=20)
	ax3.set_ylim([0.0,0.005])
	ax4 = ax2.twinx()
	ax4.set_ylabel(r'$N$' + ' [s' + r'$^{-1}$' + ']', fontsize=20)
	ax4.set_ylim([0.0,0.005])
	ax1.xaxis.set_major_locator(ticker.NullLocator())
	plt.xticks(fontsize=20)
	plt.yticks(fontsize=20)
	ax1.spines['top'].set_color('white')
	ax2.spines['top'].set_color('white')
	ax3.spines['top'].set_color('white')
	ax4.spines['top'].set_color('white')
	ax1.text(0.1,2, 'Uranus' , fontsize=32)
	ax2.text(0.1,2, 'Neptune', fontsize=32)
	if(black):
		c = ['r', 'c', 'orange', 'aliceblue', 'm']
	else:
		c = ['r', 'c', 'orange', 'b', 'm']
	spacing = 0.6
	fs = 16 
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		#print(files[i])
		#C = cv.inertia_check(rho, r, m[-1], R)
		if(i == 0):
			N_added = 8
			r_int = r[0]
			for j in range(N_added):
				r    = np.insert(r,    0, (j / N_added) * r_int)
				rho  = np.insert(rho,  0, rho[1])
				N_sq = np.insert(N_sq, 0, N_sq[1])
		if(i == 5):
			#ax2.scatter(r / R, rho, s=8)
			ax2.plot(r / R, rho, linewidth=4)
			ax4.plot(r / R, np.sqrt(N_sq), linestyle='dotted', linewidth=2)
		else:
			#ax1.scatter(r / R, rho, s=8, label=sources[i])
			ax1.plot(r / R, rho, linewidth=4, color=c[i])
			ax1.text(0.2, y_max - (i * spacing) - 0.5, sources[i],  color=c[i], verticalalignment='center', fontsize=fs)
			ax3.plot(r / R, np.sqrt(N_sq), color=c[i], linestyle='dotted', linewidth=2)
	#ax1.legend(loc=9, fontsize=16)
	#ax1.legend(bbox_to_anchor=[0.5, 0.95], loc='center', fontsize=20)
	if(black):
		axes = [ax1, ax2, ax3, ax4]
		for ax in axes:
			ax.spines['top'].set_visible(False)
			#ax.spines['right'].set_visible(False)
			ax.spines['bottom'].set_linewidth(0.5)
			ax.spines['left'].set_linewidth(0.5)
			ax.spines['bottom'].set_color('white')
			ax.spines['left'].set_color('white')
			ax.spines['right'].set_color('white')
			ax.title.set_color('white')
			ax.yaxis.label.set_color('white')
			ax.xaxis.label.set_color('white')
			ax.tick_params(axis='x', colors='white')
			ax.tick_params(axis='y', colors='white')
			ax.tick_params(axis='both', direction='in')
			ax.get_xaxis().tick_bottom()
		ax1.get_yaxis().tick_left()
		ax2.get_yaxis().tick_left()
		ax3.get_yaxis().tick_right()
		ax4.get_yaxis().tick_right()
		ax1.figure.savefig('UN_density_and_BV_profiles_black.png', facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax1.figure.savefig('UN_density_and_BV_profiles.pdf')
	plt.clf()

#Define the Gaussian function
def gauss(x,a,sigma):
	return a * np.exp(-(x)**2/(2*sigma**2))

def double_gauss(x, a, b, sigma1, sigma2):
	return (a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2)))

def gauss_plus_line(x, a, sigma, m, b):
	return (a * np.exp(-(x)**2/(2*sigma**2))) + (m * x) + b

def double_gauss_times_one_minus_r(x, a, b, sigma1, sigma2):
	return ((1 - x) * (((a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2)))) * ((1 - (1 * np.tanh(15 * (x - 0.8)))) / 2))) + 5.0E-04

def double_tanh_single_gauss_one_minus_r(x, a, sigma, b, c, ab=0.78, cb=0.16):
	#return (1 - x) * ((b * (1 - (1 * np.tanh(25 * (x - ab))))) + (c * (1 - (1 * np.tanh(15 * (x - cb)))))) * (a * np.exp(-(x)**2/(2*sigma**2)))
	return ((b * (1 - (1 * np.tanh(25 * (x - ab))))) + (c * (1 - (1 * np.tanh(15 * (x - cb)))))) * (a * np.exp(-(x)**2/(2*sigma**2)))

def fit_Net16(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(20 * (x - 0.803049915212053)))) + (c * (1 - np.tanh(20 * (x - 0.1795737879803282))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (2.1123524917657734 * (1 - x))

def fit_Net13s(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(20 * (x - 0.772460027704274)))) + (c * (1 - np.tanh(20 * (x - 0.175965042838437))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.8470085775875662 * (1 - x))

def fit_Net13f(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(20 * (x - 0.757626895535424)))) + (c * (1 - np.tanh(20 * (x - 0.14878409471968446))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.6383646099853342 * (1 - x))

def fit_Sch(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(20 * (x - 0.7490239925459978)))) + (c * (1 - np.tanh(20 * (x - 0.07178373865247327))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (2.471561135633503 * (1 - x))

def fit_Beth(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(20 * (x - 0.7842305293814895)))) + (c * (1 - np.tanh(20 * (x - 0.1813453920113811))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.7356022560861448 * (1 - x))

#def double_gauss_plus_line(x, a, b, sigma1, sigma2, m, c):
#	return (a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2))) + (m * x) + c
#
#def triple_gauss(x, a, b, c, sigma1, sigma2, sigma3):
#	return (a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2))) + (c * np.exp(-(x)**2/(2*sigma3**2)))

def find_boundaries(r, rho, R):
	for j in range(1, len(r)):
		if((rho[j-1] / rho[j]) > 1.5):
			if((r[j] / R) < 0.5): # core boundary
				cb = ((r[j] / R) + (r[j-1] / R)) / 2
			else: # atmosphere boundary
				ab = ((r[j] / R) + (r[j-1] / R)) / 2

				# slope of atmosphere
				delta_rho = rho[j] - rho[-1]
				delta_r = (r[-1] - r[j]) / R 
				slope = delta_rho / delta_r 
				#print('density at bottom of atmosphere: ' + str(rho[j]))
				#print('density at top of atmosphere: ' + str(rho[-1]))
				#print('slope: ' + str(slope))
	return ab, cb 

def plot_density(planet='Uranus', black=False, rowspan=10, colspan=10):
	if(planet == 'Uranus'):
		#files = ['U_modelref.gyre', 'U1modelref.gyre', 'U2modelref.gyre', 'U3modelref.gyre', 'Uimodelref.gyre', 'Uranus_CM_model.gyre']
		#sources = ['Nettelmann et al. 2016', 'Nettelmann et al. 2013 slow', 'Nettelmann et al. 2013 fast', 'Scheibe et al. 2019', 'Bethkenhagen et al. 2018', 'Mankovich']
		files = ['U_modelref.gyre', 'U1modelref.gyre', 'U2modelref.gyre', 'U3modelref.gyre', 'Uimodelref.gyre']
		sources = ['Nettelmann et al. 2016', 'Nettelmann et al. 2013 slow', 'Nettelmann et al. 2013 fast', 'Scheibe et al. 2019', 'Bethkenhagen et al. 2018']
	elif(planet == 'Neptune'):
		#files = ['N1model.gyre', 'N2model.gyre', 'N3model.gyre']
		#sources = ['Nettelmann et al. 2013 slow', 'Nettelmann et al. 2013 fast', 'Scheibe et al. 2019']
		files = ['neptune1637260138327702.gyre']
		sources = ['Mankovich']
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	fig, ax = pa.title_and_axes(fig, ax, 'Density', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', 0, 1, 0)
	#fig, ax = pa.title_and_axes(fig, ax, 'Density', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho$' + ' [g cm' + r'$^{-3}$' + ']', 0.7, 1, 0, 1) # to look at atmosphere
	for i in range(len(files)):
		r, m, P, T, rho, Delta, N_sq, Gamma_1, Delta_ad, delta = read_model(files[i])
		R = np.amax(r)
		ax.scatter(r / R, rho, s=10, label=sources[i])
		if(i < 5):	

			# find boundaries
			#ab, cb = find_boundaries(r, rho, R)

			# Gaussian fit
			x = r / R 
			y = rho

			#parameters, covariance = curve_fit(gauss, x, y)  
			#fit_y = gauss(x, parameters[0], parameters[1])
			#ss_res = np.sum((y - fit_y) ** 2) # residual sum of squares
			#ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
			#r2 = 1 - (ss_res / ss_tot) # r-squared
			#ax.plot(x, fit_y, c='k', linestyle='dashed', label='fit 1: ' + r'$R^2=$' + str(r2))
#
			#parameters, covariance = curve_fit(double_gauss_times_one_minus_r, x, y)  
			#fit_y = double_gauss_times_one_minus_r(x, parameters[0], parameters[1], parameters[2], parameters[3])

			#parameters, covariance = curve_fit(double_tanh_single_gauss_one_minus_r, x, y)  
			#fit_y = double_tanh_single_gauss_one_minus_r(x, parameters[0], parameters[1], parameters[2], parameters[3])

			if(i == 0):
				parameters, covariance = curve_fit(fit_Net16, x, y)  
				fit_y = fit_Net16(x, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(i == 1):
				parameters, covariance = curve_fit(fit_Net13s, x, y)  
				fit_y = fit_Net13s(x, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(i == 2):
				parameters, covariance = curve_fit(fit_Net13f, x, y)  
				fit_y = fit_Net13f(x, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(i == 3):
				parameters, covariance = curve_fit(fit_Sch, x, y)  
				fit_y = fit_Sch(x, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(i == 4):
				parameters, covariance = curve_fit(fit_Beth, x, y)  
				fit_y = fit_Beth(x, parameters[0], parameters[1], parameters[2], parameters[3])

			print(sources[i], parameters)

			#parameters, covariance = curve_fit(double_gauss, x, y)  
			#fit_y = double_gauss(x, parameters[0], parameters[1], parameters[2], parameters[3])
			ss_res = np.sum((y - fit_y) ** 2) # residual sum of squares
			ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
			r2 = 1 - (ss_res / ss_tot) # r-squared
			ax.plot(x, fit_y, linestyle='dashed', label=sources[i] + ' fit : ' + r'$R^2=$' + str(np.round(r2, decimals=4)))
#
			#parameters, covariance = curve_fit(gauss_plus_line, x, y)  
			#fit_y = gauss_plus_line(x, parameters[0], parameters[1], parameters[2], parameters[3])
			#ss_res = np.sum((y - fit_y) ** 2) # residual sum of squares
			#ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
			#r2 = 1 - (ss_res / ss_tot) # r-squared
			#ax.plot(x, fit_y, c='b', linestyle='dashed', label='fit 3: ' + r'$R^2=$' + str(r2))

			#parameters, covariance = curve_fit(double_gauss_plus_line, x, y)  
			#fit_y = double_gauss_plus_line(x, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
			#ss_res = np.sum((y - fit_y) ** 2) # residual sum of squares
			#ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
			#r2 = 1 - (ss_res / ss_tot) # r-squared
			#ax.scatter(x, fit_y, c='g', label='fit 4: ' + r'$R^2=$' + str(r2))
#
			#parameters, covariance = curve_fit(triple_gauss, x, y)  
			#fit_y = triple_gauss(x, parameters[0], parameters[1], parameters[2], parameters[3], parameters[4], parameters[5])
			#ss_res = np.sum((y - fit_y) ** 2) # residual sum of squares
			#ss_tot = np.sum((y - np.mean(y)) ** 2) # total sum of squares
			#r2 = 1 - (ss_res / ss_tot) # r-squared
			#ax.scatter(x, fit_y, c='chartreuse', label='fit 5: ' + r'$R^2=$' + str(r2))

			##plt.plot(xdata, ydata, 'o', label='data')
			#ax.scatter(xdata / R, fit_y, c='k', label='fit')

			#n = len(x)  
			#mean = np.sum(x*y)/n 
			#sigma = np.sum(y*(x-mean)**2)/n
#
			#popt,pcov = curve_fit(gaus,x,y,p0=[1,mean,sigma])
			#
			#plt.scatter(x ,gaus(x,*popt), c='k' ,label='fit')

			#initial_guess = [1,20,2,0]
			#initial_guess = [1,20,2,0]
			#popt, pcov = curve_fit(gaussian_func, x, y, p0=initial_guess)
			
			#a = np.max(rho)
			#sigma = 0.5
		
			#xplot = np.linspace(0, 1, 1000)
			#yplot = a * np.exp(-(xplot)**2/(2*sigma**2))
			#ax.scatter(xplot, yplot, c='k') 

			
			#ax.scatter(xplot,gaussian_func(xplot,*popt), c='k', label='fit')



	ax.legend(loc=0)
	pa.save_and_clear_plot(fig, ax, filename=planet + '_density_profiles', black=black)

def quick_plot(x, y, title=''):
	fig, ax = pa.initialize_plot(False)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '')
	ax.scatter(x, y)
	pa.save_and_clear_plot(fig, ax, filename='quick_plot' + title, black=False)

def quick_plot2(x1, x2, y, title=''):
	fig, ax = pa.initialize_plot(False)
	fig, ax = pa.title_and_axes(fig, ax, '', '', '')
	ax.scatter(x1, y)
	ax.scatter(x2, y)
	pa.save_and_clear_plot(fig, ax, filename='quick_plot' + title, black=False)

#def calculate_moment_of_inertia():
	#filenames = 

#plot_models()
#plot_models2(black=True)
#plot_derivatives()
#plot_density(planet='Uranus')
plot_density(planet='Neptune')
#plot_density_and_BV(black=False)
#plot_density_and_BV(planet='Uranus',  black=True)
#plot_density_and_BV(planet='Neptune', black=True)
plot_density_and_BV_for_paper(black=True)
#plot_boundaries()
#calculate_moment_of_inertia()