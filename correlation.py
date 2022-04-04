# Correlation 
#
# Author: Joseph A'Hearn
# Created 03/09/2022
#
# This program tests for correlation between data and predictions from a few different models

import numpy as np 
import planet_ring_and_moon_data as prm 
import rings_and_resonances as rr  
import scipy.signal as sig 
import matplotlib.pyplot as plt 
import plot_assistant as pa 

def truncate_radscan(r_radscan, i_over_f_radscan, R, r_truncation=1.94):
	truncation_idx = int(((r_truncation - (1000 * r_radscan[0] / R)) / ((1000 * r_radscan[-1] / R) - (1000 * r_radscan[0] / R))) * len(r_radscan))
	r_truncated = r_radscan[:truncation_idx]
	i_truncated = i_over_f_radscan[:truncation_idx]
	return r_truncated, i_truncated

def quadratic_fit(r_truncated, i_truncated):
	# The method returns the polynomial coefficients ordered from low to high
	polynomial_coeff = np.polyfit(r_truncated, i_truncated, 2)
	#print("polynomial_coeff...", polynomial_coeff)
	return polynomial_coeff[0]*r_truncated**2 + polynomial_coeff[1]*r_truncated + polynomial_coeff[2]

def find_local_maxima(x, y):
	r_features = []
	#peak_idx, _ = sig.find_peaks(y, distance=15)
	peak_idx, _ = sig.find_peaks(y, distance=3, prominence=0.0045 / 1000)
	for i in range(len(peak_idx)):
		r_features.append(x[peak_idx[i]])
	prominences = sig.peak_prominences(y, peak_idx)[0]
	return r_features, prominences

def fourier_spectrum(t, y_t, freq_start=1.0, freq_end=1800.0, N=2**12):
	sample_freq = np.linspace(freq_start, freq_end, N)
	y_1_omega = []
	y_2_omega = []
	amplitude_fourier = []
	period_by_fourier = []
	for i in range(N):
		y_1_omega.append(np.sum(y_t * np.cos(sample_freq[i] * t)))
		y_2_omega.append(np.sum(y_t * np.sin(sample_freq[i] * t)))
		amplitude_fourier.append(2 * np.sqrt(y_1_omega[i]**2 + y_2_omega[i]**2) / len(t))
		period_by_fourier.append(2 * np.pi / sample_freq[i])
	return amplitude_fourier, sample_freq, period_by_fourier

def distance_to_nearest_feature():
	return r_res2feat

def distance_to_nearest_resonance(r_features, r_lmn, N_dim, r_min, r_max, m_min=2):
	# first, just make a list of the radial locations of the Lindblad resonances listed in all models
		# there may be a different total number from each model
		# maybe split into multiple lists: list of only l=m; list of l=m and l-m=2; list of l=m, l-m=2, and l-m=4

	r_Lres = []
	l_Lres = []
	m_Lres = []
	j_Lres = []
	idx_model_end = [] # to more easily separate the models later
	idx_model_end.append(0)
	for j in range(N_dim):
		for l in range(len(r_lmn)):
			for m in range(len(r_lmn[l])):
				if(((l - m) % 2) == 0): # only Lindblad resonances
					if(r_min < r_lmn[l][m][j] < r_max):
						r_Lres.append(r_lmn[l][m][j])
						l_Lres.append(l)
						m_Lres.append(m)
						j_Lres.append(j)
		idx_model_end.append(len(r_Lres)-1)

	r_feat2res = np.zeros((N_dim, len(r_features)))
	for i in range(len(r_features)):
		for j in range(N_dim):
			r_feat2res[j][i] = np.amin(r_features[i] - r_Lres[idx_model_end[j]:idx_model_end[j+1]])
			
				
			# incomplete. The above will likely need to be modified


	############

		#distance_between = np.absolute(r_features[i] - x2)
		#r_feat2res = np.amin(distance_between)

	return r_feat2res

def sum_of_squares(x):
	ss = 0
	for i in range(len(x)):
		ss += x[i]**2
	return ss 

def plot_features(planet, R, r_truncated, i_truncated, r_features, prominences, y_residuals, x_min=1.57, x_max=1.942, y_min=-0.03, y_max=0.15, rowspan=10, colspan=10, black=False):
	fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	#fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', 'mean ring normal quadratic fit residuals (' + r'$\times 1000$' + ')', x_min, x_max, y_min, y_max)
	fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', 'y', x_min, x_max, y_min, y_max)
	#ax.plot(1000 * r_truncated / R, i_truncated * 1000, c='k')
	for i in range(len(r_features)):
		y_top = y_max # default value 
		# find y_residuals idx closest to prominences i 
		for j in range(1, len(y_residuals)):
			if((np.sign(r_truncated[j] - r_features[i])) != np.sign(r_truncated[j-1] - r_features[i])):
				y_top = np.fmax(prominences[i], y_residuals[j]) * 1000
		ax.plot([1000 * r_features[i] / R, 1000 * r_features[i] / R], [y_min, y_top], linestyle='dashdot', color='r')
		ax.scatter(1000 * r_features[i] / R, prominences[i] * 1000, c='c')
	#ax.plot(1000 * r_truncated / R, y_fit * 1000, c='c')
	ax.plot(1000 * r_truncated / R, (y_residuals * 1000), c='g')
	ax.plot([x_min, x_max], [0, 0], linestyle='dotted', color='gray')
	ax.figure.savefig('feature_identification.pdf')
	plt.clf()

def plot_fourier(A, f, tau, rowspan=10, colspan=10, black=False):
	fig = plt.figure(figsize=(16,8))
	ax1 = plt.subplot2grid((10,10),(0, 0), rowspan=int(rowspan/2), colspan=colspan)
	fig, ax1 = pa.title_and_axes(fig, ax1, '', r'$f$' + ' [m' + r'$^{-1}$' + ']', 'Amplitude', f[0], f[-1], 0, np.amax(A) * 1.01)
	ax1.plot(f, A, c='r')
	
	ax2 = plt.subplot2grid((10,10),(5, 0), rowspan=int(rowspan/2), colspan=colspan)
	fig, ax2 = pa.title_and_axes(fig, ax2, '', r'$\lambda$' + ' [m]', 'Amplitude', np.amin(tau), np.amax(tau), 0, np.amax(A) * 1.01)
	ax2.plot(f, tau, c='r')
	ax2.figure.savefig('fourier_spectrum.pdf')
	plt.clf()

def plot_results():
	print('Plotting results...')

def main(filenames=[], planet='Uranus'):
	R, M, Omega, J2, J4, J6 = prm.import_planet(planet)
	# Run rings_as_resonances.py
	#r_rings, ring_names = prm.import_rings(planet)
	r_radscan, i_over_f_radscan = prm.import_radscan()
	r_truncated, i_truncated = truncate_radscan(r_radscan, i_over_f_radscan, R)
	y_fit = quadratic_fit(r_truncated, i_truncated)
	y_residuals = i_truncated - y_fit

	if(len(filenames) > 0):
		r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell, emm = rr.import_r_lmn(filenames, m_max=25)

	r_features, prominences = find_local_maxima(r_truncated, y_residuals)
	print(str(len(r_features)) + ' features found!')
	#plot_features(planet, R, r_truncated, i_truncated, r_features, prominences, y_residuals)
	A, f, tau = fourier_spectrum(r_truncated, y_residuals)
	plot_fourier(A, f, tau)
	
	if(len(filenames) > 0):
		r_feat2res = distance_to_nearest_resonance(r_features, r_lmn, len(filenames), r_truncated[0], r_truncated[-1])
		ss_feat2res = sum_of_squares(r_feat2res)
		for i in range(len(filenames)):
			print(filenames[i] + ': ' + str(ss_feat2res))


	#	for i in range(len(r_lmn[j])):
	#		r_res2feat = distance_to_nearest_feature(r_lmn[j], r_features)
	#	ss_res2feat = sum_of_squares(r_res2feat)
#
	#plot_results(ss_feat2res, ss_res2feat, filenames)

main()
#main(filenames=['Uranus_frequencies_thin1.txt', 'Uranus_frequencies_medium1.txt', 'Uranus_frequencies_thick1.txt', 'Uranus_frequencies_shallow4.txt', 'Uranus_frequencies_adiabatic1.txt'])