# Rings and resonances 
#
# Author: Joseph A'Hearn
# Created 02/02/2021
#
# This program plots the semi-major axes of the rings and of some resonances in the systems of the outer planets

import numpy as np 
import matplotlib.pyplot as plt 
import plot_assistant as pa 
import matplotlib.ticker as ticker
import get_gyre_output as gg 
import beta_check as bc 
import matrix as mx 
from collections import Counter
import matplotlib.markers as mk 
import planet_ring_and_moon_data as prm 
import compare_models as cm 
#print('matplotlib: {}'.format(matplotlib.__version__)) #---> matplotlib: 2.0.2 before, 3.3.4 now

def import_frequencies(M, R, planet, G=6.6743E-11):
	data = np.loadtxt('frequencies_' + planet + '.txt', skiprows=2)
	l = data[:,0] # oscillation degree ell 	
	#m = data[:,1] # oscillation order m 
	omega = data[:,2] # dimensionless units 
	beta = data[:,3] # splitting coefficient times 4pi

	freq = np.sqrt(G * M / (R**3)) * omega # multiplying the frequency by the scaling factor to get units of rad/s 
	beta /= (4*np.pi) # dividing by 4pi to normalize beta (should be between 0 and 1)

	return l, freq, beta  

def frequency_corrections(freq, m, beta, Omega): # from Marley & Porco 1993
	# freq input is the oscillation frequency on a non-rotating body
	# m is the oscillation order
	# Omega is the rotation frequency of the planet
	return (freq + (m * beta * Omega)) / m # what is being returned is the pattern frequency Omega_pat

def mean_motion(M, R, J2, J4, J6, r, G=6.6743E-11):
	b = R / r
	nsq = ((G * M) / (r**3)) * (1.0 + (1.5 * J2 * (b**2)) - (1.875 * J4 * (b**4)) + (2.1875 * J6 * (b**6)))
	return np.sqrt(nsq)

def horizontal_epicyclic(M, R, J2, J4, J6, r, G=6.6743E-11): 
	b = R / r 
	k1 = 1 - ((3 / 4) * J2 * (b**2)) + ((45 / 16) * J4 * (b**4)) - ((175 / 32) * J6 * (b**6))
	k2 = (- ((9 / 32) * (J2**2) * (b**4))) + ((135 / 64) * J2 * J4 * (b**6)) - ((27 / 128) * (J2**3) * (b**6)) 
	return np.sqrt((G * M) / (r**3)) * (k1 + k2)

def Lindblad_resonance_location(Omega_pat, Omega, kappa, m, r):
	func = (m * (Omega - Omega_pat)) + kappa
	r_res = 0
	for i in range(1, len(r)):
		if(np.sign(func[i]) != np.sign(func[i-1])):
			r_res = (r[i] + r[i-1]) / 2
	return r_res

def vertical_resonance_location(Omega_pat, Omega, kappa, m, r):
	mu = np.sqrt((2 * (Omega**2)) - (kappa**2))
	func = (m * (Omega - Omega_pat)) + mu 
	r_res = 0
	for i in range(1, len(r)):
		if(np.sign(func[i]) != np.sign(func[i-1])):
			r_res = (r[i] + r[i-1]) / 2
	return r_res

def find_corotation_resonances(filenames=[], planet='Uranus'):
	R, M, Omega, J2, J4, J6 = prm.import_planet(planet)
	r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell, emm = import_r_lmn(filenames, m_max=25) # only need Omega_pat
	r_min = R 
	r_max = R * 2.75
	r = np.linspace(r_min, r_max, 10000000)
	n = mean_motion(M, R, J2, J4, J6, r)
	r_c = []
	rce = []
	m_c = []
	for l in range(m_max + 1):
		m = int(l) # only doing Lindblad resonances with l = m
		for i in range(1, len(n)):
			if(np.sign(n[i] - Omega_pat[l][m]) != np.sign(n[i-1] - Omega_pat[l][m])):
				r_res = (r[i] + r[i-1]) / 2000 # finding the mean and converting to km
				print('l = ' + str(l) + ', m = ' + str(m) + ' ---> R_c = ' + str(np.round(r_res, decimals=1)) + ' km')
				r_c.append(r_res)
				m_c.append(m)
				record_data(planet, l, m, Omega_pat[l][m], r_res, Omega_err[l][m],  ((r_err[l][m] - r_lmn[l][m]) / 1000) + r_res)
			if(np.sign(n[i] - Omega_err[l][m]) != np.sign(n[i-1] - Omega_err[l][m])):
				rcerr = (r[i] + r[i-1]) / 2000 # finding the mean and converting to km
				print('l = ' + str(l) + ', m = ' + str(m) + ' ---> R_c_err = ' + str(np.round(rcerr, decimals=1)) + ' km')
				rce.append(rcerr) 
	for i in range(len(r_c)):
		print(str(m_c[i]) + ' &' + str(int(r_c[i])) + ' + ' + str(int(rce[i] - r_c[i])) + '\\')
	#cm.quick_plot(r_c, m_c, '_corotation_resonances')



def resonance_location(Omega_pat, M, R, r_min, r_max, lminusm, J2, J4, J6, m):
	r = np.linspace(R * r_min, R * r_max, 1000000)
	Omega = mean_motion(M, R, J2, J4, J6, r)
	kappa = horizontal_epicyclic(M, R, J2, J4, J6, r)
	if(lminusm % 2 == 0): 
		r_res = Lindblad_resonance_location(Omega_pat, Omega, kappa, m, r)
	else: 
		r_res = vertical_resonance_location(Omega_pat, Omega, kappa, m, r)
	return r_res 

def convert_to_pattern_frequency(P): # input is a pattern period in minutes
	return 2 * np.pi / (P * 60) # output is in rad/s

def convert_from_deg_day_to_rad_s(Omega):
	return Omega * np.pi / (180 * 8.64E+04)

def organize_data(ell, omega, beta):
	l_max = np.amax(ell) + 1 # first dimension of the 2D arrays
	# need to figure out the second dimenstion, the max number of n=0 modes for a single l
	c = Counter(ell)
	mc = c.most_common(1)#[1] # most common returns: [0] the 1 element that is most common (because 1 is in parentheses, and [1] the frequency of that element, which is what we want to know)
	N_l_max = mc[0][1]
	print('N_l_max: ' + str(N_l_max))
	#if(N_l_max == 1):
	#	l = np.zeros(l_max) 
	#	o = np.zeros(l_max)
	#	b = np.zeros(l_max)
	#	for i in range(l_max):
	#		l[i] = i
	#		for j in range(len(ell)):
	#			if(ell[j] == i):
	#				o[i] = omega[j]
	#				b[i] = beta[j]
	#else:
	l = np.zeros((l_max, N_l_max)) 
	o = np.zeros((l_max, N_l_max))
	b = np.zeros((l_max, N_l_max))
	for i in range(l_max):
		k = 0 # index of the second dimension
		for j in range(len(ell)):
			if(ell[j] == i):
				l[i][k] = i
				o[i][k] = omega[j]
				b[i][k] = beta[j]
				k += 1

	return l, o, b 

def xi_plots(path, filename='summary.h5', beta_plots=True, g_modes=False):
	pathname = path + filename
	j, ell, omega, beta = gg.relevant_data(pathname)
	tol = 0.01
	for i in range(len(j)):
		print(j[i], ell[i], omega[i], beta[i])
		beta_trapz = bc.main(j[i], ell[i], 0, path + 'eigs/j', plots=beta_plots, g_modes=g_modes)
		beta_error = np.absolute(beta[i] - beta_trapz) / beta[i]
		if(beta_error > tol):
			print('beta error: ' + str(beta_error))

def get_r_lmn(Omega_pat, Omega_err, r_lmn, r_err, freq, m, l, ell, beta, Omega, error, r_min, r_max, m_max, J2, J4, J6, M, R, dim_switch, gyre_data):
	if(gyre_data):
		if(dim_switch == False):
			Omega_pat[l][m] = frequency_corrections(freq[l], m, beta[l], Omega)
			domega = freq[l] * error 
			Omega_err[l][m] = frequency_corrections(freq[l] - domega, m, beta[l], Omega)		
		else:
			for i in range(N_dim):
				Omega_pat[l][m][i] = frequency_corrections(freq[l][i], m, beta[l][i], Omega) 
				Omega_err[l][m][i] = frequency_corrections(freq[l][i] * (1 - error), m, beta[l][i], Omega)
	else:
		enhancement_factor = np.sqrt(m / (m + 1))
		Omega_pat[l][m] = ((enhancement_factor * Omega_dyn * np.sqrt(l)) + (m * Omega)) / m
	if(dim_switch == False):
		r_lmn[l][m] = resonance_location(Omega_pat[l][m], M, R, r_min, r_max, l - m, J2, J4, J6, m)
		r_err[l][m] = resonance_location(Omega_err[l][m], M, R, r_min, r_max, l - m, J2, J4, J6, m)
		print(l, m, np.round(Omega_pat[l][m] * 180 * 86400 / np.pi), np.round(r_lmn[l][m] / 1000)) 
	else:
		for i in range(N_dim):
			r_lmn[l][m][i] = resonance_location(Omega_pat[l][m][i], M, R, r_min, r_max, l - m, J2, J4, J6, m)
			r_err[l][m][i] = resonance_location(Omega_err[l][m][i], M, R, r_min, r_max, l - m, J2, J4, J6, m)
			print(l, m, i, np.round(Omega_pat[l][m][i] * 180 * 86400 / np.pi), np.round(r_lmn[l][m][i] / 1000))
	return r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell

def calculate_r_lmn(planet, gyre_data, pathname, path, m_min, r_min, r_max, g_modes, beta_plots=True, G=6.6743E-11):
	R, M, Omega, J2, J4, J6 = prm.import_planet(planet)
	Omega_dyn = np.sqrt(G * M / (R**3))
	#print('omega_dyn: ' + str(Omega_dyn) + ' rad/s')
	error = (Omega / Omega_dyn)**2
	print('error: ' + str(error))
	if(gyre_data):
		if(g_modes):
			if((path == '/Users/josephahearn/Downloads/work/84_U_shallow/') or (path == '/Users/josephahearn/Downloads/work/92_shallow_including_1/')):
				j, ell, omega, beta = gg.data_for_gmodes2(pathname) # Use this for the shallow model
			else:
				j, ell, omega, beta = gg.data_for_gmodes(pathname)
		else:
			j, ell, omega, beta = gg.relevant_data(pathname)
			#j, ell, omega, beta = gg.relevant_data2(pathname) # Use this for the shallow model
		
		tol = 0.01
		for i in range(len(j)):
			print(j[i], ell[i], omega[i], beta[i])
			if(g_modes):
				beta_trapz = bc.main(j[i], ell[i], -1, path + 'eigs/j', plots=beta_plots, g_modes=g_modes)
			else:
				beta_trapz = bc.main(j[i], ell[i], 0, path + 'eigs/j', plots=beta_plots, g_modes=g_modes)
			beta_error = np.absolute(beta[i] - beta_trapz) / beta[i]
			if(beta_error > tol):
				print('beta error: ' + str(beta_error))
		# function that reorders ell, omega, and beta, so that ell counts from 0 to l_max
		ell, omega, beta = organize_data(ell, omega, beta)
		N_dim = len(ell[0])
		if(N_dim > 1):
			dim_switch = True 
		else:
			dim_switch = False

		# from this point on, I've added another dimension to ell, omega, and beta, which are now arrays of arrays 
		if(dim_switch):
			freq = np.zeros((len(omega), len(omega[0]))) 
			for i in range(len(freq)):
				freq[i] = Omega_dyn * omega[i] # multiplying the frequency by the scaling factor to get units of rad/s
		else:
			freq = Omega_dyn * omega # multiplying the frequency by the scaling factor to get units of rad/s
	else:
		dim_switch = False
		N = 26 #12
		ell = np.zeros(N)
		for i in range(N):
			ell[i] = i
	
	m_max = int(np.amax(ell))

	if(dim_switch == False):
		Omega_pat = np.zeros((m_max + 1, m_max + 1)) # the first two rows and columns will be zeros for l=m=0,1 
		Omega_err = np.zeros((m_max + 1, m_max + 1))
		r_lmn = np.zeros((m_max + 1, m_max + 1))
		r_err = np.zeros((m_max + 1, m_max + 1))
	else: 
		Omega_pat = np.zeros((m_max + 1, m_max + 1, N_dim))
		Omega_err = np.zeros((m_max + 1, m_max + 1, N_dim))
		r_lmn = np.zeros((m_max + 1, m_max + 1, N_dim))
		r_err = np.zeros((m_max + 1, m_max + 1, N_dim))

	if(g_modes):
		for l in range(1, m_max + 1):
			m = l
			r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell = get_r_lmn(Omega_pat, Omega_err, r_lmn, r_err, freq, m, l, ell, beta, Omega, error, r_min, r_max, m_max, J2, J4, J6, M, R, dim_switch, gyre_data)
	else:
		for l in range(2, m_max + 1):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1):
				r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell = get_r_lmn(Omega_pat, Omega_err, r_lmn, r_err, freq, m, l, ell, beta, Omega, error, r_min, r_max, m_max, J2, J4, J6, M, R, dim_switch, gyre_data)
	return r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell

def import_r_lmn(filenames, m_max=15):
	# initialize arrays
	Omega_pat = np.zeros((m_max + 1, m_max + 1, len(filenames)))
	r_lmn     = np.zeros((m_max + 1, m_max + 1, len(filenames)))
	Omega_err = np.zeros((m_max + 1, m_max + 1, len(filenames)))
	r_err     = np.zeros((m_max + 1, m_max + 1, len(filenames)))

	extra = 0
	for i in range(len(filenames)):
		data = np.loadtxt(filenames[i])
		ell =   data[:,0]
		emm =   data[:,1]
		Ompat = data[:,2] 
		rlmn  = data[:,3] 
		Omerr = data[:,4]
		rerr  = data[:,5]
		# concatenate onto master array
		# insert into arrays
		for j in range(len(ell)):
			l = int(ell[j])
			m = int(emm[j])
			if(Omega_pat[l][m][i] == 0):
				Omega_pat[l][m][i] = Ompat[j] 
				r_lmn[l][m][i] = rlmn[j]
				Omega_err[l][m][i] = Omerr[j]
				r_err[l][m][i] = rerr[j]
			else:
				# add a new matrix level to Omega_pat and r_lmn
				extra += 1
				Omega_pat = np.reshape(Omega_pat, (m_max + 1, m_max + 1, len(filenames) + extra))
				r_lmn = np.reshape(r_lmn, (m_max + 1, m_max + 1, len(filenames) + extra))
				Omega_pat[l][m][i+extra] = Ompat[j] 
				r_lmn[l][m][i+extra] = rlmn[j]
				Omega_err[l][m][i+extra] = Omerr[j] 
				r_err[l][m][i+extra] = rerr[j]

	# convert to SI units
	for i in range(len(r_lmn)):
		for j in range(len(r_lmn[i])):
			r_lmn[i][j] = r_lmn[i][j] * 1000 # converting from km to m 
			Omega_pat[i][j] = Omega_pat[i][j] * np.pi / (180 * 86400) # converting from deg/day to rad/s
			r_err[i][j] = r_err[i][j] * 1000 # converting from km to m 
			Omega_err[i][j] = Omega_err[i][j] * np.pi / (180 * 86400) # converting from deg/day to rad/s

	N_dim = 1 
	if(len(filenames) > 1):
		N_dim += len(filenames)
	else:
		N_dim += extra
	if(N_dim > 1):
		dim_switch = True 
	else:
		dim_switch = False
	m_max = int(np.amax(ell))

	return r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell, emm

def calculate_mean_uncertainty(filename):
	data = np.loadtxt(filename)
	r_unc = data[:]
	print(np.mean(r_unc))

def calculate_torque_over_surface_mass_density(J2, Omega, R, r, m, l): # unfinished
	D_L = (-(3 - ((9 / 2) * J2 * ((R / r)**2))) * (Omega**2) * (1 + m)) + ((21 / 2) * J2 * ((R / r)**2) * (Omega**2))
	D_V = (-(3 + ((9 / 2) * J2 * ((R / r)**2))) * (Omega**2) * (1 - m)) - ((21 / 2) * J2 * ((R / r)**2) * (Omega**2))

	P = legendre_polynomial # P(1)
	integrand = rho_prime * (r**(l + 2))
	integral = 1 # from 0 to R, dr

	Phi_prime = (G / (r**(l + 1))) * P * ((-1)**(m)) * np.sqrt(4 * np.pi * math.factorial(l - m) / (((2 * l) + 1) * math.factorial(l + m))) * integral  
	dPhi_prime_dtheta = 1

	T_L_over_Sigma =  -(m * (np.pi**2) / D_L) * ((2 * m) + l + 1) * (Phi_prime**2)
	T_V_over_Sigma =   (m * (np.pi**2) / D_V) * (dPhi_prime_dtheta**2)

	print(l, m, T_L_over_Sigma)
	print(l, m, T_V_over_Sigma)

def record_data(planet, l, m, Omegpat, rlmn, Omegaerr, rerr):
	with open(planet + '_frequencies.txt', 'a') as myfile:
		myfile.write(str(l) + '\t')
		myfile.write(str(m) + '\t')
		myfile.write('%.2f' % (Omegpat)) # frequency in degrees per day 
		myfile.write('\t\t\t%.1f' % (rlmn)) # resonance location in km 
		myfile.write('\t\t%.2f' % (Omegaerr)) # error limit frequency in degrees per day 
		myfile.write('\t\t\t%.1f' % (rerr)) # error limit resonance location in km
		myfile.write('\n')

def latex_friendly_output(planet, l, m, Omegpat, rlmn, Omegaerr, rerr):
	with open(planet + '_latex_table_data_' + str(l-m) + '.txt', 'a') as myfile:
		myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(Omegpat), decimals=1)) + ' - ' + str(np.round(np.absolute(float(Omegpat) - float(Omegaerr)), decimals=1)) + ' &' + str(int(np.round(float(rlmn), decimals=0))) + ' + ' + str(int(np.round(np.absolute(float(rlmn) - float(rerr)), decimals=0)))   +  ' \\\\\n')

def latex_friendly_output2(planet, l, m, Omegpat, rlmn, Omegaerr, rerr):
	minOmega = np.amin(Omegaerr)
	maxOmega = np.amax(Omegpat) 
	minr = np.amin(rlmn)
	maxr = np.amax(rerr)

	r_min = 41838 # 6 ring
	r_max = 42572 # 4 ring
	if((r_min < minr < r_max) or (r_min < maxr < r_max) or ((minr < r_min) and (maxr > r_max))):
		if((l - m) % 2 == 0):
			with open(planet + '_latex_table_data_OLR_654.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(minOmega), decimals=1)) + '-' + str(np.round(np.absolute(float(maxOmega)), decimals=1)) + ' &' + str(int(np.round(float(minr), decimals=0))) + '-' + str(int(np.round(np.absolute(float(maxr)), decimals=0)))   +  ' \\\\\n')
		else:
			with open(planet + '_latex_table_data_OVR_654.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(minOmega), decimals=1)) + '-' + str(np.round(np.absolute(float(maxOmega)), decimals=1)) + ' &' + str(int(np.round(float(minr), decimals=0))) + '-' + str(int(np.round(np.absolute(float(maxr)), decimals=0)))   +  ' \\\\\n')
	r_min = 44719 # alpha ring
	r_max = 45661 # beta ring
	if((r_min < minr < r_max) or (r_min < maxr < r_max) or ((minr < r_min) and (maxr > r_max))):
		if((l - m) % 2 == 0):
			with open(planet + '_latex_table_data_OLR_alphabeta.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(minOmega), decimals=1)) + '-' + str(np.round(np.absolute(float(maxOmega)), decimals=1)) + ' &' + str(int(np.round(float(minr), decimals=0))) + '-' + str(int(np.round(np.absolute(float(maxr)), decimals=0)))   +  ' \\\\\n')
		else: 
			with open(planet + '_latex_table_data_OVR_alphabeta.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(minOmega), decimals=1)) + '-' + str(np.round(np.absolute(float(maxOmega)), decimals=1)) + ' &' + str(int(np.round(float(minr), decimals=0))) + '-' + str(int(np.round(np.absolute(float(maxr)), decimals=0)))   +  ' \\\\\n')
def latex_friendly_output3(planet, l, m, Omegpat, rlmn, Omegaerr, rerr):
	r_min = 40900 # Galle inner edge
	r_max = 42900 # Galle outer edge
	if((r_min < rlmn < r_max) or (r_min < rerr < r_max)):
		if((l - m) % 2 == 0):
			with open(planet + '_latex_table_data_OLR.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(Omegaerr), decimals=1)) + '-' + str(np.round(np.absolute(float(Omegpat)), decimals=1)) + ' &' + str(int(np.round(float(rlmn), decimals=0))) + '-' + str(int(np.round(np.absolute(float(rerr)), decimals=0)))   +  ' \\\\\n')
		else:
			with open(planet + '_latex_table_data_OVR.txt', 'a') as myfile:
				myfile.write(str(l) + ' &' + str(m) + ' &' + str(np.round(float(Omegaerr), decimals=1)) + '-' + str(np.round(np.absolute(float(Omegpat)), decimals=1)) + ' &' + str(int(np.round(float(rlmn), decimals=0))) + '-' + str(int(np.round(np.absolute(float(rerr)), decimals=0)))   +  ' \\\\\n')

def include_rings_and_moons(fig, ax, planet, units, r_min, r_max, R, m_min, m_max, c, g_modes, black):
	r_rings, ring_names = prm.import_rings(planet)
	r_moons, moon_names = prm.import_moons(planet)
	r_inner, r_outer, ringbnames = prm.import_broad_rings(planet)
	if(units == 'R_p'):
		for i in range(len(r_rings)):
			r_rings[i] /= R
		for i in range(len(r_moons)):
			r_moons[i] /= R
		for i in range(len(r_inner)):
			r_inner[i] /= R
			r_outer[i] /= R

	#h_offset = 0.0085
	if(g_modes):
		if(planet == 'Uranus'):
			#h_offset = 0.007
			h_offset = 0.015
		else:
			h_offset = 0.02
	else:
		h_offset = 0.02
	print('Adding rings...')
	for i in range(len(r_rings)):
		if(r_min < r_rings[i] < r_max):
			ax.plot([r_rings[i], r_rings[i]], [0, m_max + 1], linestyle='dashed', color=c)
			#ax.text(r_rings[i] - 0.0085, (m_max + m_min) / 2, ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
			if(g_modes == False):
				if(i < 5):
					if(planet == 'Uranus'):
						ax.text(r_rings[i] - h_offset, (m_max + m_min) * 0.75 + (i), ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
					else:
						ax.text(r_rings[i] - h_offset, (m_max + m_min) / 2, ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
				else:
					ax.text(r_rings[i] - h_offset, ((m_max + m_min) / 4) + i, ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
			else:
				if(planet == 'Uranus'):
					ax.text(r_rings[i] - h_offset, (m_max + m_min) * 0.75 + (0.53 * i), ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
				else:
					if(i < 2):
						ax.text(r_rings[i] - h_offset, (m_max + m_min) / 2, ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
					else:
						ax.text(r_rings[i] - h_offset, (m_max + m_min) / 5, ring_names[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)
	
	if(black):
		c1 = ['dimgray', 'mediumblue']
	else:
		c1 = ['lightgray', 'powderblue']

	for i in range(len(r_inner)):
		if((r_min < r_inner[i] < r_max) or (r_min < r_outer[i] < r_max)):
			ax.fill_between([r_inner[i], r_outer[i]], 0, 1, color=c1[i], alpha=0.5, transform=ax.get_xaxis_transform())
			ax.text((r_inner[i] + r_outer[i]) / 2, (m_max + m_min) / 2, ringbnames[i], verticalalignment='bottom', rotation='vertical', fontsize=18, color=c)

	print('Adding moons...')
	for i in range(len(r_moons)):
		if(r_min < r_moons[i] < r_max):
			if((planet == 'Uranus') and (g_modes) and (i > 1)):
				ax.plot([r_moons[i], r_moons[i]], [0, 15], linestyle='dashdot', color='m')
			else:
				ax.plot([r_moons[i], r_moons[i]], [0, m_max + 1], linestyle='dashdot', color='m')
			ax.text(r_moons[i] - h_offset, 1 * (m_max + m_min) / 3, moon_names[i], verticalalignment='bottom', rotation='vertical', fontsize=16, color='m')


	return fig, ax 

def radscan_plot(filenames=[], planet='Uranus', units='R_p', rowspan=10, colspan=10, black=False):
	R, M, Omega, J2, J4, J6 = prm.import_planet(planet)
	r_rings, ring_names = prm.import_rings(planet)
	r_moons, moon_names = prm.import_moons(planet)
	r_radscan, i_over_f_radscan = prm.import_radscan()
	#r_min = r_radscan[0] / R
	#r_max = r_radscan[-1] / R
	r_min = 1.57
	r_max = 1.97
	y_min = 0.125
	y_max = 0.35
	r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell, emm = import_r_lmn(filenames, m_max=25)
	if(units == 'R_p'):
		print('Converting radial units to planetary radii...')
		#a_sr /= R
		for i in range(len(r_rings)):
			r_rings[i] /= R
		for i in range(len(r_moons)):
			r_moons[i] /= R
		for i in range(len(r_lmn)):
			if(dim_switch == False):
				r_lmn[i] /= R 
				r_err[i] /= R
			else:
				for j in range(len(r_lmn[0])):
					r_lmn[i][j] /= R
					r_err[i][j] /= R
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
		cbw = 'w'
	else:
		fig = plt.figure(figsize=(16,8))
		cbw = 'k'
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', 'mean ring normal ' + r'$I/F$' + ' (' + r'$\times 1000$' + ')', r_min, r_max, y_min, y_max)
	ax2 = ax.twiny()
	ax2.set_xlabel(r'$r$ [$10^3$ km]', fontsize=18, color=cbw)
	ax2.set_xlim([r_min * R / 1.0E+06, r_max * R / 1.0E+06])
	ax.tick_params(axis='x', labelsize=18)
	ax.tick_params(axis='y', labelsize=18)
	ax.plot(1000 * r_radscan / R, i_over_f_radscan * 1000 , color='c')
	#ymin, ymax = ax.get_ylim()
	offset = 0.0045
	r_offset = -0.0035
	y_rings = [0.205, 0.21, 0.205, 0.21, 0.195, 0.2, 0.18, 0.2, 0.25, 0.25]
	for i in range(len(r_rings)):
		if(r_min < r_rings[i] < r_max):
			#ax.plot([r_rings[i], r_rings[i]], [0, m_max + 1], linestyle='dashed', color=cbw)
			ax.text(r_rings[i] + r_offset, y_rings[i], ring_names[i], verticalalignment='bottom', fontsize=18, color=cbw)

	for i in range(len(r_moons)):
		if(r_min < r_moons[i] < r_max):
			#ax.plot([r_moons[i], r_moons[i]], [0, m_max + 1], linestyle='dashdot', color='m')
			ax.text(r_moons[i] + r_offset, 0.25, moon_names[i], verticalalignment='bottom', rotation='vertical', fontsize=16, color='m')

	v_factor = 0.01
	for i in range(len(r_lmn)):
		for j in range(len(r_lmn[i])):
			if(r_min < r_lmn[i][j] < r_max):
				if(((i - j) % 2) == 0): # only Lindblad resonances
					if(i == j):
						c = 'r'
					elif(i - j == 2):
						c = 'g'
					elif(i - j == 4):
						if(black):
							c = 'y'
						else:
							c = 'b'
					ax.plot([r_lmn[i][j], r_lmn[i][j]], [y_min, y_max - 0.015 - (v_factor * (i - j))], linestyle='dashdot', color=c)
					ax.text(r_lmn[i][j], 0.335 - (v_factor * (i - j)), str(j), horizontalalignment='center', fontsize=12, color=c)
	
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
		ax2.get_xaxis().tick_top()
		ax2.get_yaxis().tick_left()
		ax.figure.savefig(planet + '_system_radscan.png', dpi=400, facecolor=fig.get_facecolor(), transparent=True)
	else:
		ax.figure.savefig(planet + '_system_radscan.pdf')
	plt.clf()


def make_plot(planet='Uranus', path='/Users/josephahearn/Downloads/work/', filename='summary.h5', filenames=[], gyre_data=True, r_lmn_source='calculate', units='R_p', ell_max=10, m_min=2, r_min=1.5, r_max=2.05, rowspan=10, colspan=10, black=False, write_to_file=True, table_for_latex=False, compare_to_literature=False, add_rad_scan=False, rings_and_moons=True, separate_L_and_v=False, latex_tables_in_text=False, g_modes=False, beta_plots=False, errorbars='one', G=6.6743E-11, dpi=100):
	pathname = path + filename
	R, M, Omega, J2, J4, J6 = prm.import_planet(planet)
	a_sr = (G * M / (Omega**2))**(1/3)
	
	if(add_rad_scan == True):
		r_radscan, i_over_f_radscan = prm.import_radscan()
	
	if(r_lmn_source == 'calculate'):
		print('Calculating resonance locations...')
		r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell = calculate_r_lmn(planet, gyre_data, pathname, path, m_min, r_min, r_max, g_modes, beta_plots)
	elif(r_lmn_source == 'import'):
		write_to_file = False
		r_lmn, r_err, Omega_pat, Omega_err, dim_switch, m_max, ell, emm = import_r_lmn(filenames, m_max=25)

	if(units == 'R_p'):
		print('Converting radial units to planetary radii...')
		a_sr /= R
		for i in range(len(r_lmn)):
			if(dim_switch == False):
				r_lmn[i] /= R 
				r_err[i] /= R
			else:
				for j in range(len(r_lmn[0])):
					r_lmn[i][j] /= R
					r_err[i][j] /= R
	
	if(separate_L_and_v == False):
		print('Initalizing plot...')
		if (black == True):
			fig = plt.figure(facecolor='black', figsize=(16,8))
			cbw = 'w'
		else:
			fig = plt.figure(figsize=(16,8))
			cbw = 'k'
		ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
		#m_max = 11 # just for Saturn comparison plot
		#m_max = 14 # just for zoomed in plot
		fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$m$', r_min, r_max, m_min - 1, m_max + 1)
		#fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$m$', r_min, r_max, m_min - 1, 11) used for Saturn comparison plot
		ax.plot([a_sr, a_sr], [m_min - 1, m_max + 1], linestyle='dashed', color=cbw)
		ax.text(a_sr, m_max / 2, 'synchronous orbit', verticalalignment='top', rotation='vertical', color=cbw)
		ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
	
		if(rings_and_moons == True):
			fig, ax = include_rings_and_moons(fig, ax, planet, units, r_min, r_max, R, m_min, m_max, cbw, g_modes, black)
	
		ax2 = ax.twiny()
		ax2.set_xlabel(r'$r$ [$10^3$ km]', fontsize=18, color=cbw)
		ax2.set_xlim([r_min * R / 1.0E+06, r_max * R / 1.0E+06])
		ax.tick_params(axis='x', labelsize=18)
		ax.tick_params(axis='y', labelsize=18)
	
		if(add_rad_scan == True):
			ax3 = ax.twinx()
			ax3.set_ylabel('mean ring normal ' + r'$I/F$' + ' (' + r'$\times 1000$' + ')', fontsize=18, color='c')
			ax3.plot(1000 * r_radscan / R, i_over_f_radscan * 1000, color='c')
			
		print('Scattering resonant locations...')
		# colors
		if(black == True):
			c = ['r', 'c', 'orange', 'aliceblue', 'y', 'm']
		else:
			c = ['r', 'c', 'orange', 'b', 'y', 'm']
		# markers
		if(len(filenames) == 2): 
			if(planet == 'Saturn'): # used for Saturn comparison plot
				markers = ['s', 'D'] 
				fill = ['full', 'none'] 
			else:
				markers = ['o','o']
				fill = ['full', 'none'] 
		if(len(filenames) == 3): # used for triplet of Aramona & Mankovich models plot
			markers = ['>', 's', '<']
			fill = ['none', 'full', 'none']
		elif(len(filenames) == 4): # used for all Mankovich-related plots
			markers = ['>', 's', '<', 'D']
			fill = ['none', 'full', 'none', 'full']
		elif(len(filenames) == 5): # used for plot of all Uranus models that run well
			markers = ['>', 's', '<', 'D', '*']
			fill = ['none', 'full', 'none', 'full', 'none']
		
		uncertainty = []
		for l in range(len(r_lmn)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1):
				#print(len(uncertainty))
				if(dim_switch == False):
					if(r_min < r_lmn[l][m] < r_max):
						if(errorbars == 'one'):
							jump = 0
						else:
							if(((l - m) % 3) == 0):
								jump = 0
							elif(((l - m) % 3) == 1):
								jump = 0.25
							else:
								jump = 0.5
						#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m] * 180 * 86400 / np.pi) + ' deg/day, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
						#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m]) + ' rad/s, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
						ax.scatter(r_lmn[l][m], m + jump, color=c[l-m])
						#ax.scatter(r_lmn[l][m], m, color=c[l-m], marker=mk.MarkerStyle('s', 'full'))
						#ax.scatter(r_lmn[l][m], m, color=c[l-m], marker=mk.MarkerStyle('*', 'none'))
						if(errorbars == 'all'):
							ax.errorbar(r_lmn[l][m], m + jump, xerr=np.absolute(r_lmn[l][m] - r_err[l][m]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
						elif(errorbars == 'one'):
							uncertainty.append(np.absolute(r_lmn[l][m] - r_err[l][m]))
						#if((l == 2) and (m == 2)):
						#	#asymmetric_error = [error, error**2]
						#	#print(np.shape(asymmetric_error))
						#	ax.errorbar(r_lmn[l][m], m, xerr=np.absolute(r_lmn[l][m] - r_err[l][m]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
						if(write_to_file == True):
							record_data(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000) # frequency in degrees per day, resonance location in km 
						if(table_for_latex == True):
							latex_friendly_output(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
						if(latex_tables_in_text == True):
							latex_friendly_output3(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
				else:
					verticaloffset = 0.1
					N_dim = len(r_lmn[0][0])
					for i in range(N_dim):
						#print(l, m, i)
						if(r_min < r_lmn[l][m][i] < r_max):
							if(((l - m) % 2) == 0):
								jump = 0
							else:
								jump = 0.25
							if(N_dim < 3):
								jump = 0
							#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m] * 180 * 86400 / np.pi) + ' deg/day, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
							#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m]) + ' rad/s, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
							#ax.scatter(r_lmn[l][m][i], m, color=c[l-m], marker=mk.MarkerStyle(markers[i], fillstyle=fill[i]))
							#fc = c[l-m]
							if((i == 1) or (i == 3)):
								fc = c[l-m]
							else:
								fc = 'none'
							if(len(filenames) == 2): 
								if(planet == 'Saturn'): # for Saturn comparison plot
									if(((l == m) or (l - m == 1)) and (m < 11)): 
										if(i == 0):
											fc = c[l-m]
										else:
											fc = 'none'
								else:
									if(i == 0):
										fc = c[l-m]
									else:
										fc = 'none'
							ax.scatter(r_lmn[l][m][i], m + jump + (verticaloffset * i), marker=mk.MarkerStyle(markers[i], fillstyle=fill[i]), facecolors=fc, edgecolors=c[l-m])
							if(errorbars == 'all'):
								ax.errorbar(r_lmn[l][m][i], m + jump + (verticaloffset * i), xerr=np.absolute(r_lmn[l][m][i] - r_err[l][m][i]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
							elif(errorbars == 'one'):
								uncertainty.append(np.absolute(r_lmn[l][m][i] - r_err[l][m][i]))
							#ax.scatter(r_lmn[l][m][i], m, facecolors=c[l-m], edgecolors=c[l-m])
							#ax.scatter(r_lmn[l][m][i], m + (verticaloffset * i), color=c[l-m]) # default
							#print('datum:')
							#print(l, m, i, r_lmn[l][m][i])
							#if((l == 2) and (m == 2) and (i == 0)):
								#ax.errorbar(r_lmn[l][m][i], m, xerr=error**2, color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
								#ax.errorbar(r_lmn[l][m][i], m, xerr=np.absolute(r_lmn[l][m][i] - r_err[l][m][i]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
							
							if(write_to_file == True):
								record_data(planet, l, m, Omega_pat[l][m][i] * 180 * 86400 / np.pi, r_lmn[l][m][i] * R / 1000) # frequency in degrees per day, resonance location in km 
							if(table_for_latex == True):
								latex_friendly_output(planet, l, m, Omega_pat[l][m][i] * 180 * 86400 / np.pi, r_lmn[l][m][i] * R / 1000, Omega_err[l][m][i] * 180 * 86400 / np.pi, r_err[l][m][i] * R / 1000)
							if(latex_tables_in_text == True):
								if(i == 0):
									latex_friendly_output2(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
		fs = 20
		indent = 0
		y_jump = 0 
		if(g_modes):
			indent = 0.74
			y_jump = 6
		if(errorbars != 'all'):
			# uncertainty values are in planet radii
			max_uncertainty = np.amax(uncertainty)
			min_uncertainty = np.amin(uncertainty)
			meanuncertainty = np.mean(uncertainty)
			print('max uncertainty: '  + str(max_uncertainty * R / 1000) + ' km') 
			print('mean uncertainty: ' + str(meanuncertainty * R / 1000) + ' km')
			print('min uncertainty: '  + str(min_uncertainty * R / 1000) + ' km')
			#ax.plot([r_min + 0.01, r_min + 0.01 + max_uncertainty], [(m_min + m_max) / 2, (m_min + m_max) / 2], color=cbw)
			#ax.text(r_min + 0.01, (m_min + m_max) / 2, 'max error', color=cbw)
			ax.plot([r_min + 0.01 + indent, r_min + 0.01 + indent + meanuncertainty], [((m_min + m_max) / 2) - 1 + y_jump, ((m_min + m_max) / 2) - 1 + y_jump], color=cbw)
			ax.text(r_min + 0.02 + indent + meanuncertainty, ((m_min + m_max) / 2) - 1 + y_jump, 'mean error', color=cbw, fontsize=fs, verticalalignment='center')

			# arrow on the right side of the uncertainty bar 
			x_arrow_offset = 0.002
			y_arrow_offset = 0.1
			ax.plot([r_min + 0.01 + meanuncertainty + indent - x_arrow_offset, r_min + 0.01 + meanuncertainty + indent], [((m_min + m_max) / 2) - 1 + y_jump + y_arrow_offset, ((m_min + m_max) / 2) - 1 + y_jump], color=cbw)
			ax.plot([r_min + 0.01 + meanuncertainty + indent - x_arrow_offset, r_min + 0.01 + meanuncertainty + indent], [((m_min + m_max) / 2) - 1 + y_jump - y_arrow_offset, ((m_min + m_max) / 2) - 1 + y_jump], color=cbw)

			#ax.plot([r_min + 0.01, r_min + 0.01 + min_uncertainty], [((m_min + m_max) / 2) - 2, ((m_min + m_max) / 2) - 2], color=cbw)
			#ax.text(r_min + 0.01, ((m_min + m_max) / 2) - 2, 'min error', color=cbw)
		print('Finishing up...')						
		#ax.legend(loc=4, bbox_to_anchor=(0.865, 0.01), fontsize=16)
		if(m_max < 16):
			spacing = 0.6 # for Saturn comparison plot
		else:
			#spacing = 1.2
			spacing = 1.5
		#offset = 0.105 for Saturn comparison plot
		offset = 0.38 #for zoomed out
		if((planet == 'Uranus') and (g_modes)):
			offset = -0.3 # for Uranus g-mode plot
		offset = 0.3
		#offset = 0.01
		nudge = 0.007
		
		if(g_modes == False):
			for i in range(len(c)):
			#for i in range(2): used for Saturn comparison plot
				if(i == 0):
					#ax.text(r_min + 0.01, m_max - 0.5, r'$\ell = m$', color=c[0], fontsize=fs)
					ax.text(r_min + 0.01, m_max, r'$\ell = m$', color=c[0], fontsize=fs)
				else:
					#print('No effect')
					#ax.text(r_min + 0.01, m_max - 0.5 - (i * spacing), r'$\ell - m = $' + str(i), color=c[i], fontsize=fs)
					ax.text(r_min + 0.01, m_max - (i * spacing), r'$\ell - m = $' + str(i), color=c[i], fontsize=fs)
	
		# legend for Saturn comparison plot
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - 0.5,                 marker=mk.MarkerStyle(markers[0], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - 0.5 - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='none'), facecolors='none', edgecolors=cbw)	
		#ax.text(((r_min + r_max) / 2) - offset, m_max - 0.5 - (0 * spacing),       'current rotation rate',    color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - 0.5 - (1 * spacing), '80% of current rotation rate',  color=cbw, verticalalignment='center', fontsize=fs)
	
		# legend for the thin, medium, and thick models
		if(g_modes == True):
			cbw = 'r'
		if(len(filenames) == 2):
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max,                 marker=mk.MarkerStyle(markers[0], fillstyle='full'), facecolors=cbw,    edgecolors=cbw)
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax.text(((r_min + r_max) / 2) - offset, m_max,                 'Lindblad',   color=cbw, verticalalignment='center', fontsize=fs)
			ax.text(((r_min + r_max) / 2) - offset, m_max - spacing,       'corotation', color=cbw, verticalalignment='center', fontsize=fs)
			
		if(len(filenames) > 2):
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max,                 marker=mk.MarkerStyle(markers[0], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='full'), facecolors=cbw,    edgecolors=cbw)
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (2 * spacing), marker=mk.MarkerStyle(markers[2], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax.text(((r_min + r_max) / 2) - offset, m_max,                 'thin',   color=cbw, verticalalignment='center', fontsize=fs)
			ax.text(((r_min + r_max) / 2) - offset, m_max - spacing,       'medium', color=cbw, verticalalignment='center', fontsize=fs)
			ax.text(((r_min + r_max) / 2) - offset, m_max - (2 * spacing), 'thick',  color=cbw, verticalalignment='center', fontsize=fs)
	
		if(len(filenames) > 3):
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (3 * spacing), marker=mk.MarkerStyle(markers[3], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
			ax.text(((r_min + r_max) / 2) - offset, m_max - (3 * spacing), 'shallow',  color=cbw, verticalalignment='center', fontsize=fs)
		if(len(filenames) > 4):
			ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (4 * spacing), marker=mk.MarkerStyle(markers[4], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax.text(((r_min + r_max) / 2) - offset, m_max - (4 * spacing), 'Scheibe et al. (2019)',  color=cbw, verticalalignment='center', fontsize=fs)
		#cbw = 'chartreuse'
		# legend for Mark Marley predictions
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max,                 marker=mk.MarkerStyle(markers[0], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (2 * spacing), marker=mk.MarkerStyle(markers[2], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.text(((r_min + r_max) / 2) - offset, m_max,                 'model 6', color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - spacing,       'model 5', color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - (2 * spacing), 'model 4', color=cbw, verticalalignment='center', fontsize=fs)
	
		#ax.text(r_max - 0.075, m_max, r'error$=\pm$' + str(np.round(error**2, decimals=2)), color='k', fontsize=16)
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
	
			if(add_rad_scan == True):
				ax3.spines['bottom'].set_color('white')
				ax3.spines['left'].set_color('white')
				ax3.spines['top'].set_color('white')
				ax3.spines['right'].set_color('white')
				ax3.yaxis.label.set_color('white')
				ax3.xaxis.label.set_color('white')
				ax3.tick_params(axis='x', colors='white')
				ax3.tick_params(axis='y', colors='white')
				ax3.tick_params(axis='both', direction='in')
		
			ax.figure.savefig(planet + '_system.png', dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
		else:
			#ax.figure.savefig(planet + '_system.png', dpi=dpi)
			ax.figure.savefig(planet + '_system.pdf')
		plt.clf()
	
	elif(separate_L_and_v == True):
		print('Initalizing plot...')
		if (black == True):
			fig = plt.figure(facecolor='black', figsize=(24,8))
			cbw = 'w'
		else:
			fig = plt.figure(figsize=(24,8))
			cbw = 'k'
		ax1 = plt.subplot2grid((12,12),(0, 0), rowspan=rowspan, colspan=int(colspan / 2), fig=fig)
		ax2 = plt.subplot2grid((12,12),(0, 6), rowspan=rowspan, colspan=int(colspan / 2), fig=fig)
		#m_max = 11 # just for Saturn comparison plot
		#m_max = 14 # just for zoomed in plot
		fig, ax1 = pa.title_and_axes(fig, ax1, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$m$', r_min, r_max, m_min - 1, m_max + 1)
		fig, ax2 = pa.title_and_axes(fig, ax2, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$m$', r_min, r_max, m_min - 1, m_max + 1)
		#fig, ax = pa.title_and_axes(fig, ax, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$m$', r_min, r_max, m_min - 1, 11) used for Saturn comparison plot
		if(r_min < a_sr < r_max):
			ax.plot([a_sr, a_sr], [m_min - 1, m_max + 1], linestyle='dashed', color=cbw)
			ax.text(a_sr, m_max / 2, 'synchronous orbit', verticalalignment='top', rotation='vertical', color=cbw)
		ax1.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
		ax2.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
	
		if(rings_and_moons == True):
			fig, ax1 = include_rings_and_moons(fig, ax1, planet, units, r_min, r_max, R, m_min, m_max, cbw, g_modes, black)
			fig, ax2 = include_rings_and_moons(fig, ax2, planet, units, r_min, r_max, R, m_min, m_max, cbw, g_modes, black)
	
		ax3 = ax1.twiny()
		ax3.set_xlabel(r'$r$ [$10^3$ km]', fontsize=18, color=cbw)
		ax3.set_xlim([r_min * R / 1.0E+06, r_max * R / 1.0E+06])
		
		ax4 = ax2.twiny()
		ax4.set_xlabel(r'$r$ [$10^3$ km]', fontsize=18, color=cbw)
		ax4.set_xlim([r_min * R / 1.0E+06, r_max * R / 1.0E+06])
		
		ax1.tick_params(axis='x', labelsize=18)
		ax1.tick_params(axis='y', labelsize=18)
		ax2.tick_params(axis='x', labelsize=18)
		ax2.tick_params(axis='y', labelsize=18)
	
		if(add_rad_scan == True):
			ax5 = ax1.twinx()
			ax5.set_ylabel('mean ring normal ' + r'$I/F$' + ' (' + r'$\times 1000$' + ')', fontsize=18, color='c')
			ax5.plot(1000 * r_radscan / R, i_over_f_radscan * 1000, color='c')
			ax6 = ax2.twinx()
			ax6.set_ylabel('mean ring normal ' + r'$I/F$' + ' (' + r'$\times 1000$' + ')', fontsize=18, color='c')
			ax6.plot(1000 * r_radscan / R, i_over_f_radscan * 1000, color='c')
			
		print('Scattering resonant locations...')
		# colors
		if(black == True):
			c = ['r', 'c', 'orange', 'aliceblue', 'y', 'm']
		else:
			c = ['r', 'c', 'orange', 'b', 'y', 'm']
		# markers
		if(len(filenames) == 2): # used for Saturn comparison plot
			markers = ['s', 'D'] 
			fill = ['full', 'none'] 
		if(len(filenames) == 3): # used for triplet of Aramona & Mankovich models plot
			markers = ['>', 's', '<']
			fill = ['none', 'full', 'none']
		elif(len(filenames) == 4): # used for all Mankovich-related plots
			markers = ['>', 's', '<', 'D']
			fill = ['none', 'full', 'none', 'full']
		elif(len(filenames) == 5): # used for plot of all Uranus models that run well
			markers = ['>', 's', '<', 'D', '*']
			fill = ['none', 'full', 'none', 'full', 'none']
		
		uncertainty = []
		for l in range(len(r_lmn)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1):
				if(dim_switch == False):
					if(r_min < r_lmn[l][m] < r_max):
						if(errorbars == 'one'):
							jump = 0
						else:
							if(((l - m) % 3) == 0):
								jump = 0
							elif(((l - m) % 3) == 1):
								jump = 0.25
							else:
								jump = 0.5
						#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m] * 180 * 86400 / np.pi) + ' deg/day, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
						#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m]) + ' rad/s, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
						
						if(((l - m) % 2) == 0):
							ax1.scatter(r_lmn[l][m], m + jump, color=c[l-m])
						else:
							ax2.scatter(r_lmn[l][m], m + jump, color=c[l-m])
						#ax.scatter(r_lmn[l][m], m, color=c[l-m], marker=mk.MarkerStyle('s', 'full'))
						#ax.scatter(r_lmn[l][m], m, color=c[l-m], marker=mk.MarkerStyle('*', 'none'))

						if(errorbars == 'all'):
							if(((l - m) % 2) == 0):
								ax1.errorbar(r_lmn[l][m], m + jump, xerr=np.absolute(r_lmn[l][m] - r_err[l][m]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
							else:
								ax2.errorbar(r_lmn[l][m], m + jump, xerr=np.absolute(r_lmn[l][m] - r_err[l][m]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
						elif(errorbars == 'one'):
							uncertainty.append(np.absolute(r_lmn[l][m] - r_err[l][m]))
						#if((l == 2) and (m == 2)):
						#	#asymmetric_error = [error, error**2]
						#	#print(np.shape(asymmetric_error))
						#	ax.errorbar(r_lmn[l][m], m, xerr=np.absolute(r_lmn[l][m] - r_err[l][m]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
						if(write_to_file == True):
							record_data(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000) # frequency in degrees per day, resonance location in km 
						if(table_for_latex == True):
							latex_friendly_output(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
						if(latex_tables_in_text == True):
							latex_friendly_output3(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
				else:
					verticaloffset = 0.1
					N_dim = len(r_lmn[0][0])
					for i in range(N_dim):
						#print(l, m, i)
						if(r_min < r_lmn[l][m][i] < r_max):
							jump = 0
							#if(((l - m) % 2) == 0):
							#	jump = 0
							#else:
							#	jump = 0.25
							#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m] * 180 * 86400 / np.pi) + ' deg/day, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
							#print('l=' + str(int(ell[i])) + ', m=' + str(m) + ', Omega_pat = ' + str(Omega_pat[i][m]) + ' rad/s, r = ' + str(r_lmn[i][m] * R / 1000) + ' km')
							#ax.scatter(r_lmn[l][m][i], m, color=c[l-m], marker=mk.MarkerStyle(markers[i], fillstyle=fill[i]))
							#fc = c[l-m]
							if((i == 1) or (i == 3)):
								fc = c[l-m]
							else:
								fc = 'none'
							if(len(filenames) == 2): # for Saturn comparison plot
								if(((l == m) or (l - m == 1)) and (m < 11)): 
									if(i == 0):
										fc = c[l-m]
									else:
										fc = 'none'
							

							if(((l - m) % 2) == 0):
								ax1.scatter(r_lmn[l][m][i], m + jump + (verticaloffset * i), marker=mk.MarkerStyle(markers[i], fillstyle=fill[i]), facecolors=fc, edgecolors=c[l-m])
							else:
								ax2.scatter(r_lmn[l][m][i], m + jump + (verticaloffset * i), marker=mk.MarkerStyle(markers[i], fillstyle=fill[i]), facecolors=fc, edgecolors=c[l-m])
							if(errorbars == 'all'):
								if(((l - m) % 2) == 0):
									ax1.errorbar(r_lmn[l][m][i], m + jump + (verticaloffset * i), xerr=np.absolute(r_lmn[l][m][i] - r_err[l][m][i]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
								else:
									ax2.errorbar(r_lmn[l][m][i], m + jump + (verticaloffset * i), xerr=np.absolute(r_lmn[l][m][i] - r_err[l][m][i]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True)
							elif(errorbars == 'one'):
								uncertainty.append(np.absolute(r_lmn[l][m][i] - r_err[l][m][i]))
							#ax.scatter(r_lmn[l][m][i], m, facecolors=c[l-m], edgecolors=c[l-m])
							#ax.scatter(r_lmn[l][m][i], m + (verticaloffset * i), color=c[l-m]) # default
							#print('datum:')
							#print(l, m, i, r_lmn[l][m][i])
							#if((l == 2) and (m == 2) and (i == 0)):
								#ax.errorbar(r_lmn[l][m][i], m, xerr=error**2, color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
								#ax.errorbar(r_lmn[l][m][i], m, xerr=np.absolute(r_lmn[l][m][i] - r_err[l][m][i]), color=c[l-m], capsize=0.1, capthick=1, xlolims=True) 
							
							if(write_to_file == True):
								record_data(planet, l, m, Omega_pat[l][m][i] * 180 * 86400 / np.pi, r_lmn[l][m][i] * R / 1000) # frequency in degrees per day, resonance location in km 
							if(table_for_latex == True):
								latex_friendly_output(planet, l, m, Omega_pat[l][m][i] * 180 * 86400 / np.pi, r_lmn[l][m][i] * R / 1000, Omega_err[l][m][i] * 180 * 86400 / np.pi, r_err[l][m][i] * R / 1000)
							if(latex_tables_in_text == True):
								if(i == 0):
									latex_friendly_output2(planet, l, m, Omega_pat[l][m] * 180 * 86400 / np.pi, r_lmn[l][m] * R / 1000, Omega_err[l][m] * 180 * 86400 / np.pi, r_err[l][m] * R / 1000)
		fs = 16
		if(errorbars != 'all'):
			max_uncertainty = np.amax(uncertainty)
			min_uncertainty = np.amin(uncertainty)
			meanuncertainty = np.mean(uncertainty)
			print('max uncertainty: '  + str(max_uncertainty * R / 1000) + ' km')
			print('mean uncertainty: ' + str(meanuncertainty * R / 1000) + ' km')
			print('min uncertainty: '  + str(min_uncertainty * R / 1000) + ' km')
			#ax.plot([r_min + 0.01, r_min + 0.01 + max_uncertainty], [(m_min + m_max) / 2, (m_min + m_max) / 2], color=cbw)
			#ax.text(r_min + 0.01, (m_min + m_max) / 2, 'max error', color=cbw)
			ax1.plot([r_min + 0.01, r_min + 0.01 + meanuncertainty], [((m_min + m_max) / 2) - 1, ((m_min + m_max) / 2) - 1], color=cbw)
			ax1.text(r_min + 0.02 + meanuncertainty, ((m_min + m_max) / 2) - 1, 'mean error', color=cbw, fontsize=fs, verticalalignment='center')
			#ax.plot([r_min + 0.01, r_min + 0.01 + min_uncertainty], [((m_min + m_max) / 2) - 2, ((m_min + m_max) / 2) - 2], color=cbw)
			#ax.text(r_min + 0.01, ((m_min + m_max) / 2) - 2, 'min error', color=cbw)
			# arrow on the right side of the uncertainty bar 
			x_arrow_offset = 0.002
			y_arrow_offset = 0.1
			ax1.plot([r_min + 0.01 + meanuncertainty - x_arrow_offset, r_min + 0.01 + meanuncertainty], [((m_min + m_max) / 2) - 1 + y_arrow_offset, ((m_min + m_max) / 2) - 1], color=cbw)
			ax1.plot([r_min + 0.01 + meanuncertainty - x_arrow_offset, r_min + 0.01 + meanuncertainty], [((m_min + m_max) / 2) - 1 - y_arrow_offset, ((m_min + m_max) / 2) - 1], color=cbw)

		print('Finishing up...')						
		#ax.legend(loc=4, bbox_to_anchor=(0.865, 0.01), fontsize=16)
		if(m_max < 16):
			spacing = 0.6 # for Saturn comparison plot
		else:
			spacing = 1.2
		#offset = 0.105 for Saturn comparison plot
		#offset = 0.3 #for zoomed out
		offset = 0.25
		#offset = 0.1
		#offset = 0.01
		#nudge = 0.007
		nudge = 0.01
		
		if(g_modes == False):
			for i in range(len(c)):
			#for i in range(2): used for Saturn comparison plot
				if(i == 0):
					#ax.text(r_min + 0.01, m_max - 0.5, r'$\ell = m$', color=c[0], fontsize=fs)
					ax1.text(r_min + 0.01, m_max, r'$\ell = m$', color=c[0], fontsize=fs)
				else:
					#print('No effect')
					#ax.text(r_min + 0.01, m_max - 0.5 - (i * spacing), r'$\ell - m = $' + str(i), color=c[i], fontsize=fs)
					if((i % 2) == 0):
						ax1.text(r_min + 0.01, m_max - (i * (spacing / 2)), r'$\ell - m = $' + str(i), color=c[i], fontsize=fs)
					else:
						ax2.text(r_min + 0.01, m_max - ((i - 1) * (spacing / 2)), r'$\ell - m = $' + str(i), color=c[i], fontsize=fs)
	
		# legend for Saturn comparison plot
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - 0.5,                 marker=mk.MarkerStyle(markers[0], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - 0.5 - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='none'), facecolors='none', edgecolors=cbw)	
		#ax.text(((r_min + r_max) / 2) - offset, m_max - 0.5 - (0 * spacing),       'current rotation rate',    color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - 0.5 - (1 * spacing), '80% of current rotation rate',  color=cbw, verticalalignment='center', fontsize=fs)
	
		# legend for the thin, medium, and thick models
		if(len(filenames) > 2):
			ax1.scatter(((r_min + r_max) / 2) - offset - nudge, m_max,                 marker=mk.MarkerStyle(markers[0], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax1.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='full'), facecolors=cbw,    edgecolors=cbw)
			ax1.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (2 * spacing), marker=mk.MarkerStyle(markers[2], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax1.text(((r_min + r_max) / 2) - offset, m_max,                 'thin',   color=cbw, verticalalignment='center', fontsize=fs)
			ax1.text(((r_min + r_max) / 2) - offset, m_max - spacing,       'medium', color=cbw, verticalalignment='center', fontsize=fs)
			ax1.text(((r_min + r_max) / 2) - offset, m_max - (2 * spacing), 'thick',  color=cbw, verticalalignment='center', fontsize=fs)
	
		if(len(filenames) > 3):
			ax1.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (3 * spacing), marker=mk.MarkerStyle(markers[3], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
			ax1.text(((r_min + r_max) / 2) - offset, m_max - (3 * spacing), 'shallow',  color=cbw, verticalalignment='center', fontsize=fs)
		if(len(filenames) > 4):
			ax1.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (4 * spacing), marker=mk.MarkerStyle(markers[4], fillstyle='none'), facecolors='none', edgecolors=cbw)
			ax1.text(((r_min + r_max) / 2) - offset, m_max - (4 * spacing), 'adiabatic',  color=cbw, verticalalignment='center', fontsize=fs)
		#cbw = 'chartreuse'
		# legend for Mark Marley predictions
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max,                 marker=mk.MarkerStyle(markers[0], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - spacing,       marker=mk.MarkerStyle(markers[1], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.scatter(((r_min + r_max) / 2) - offset - nudge, m_max - (2 * spacing), marker=mk.MarkerStyle(markers[2], fillstyle='full'), facecolors=cbw, edgecolors=cbw)
		#ax.text(((r_min + r_max) / 2) - offset, m_max,                 'model 6', color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - spacing,       'model 5', color=cbw, verticalalignment='center', fontsize=fs)
		#ax.text(((r_min + r_max) / 2) - offset, m_max - (2 * spacing), 'model 4', color=cbw, verticalalignment='center', fontsize=fs)
	
		#ax.text(r_max - 0.075, m_max, r'error$=\pm$' + str(np.round(error**2, decimals=2)), color='k', fontsize=16)

		ax1.text(r_min + 0.01, 19, 'OLR', fontsize=fs*2)
		ax2.text(r_min + 0.01, 19, 'OVR', fontsize=fs*2)

		if(black == True):
			ax1.spines['bottom'].set_color('white')
			ax1.spines['left'].set_color('white')
			ax1.spines['top'].set_color('white')
			ax1.spines['right'].set_color('white')
			ax1.yaxis.label.set_color('white')
			ax1.xaxis.label.set_color('white')
			ax1.tick_params(axis='x', colors='white')
			ax1.tick_params(axis='y', colors='white')
			ax1.tick_params(axis='both', direction='in')
			ax1.get_xaxis().tick_bottom()
			ax1.get_yaxis().tick_left()
	
			ax2.spines['bottom'].set_color('white')
			ax2.spines['left'].set_color('white')
			ax2.spines['top'].set_color('white')
			ax2.spines['right'].set_color('white')
			ax2.yaxis.label.set_color('white')
			ax2.xaxis.label.set_color('white')
			ax2.tick_params(axis='x', colors='white')
			ax2.tick_params(axis='y', colors='white')
			ax2.tick_params(axis='both', direction='in')
	
			if(add_rad_scan == True):
				ax3.spines['bottom'].set_color('white')
				ax3.spines['left'].set_color('white')
				ax3.spines['top'].set_color('white')
				ax3.spines['right'].set_color('white')
				ax3.yaxis.label.set_color('white')
				ax3.xaxis.label.set_color('white')
				ax3.tick_params(axis='x', colors='white')
				ax3.tick_params(axis='y', colors='white')
				ax3.tick_params(axis='both', direction='in')
		
			ax1.figure.savefig(planet + '_system.png', dpi=dpi, facecolor=fig.get_facecolor(), transparent=True)
		else:
			#ax.figure.savefig(planet + '_system.png', dpi=dpi)
			ax1.figure.savefig(planet + '_system.pdf')
		plt.clf()

	# comparison of resonance locations to what is found in the literature
	if(compare_to_literature == True):
		prm.lit_compare(r_lmn, ell, m_max, r_min, r_max)
	freq_error_plot = False
	if(freq_error_plot == True):
		fig = plt.figure(figsize=(16,8))
		ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
		fig, ax = pa.title_and_axes(fig, ax, '', r'$\Omega_{\rm{pat}}$' + ' [deg/day]', r'$m$', None, None, m_min - 1, m_max + 1)
		ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
		for l in range(len(Omega_pat)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1): 
				if(m == l):
					ax.scatter(Omega_pat[l][m] * 180 * 86400 / np.pi, m, c='k')
					ax.scatter(Omega_err[l][m] * 180 * 86400 / np.pi, m, c='r')
		
		ax.figure.savefig('Omega_pattern_error_estimate.pdf')
		plt.clf()
	resloc_error_plots = False
	if(resloc_error_plots == True):
		fig = plt.figure(figsize=(16,8))
		ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
		fig, ax = pa.title_and_axes(fig, ax, '', r'$r_{\rm{Lind}}$' + ' [1000 km]', r'$m$', None, None, m_min - 1, m_max + 1)
		ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))
		for l in range(len(r_err)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1): 
				if(m == l):
					ax.scatter(r_lmn[l][m] * R / 1.0E+06, m, c='k')
					ax.scatter(r_err[l][m] * R / 1.0E+06, m, c='r')
		ax.figure.savefig('OLR_location_error_estimate1.pdf')
		plt.clf()

		fig = plt.figure(figsize=(16,8))
		ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
		fig, ax = pa.title_and_axes(fig, ax, '', r'Fiducial $r_{\rm{Lind}}$' + ' [1000 km]', r'$\Delta r_{\rm{Lind}}$' + ' [1000 km]')
		for l in range(len(r_err)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1): 
				if(m == l):
					ax.scatter(r_lmn[l][m] * R / 1.0E+06, (r_err[l][m] - r_lmn[l][m]) * R / 1.0E+06, c='b')
		ax.figure.savefig('OLR_location_error_estimate2.pdf')
		plt.clf()

		fig = plt.figure(figsize=(16,8))
		ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
		fig, ax = pa.title_and_axes(fig, ax, '', r'$m$', r'$\Delta r_{\rm{Lind}}$' + ' [1000 km]', m_min - 1, m_max + 1)
		ax.xaxis.set_major_locator(ticker.MaxNLocator(integer=True))
		for l in range(len(r_err)):
			m_min_ell = np.maximum(m_min, l - 5)
			for m in range(m_min_ell, l + 1): 
				if(m == l):
					ax.scatter(m, (r_err[l][m] - r_lmn[l][m]) * R / 1.0E+06, c='b')
		ax.figure.savefig('OLR_location_error_estimate3.pdf')
		plt.clf()

#make_plot(path='/Users/josephahearn/Downloads/work/37_U_sh_up_to_15/', r_min=1.6, r_max=2.7, black=True) # good for Scheibe
#make_plot(path='/Users/josephahearn/Downloads/work/37_U_sh_up_to_15/', r_min=1.0, r_max=2.7, black=True) # good for Scheibe
#make_plot(r_min=1.0, r_max=2.1, black=True) # good for Mankovich
#make_plot(path='/Users/josephahearn/Downloads/work/36_U_CM_up_to_15/', r_min=1.0, r_max=2.1, black=True) # good for Mankovich
#make_plot(path='/Users/josephahearn/Downloads/work/36_U_CM_up_to_15/', r_min=1.0, r_max=2.7, black=True) # good for Mankovich
#make_plot(r_min=1.0, r_max=2.7, black=True) # to include both Scheibe and Mankovich with the same radial limits
#make_plot('Saturn', r_min=1.1, r_max=1.5)
#make_plot('Saturn', r_min=1.0, r_max=1.4, write_to_file=False, compare_to_literature=True)
#make_plot('Neptune', gyre_data=False, r_min=1.0, r_max=2.0, black=False)
#make_plot('Neptune', path='/Users/josephahearn/Downloads/work/41_N1/', gyre_data=True, r_min=1.0, r_max=2.0, black=False)
#make_plot('Uranus', gyre_data=False, r_min=1.0, r_max=3.0)
#make_plot('Uranus', gyre_data=False, r_min=1.0, r_max=2.1, black=True) # to show Mark Marley's predictions
#make_plot('Uranus', path='/Users/josephahearn/Downloads/work/47_U_thick4096/', gyre_data=True, r_min=1.0, r_max=2.1, black=True)
#make_plot('Uranus', filenames=['Uranus_frequencies_thin4096_fast.txt', 'Uranus_frequencies_medium4096_fast.txt', 'Uranus_frequencies_thick4096_fast.txt'], r_lmn_source='import', r_min=1.42, r_max=1.81, black=False)
#make_plot('Uranus', filenames=['Uranus_frequencies_thin4096_fast.txt', 'Uranus_frequencies_medium4096_fast.txt', 'Uranus_frequencies_thick4096_fast.txt', 'Uranus_frequencies_CM_test.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=False)
#make_plot('Uranus', filenames=['Uranus_frequencies_thin4096_fast.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=False, table_for_latex=True)
#make_plot('Uranus', filenames=['Uranus_frequencies_CM_test.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=False, table_for_latex=False)
#make_plot('Uranus', filenames=['Uranus_frequencies_medium4096_fast.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=False, table_for_latex=False)
#make_plot('Uranus', filenames=['Uranus_frequencies_thin4096_slow.txt', 'Uranus_frequencies_medium4096_slow.txt', 'Uranus_frequencies_thick4096_slow.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=True)
#make_plot('Uranus', path='/Users/josephahearn/Downloads/work/46_U_medium4096/', gyre_data=True, r_min=1.0, r_max=2.1, black=True)
#make_plot('Uranus', path='/Users/josephahearn/Downloads/work/45_U_thin4096/', gyre_data=True, r_min=1.0, r_max=2.1, black=True)
#make_plot('Saturn', gyre_data=False, r_min=1.0, r_max=2.2)
#make_plot('Saturn', gyre_data=False, r_min=1.0, r_max=1.4)
#make_plot('Jupiter', gyre_data=False, r_min=1.0, r_max=1.8)
#make_plot(r_min=1.0, r_max=2.7)
#make_plot('Saturn', path='/Users/josephahearn/Downloads/work/16_Saturn_with_BVfreq_removed/', r_min=1.0, r_max=2.0)
#make_plot('Saturn', filenames=['Saturn_frequencies_no_BV_normal_spin_rate.txt', 'Saturn_frequencies_no_BV_slow_spin_rate.txt'], r_lmn_source='import', r_min=1.19, r_max=1.53, write_to_file=False)
#make_plot(r_min=1.0, r_max=2.7)
#make_plot(r_min=1.0, r_max=2.1, path='/Users/josephahearn/Downloads/work/74_Uthin4096_all_modes_found/')
#make_plot('Neptune', r_min=1.0, r_max=2.1, black=True, table_for_latex=True)
#make_plot('Neptune', r_min=1.0, r_max=2.25, write_to_file=True, table_for_latex=True)
#make_plot('Neptune', r_min=1.6, r_max=2.85, m_min=1, write_to_file=True, table_for_latex=True, g_modes=True)
#make_plot(r_min=1.6, r_max=2.75, m_min=1, write_to_file=True, table_for_latex=True, g_modes=True)
#make_plot(r_min=1.0, r_max=2.15, write_to_file=True, table_for_latex=True)

#make_plot(r_min=1.0, r_max=2.1, path='/Users/josephahearn/Downloads/work/80_U3_d_all_modes_found/', write_to_file=True, table_for_latex=True)
#make_plot(r_min=1.6, r_max=2.75, path='/Users/josephahearn/Downloads/work/86_U_thick_including_1/', m_min=1, write_to_file=True, table_for_latex=True, g_modes=True)

#make_plot(r_min=1.0, r_max=2.05, table_for_latex=True)

#make_plot(path='/Users/josephahearn/Downloads/work/77_Umedium4096_all_modes_found/', r_min=1.0, r_max=2.1, table_for_latex=True)
#make_plot('Uranus', filenames=['Uranus_frequencies_thin4096_all_modes.txt', 'Uranus_frequencies_medium4096_all_modes_found.txt', 'Uranus_frequencies_thick4096_all_modes_found.txt', 'Uranus_frequencies_CM_only_nonzero_up_to_15.txt', 'Uranus_frequencies_U3_all_modes_found.txt'], r_lmn_source='import', r_min=1.0, r_max=2.01)
#make_plot(filenames=['Uranus_frequencies_medium4096_all_modes_found.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=True)
#make_plot(filenames=['Uranus_frequencies_U3_all_modes_found.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, black=True)
#make_plot(r_min=1.0, r_max=2.1, black=True)
#make_plot('Neptune', filenames=['Neptune_frequencies_all.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1)
#make_plot('Neptune', path='/Users/josephahearn/Downloads/work/75_N/', r_min=1.5, r_max=2.1, table_for_latex=True)
#radscan_plot(filenames=['Uranus_frequencies_medium_fmodes2.txt'], black=True)
#radscan_plot(filenames=['Uranus_frequencies_medium4096_with_error_bars.txt'])
#radscan_plot(filenames=['Uranus_frequencies_thin4096_with_error_bars.txt'])
#radscan_plot(filenames=['Uranus_frequencies_thick4096_with_error_bars.txt'])
#radscan_plot(filenames=['Uranus_frequencies_CM_with_error_bars_only_up_to_m_15.txt'])
#radscan_plot(filenames=['Uranus_frequencies_U3_with_error_bars.txt'])

#make_plot(path='/Users/josephahearn/Downloads/work/74_Uthin4096_all_modes_found/', r_min=1.0, r_max=2.35, table_for_latex=False, g_modes=True, beta_plots=True)
#make_plot(path='/Users/josephahearn/Downloads/work/77_Umedium4096_all_modes_found/', r_min=1.0, r_max=2.35, table_for_latex=False, g_modes=True, beta_plots=True)
#make_plot(path='/Users/josephahearn/Downloads/work/78_Uthick4096_all_modes_found/', r_min=1.3, r_max=2.35, table_for_latex=False, g_modes=True, beta_plots=True)
#make_plot(path='/Users/josephahearn/Downloads/work/79_U_CM_all_modes_found_weird_after_15/', r_min=1.0, r_max=2.05, table_for_latex=True)
#make_plot(path='/Users/josephahearn/Downloads/work/80_U3_d_all_modes_found/', r_min=1.0, r_max=2.05, table_for_latex=True, g_modes=True)
#make_plot(path='/Users/josephahearn/Downloads/work/84_U_shallow/', r_min=1.0, r_max=2.35, table_for_latex=False, g_modes=True, beta_plots=True)
#make_plot('Neptune', path='/Users/josephahearn/Downloads/work/75_N/', r_min=1.0, r_max=2.05, table_for_latex=True)
#make_plot('Neptune', path='/Users/josephahearn/Downloads/work/75_N/', r_min=1.65, r_max=2.75, write_to_file=True, table_for_latex=True, errorbars='all', g_modes=True, beta_plots=True)

#make_plot(filenames=['Uranus_frequencies_thin4096_with_error_bars.txt', 'Uranus_frequencies_medium4096_with_error_bars.txt', 'Uranus_frequencies_thick4096_with_error_bars.txt', 'Uranus_frequencies_CM_with_error_bars_only_up_to_m_15.txt', 'Uranus_frequencies_U3_with_error_bars.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, write_to_file=False)
#make_plot('Neptune', filenames=['Neptune_frequencies.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, write_to_file=False)
#make_plot(filenames=['Uranus_frequencies_thin4096_with_error_bars.txt'], r_lmn_source='import', r_min=1.0, r_max=2.1, write_to_file=False, separate_L_and_v=False)
#make_plot(filenames=['Uranus_frequencies_thin4096_with_error_bars.txt'], r_lmn_source='import', r_min=1.0, r_max=2.01, write_to_file=False, separate_L_and_v=True)
#make_plot(filenames=['Uranus_frequencies_thin4096_with_error_bars.txt', 'Uranus_frequencies_medium4096_with_error_bars.txt', 'Uranus_frequencies_thick4096_with_error_bars.txt', 'Uranus_frequencies_CM_with_error_bars_only_up_to_m_15.txt', 'Uranus_frequencies_U3_with_error_bars.txt'], r_lmn_source='import', r_min=1.15, r_max=2.01, write_to_file=False, separate_L_and_v=True)
#make_plot('Neptune', filenames=['Neptune_frequencies_final.txt'], r_lmn_source='import', r_min=1.3, r_max=2.25, write_to_file=False, errorbars='all')
#xi_plots(path='/Users/josephahearn/Downloads/work/77_Umedium4096_all_modes_found/')
#calculate_mean_uncertainty('uncertainties_neptune.txt')
#calculate_mean_uncertainty('uncertainties_uranus.txt')

#make_plot(filenames=['Uranus_frequencies_thin1.txt', 'Uranus_frequencies_medium1.txt', 'Uranus_frequencies_thick1.txt', 'Uranus_frequencies_shallow4.txt', 'Uranus_frequencies_adiabatic1.txt'], r_lmn_source='import', r_min=1.15, r_max=2.01, write_to_file=False, separate_L_and_v=True)
#make_plot(filenames=['Uranus_frequencies_thin1.txt', 'Uranus_frequencies_medium1.txt', 'Uranus_frequencies_thick1.txt', 'Uranus_frequencies_shallow4.txt', 'Uranus_frequencies_adiabatic1.txt'], r_lmn_source='import', r_min=1.0, r_max=2.01, write_to_file=False, black=True)
make_plot('Neptune', filenames=['Neptune_frequencies_rev_err_8.txt'], r_lmn_source='import', r_min=1.0, r_max=2.25, errorbars='all', write_to_file=False, black=True)

#make_plot(filenames=['Uranus_frequencies_thin_gmodes2.txt', 'Uranus_frequencies_medium_gmodes2.txt', 'Uranus_frequencies_thick_gmodes2.txt', 'Uranus_frequencies_shallow_gmodes2.txt'], r_lmn_source='import', r_min=1.62, r_max=2.6, m_min=1, write_to_file=False, g_modes=True, black=True)
#make_plot('Neptune', filenames=['Neptune_frequencies_gmodes2.txt', 'Neptune_frequencies_gmodes_cr3.txt'], r_lmn_source='import', r_min=1.49, r_max=2.8, m_min=1, write_to_file=False, errorbars='all', g_modes=True, black=True)


#find_corotation_resonances(filenames=['Neptune_frequencies_gmodes2.txt'], planet='Neptune')