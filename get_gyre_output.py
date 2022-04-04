# Get GYRE Output
#
# Author: Joseph A'Hearn
# Created 02/11/2021
#
# This program gets relevant output from the GYRE summary file 
#

import numpy as np 

def mr(path='/Users/josephahearn/Downloads/work/summary.h5'):
	with open(path) as f:
		content = f.readlines()
	content = [x.strip() for x in content] 
	M = float(content[3][:25]) / 1000 # gives a mass in kg
	R = float(content[3][25:]) / 100  # gives a radius in m
	return M, R

def M():
	return mr()[0]

def R():
	return mr()[1]

def output_data(path='/Users/josephahearn/Downloads/work/summary.h5'):
	data = np.loadtxt(path, skiprows=6)
	k =       data[:,0]
	l =       data[:,1] 						 
	m =       data[:,2] 
	n_pg =    data[:,3] 
	n_p =     data[:,4] 
	n_g =     data[:,5] 
	omega_r = data[:,6] 
	omega_i = data[:,7]
	beta =    data[:,8] / (4 * np.pi)
	E =       data[:,9]
	E_norm =  data[:,10]
	k = k.astype('int')
	l = l.astype('int')
	m = m.astype('int')
	n_pg = n_pg.astype('int')
	n_p = n_p.astype('int')
	n_g = n_g.astype('int')
	return k, l, m, n_pg, n_p, n_g, freq_r, freq_i, beta, E, E_norm

def relevant_data(path='/Users/josephahearn/Downloads/work/summary.h5'):
	data = np.loadtxt(path, skiprows=6)
	j =      data[:,0]
	l =      data[:,1]
	n_pg =   data[:,3] 
	omega =  data[:,6] 
	beta =   data[:,8] / (4 * np.pi)
	j = j.astype('int')
	l = l.astype('int')
	nonzero_indices = []
	for i in range(len(n_pg)):
		if(n_pg[i] != 0):
			nonzero_indices.append(i)
	j =     np.delete(j,     nonzero_indices)
	l =     np.delete(l,     nonzero_indices)
	omega = np.delete(omega, nonzero_indices)
	beta =  np.delete(beta,  nonzero_indices)
	return j, l, omega, beta

def relevant_data2(path='/Users/josephahearn/Downloads/work/summary.h5'): # Used for the shallow model
	data = np.loadtxt(path, skiprows=6)
	j =      data[:,0]
	l =      data[:,1]
	n_pg =   data[:,3] 
	n_p =    data[:,4] 
	n_g =    data[:,5] 
	omega =  data[:,6] 
	beta =   data[:,8] / (4 * np.pi)
	j = j.astype('int')
	l = l.astype('int')
	n_pg = n_pg.astype('int')
	n_p = n_p.astype('int')
	n_g = n_g.astype('int')
	nonzero_indices = []
	for i in range(len(n_pg)):
		if(l[i] < 16):
			if(n_pg[i] != 0):
				nonzero_indices.append(i)
		else:
			if((n_pg[i] == 1) and (n_g[i] == 0) and (n_p[i] == 1)):
				pass
			else:
				nonzero_indices.append(i)

	j =     np.delete(j,     nonzero_indices)
	l =     np.delete(l,     nonzero_indices)
	omega = np.delete(omega, nonzero_indices)
	beta =  np.delete(beta,  nonzero_indices)
	return j, l, omega, beta

def data_for_gmodes(path='/Users/josephahearn/Downloads/work/summary.h5'):
	data = np.loadtxt(path, skiprows=6)
	j =      data[:,0]
	l =      data[:,1]
	n_pg =   data[:,3] 
	n_p =    data[:,4] 
	n_g =    data[:,5] 
	omega =  data[:,6] 
	beta =   data[:,8] / (4 * np.pi)
	j = j.astype('int')
	l = l.astype('int')
	n_pg = n_pg.astype('int')
	n_p = n_p.astype('int')
	n_g = n_g.astype('int')
	nonzero_indices = []
	for i in range(len(n_pg)):
		#if((n_pg[i] == 1) and (n_g[i] == 0)):
		if((n_pg[i] == -1) and (n_g[i] == 1) and (n_p[i] == 0)):
			pass 
		else:
			nonzero_indices.append(i)
	j =     np.delete(j,     nonzero_indices)
	l =     np.delete(l,     nonzero_indices)
	omega = np.delete(omega, nonzero_indices)
	beta =  np.delete(beta,  nonzero_indices)
	return j, l, omega, beta

def data_for_gmodes2(path='/Users/josephahearn/Downloads/work/summary.h5'):
	data = np.loadtxt(path, skiprows=6)
	j =      data[:,0]
	l =      data[:,1]
	n_pg =   data[:,3] 
	n_p =    data[:,4] 
	n_g =    data[:,5] 
	omega =  data[:,6] 
	beta =   data[:,8] / (4 * np.pi)
	j = j.astype('int')
	l = l.astype('int')
	n_pg = n_pg.astype('int')
	n_p = n_p.astype('int')
	n_g = n_g.astype('int')
	nonzero_indices = []
	for i in range(len(n_pg)):
		if(l[i] < 16):
			if((n_pg[i] == -1) and (n_g[i] == 1) and (n_p[i] == 0)):
				pass 
			else:
				nonzero_indices.append(i)
		else:
			if(n_pg[i] != 0):
				nonzero_indices.append(i)

	j =     np.delete(j,     nonzero_indices)
	l =     np.delete(l,     nonzero_indices)
	omega = np.delete(omega, nonzero_indices)
	beta =  np.delete(beta,  nonzero_indices)
	return j, l, omega, beta


def eigs_data(j, l, n=0, path='/Users/josephahearn/Downloads/work/eigs/j'):
	if(path == '/Users/josephahearn/Downloads/work/84_U_shallow/eigs/j'):
		if(l > 15): 
		#	n = 1 # Uncomment for an f-mode search
			n = 0 # Uncomment for a g-mode search
	print('eigs data function running...')
	if(n >= 0):
		sign = '+'
	else:
		sign = '-'
	path += str(j) + '_l' + str(l) + '_n' + sign + str(abs(n)) + '.h5'
	data = np.loadtxt(path, skiprows=6)
	lambda_re    = data[:,0]
	lambda_im    = data[:,1] 
	dE_dx        = data[:,2] 
	x            = data[:,3] # fractional radius
	V_2          = data[:,4] 
	As           = data[:,5] 
	U            = data[:,6] 
	c_1          = data[:,7]
	Gamma_1      = data[:,8] 
	Omega_rot    = data[:,9]
	y_1_re       = data[:,10]
	y_1_im       = data[:,11]
	y_2_re       = data[:,12] 
	y_2_im       = data[:,13] 
	y_3_re       = data[:,14] 
	y_3_im       = data[:,15] 
	y_4_re       = data[:,16] 
	y_4_im       = data[:,17] 
	xi_r_re      = data[:,18]
	xi_r_im      = data[:,19] 
	xi_h_re      = data[:,20]
	xi_h_im      = data[:,21]
	eul_phi_re   = data[:,22]
	eul_phi_im   = data[:,23] 
	deul_phi_re  = data[:,24] 
	deul_phi_im  = data[:,25] 
	eul_P_re     = data[:,26] 
	eul_P_im     = data[:,27] 
	eul_rho_re   = data[:,28] 
	eul_rho_im   = data[:,29]
	prop_type_re = data[:,30] 
	prop_type_im = data[:,31]
	M_r          = data[:,32]
	P            = data[:,33]
	rho          = data[:,34]

	return lambda_re

def xi_data(j, l, n=0, path='/Users/josephahearn/Downloads/work/eigs/j', g_modes=False):
	if((path == '/Users/josephahearn/Downloads/work/84_U_shallow/eigs/j') or (path == '/Users/josephahearn/Downloads/work/92_shallow_including_1/eigs/j')):
		if(l > 15): 
			if(g_modes == False):
				n = 1 
			else:
				n = 0 
	if(n >= 0):
		sign = '+'
	else:
		sign = '-'
	path += str(j) + '_l' + str(l) + '_n' + sign + str(abs(n)) + '.h5'
	data = np.loadtxt(path, skiprows=6)
	x       = data[:,3]
	xi_r_re = data[:,18]
	xi_r_im = data[:,19] 
	xi_h_re = data[:,20]
	xi_h_im = data[:,21]
	eul_rho = data[:,28] 
	rho     = data[:,34]

	with open(path) as f:
		content = f.readlines()
	content = [x.strip() for x in content]
	beta = float(content[3][306:330]) / (4 * np.pi)
	R = float(content[3][356:])

	x /= x[-1]

	return x, xi_r_re, xi_h_re, rho, beta, R, eul_rho

def metadata(j, l, n=0, path='/Users/josephahearn/Downloads/work/eigs/j'):
	if(n >= 0):
		sign = '+'
	else:
		sign = '-'
		path += str(j) + '_l' + str(l) + '_n' + sign + str(abs(n)) + '.h5'
	with open(path) as f:
		content = f.readlines()
	content = [x.strip() for x in content] 
	#print(content[3])
	n = int(content[3][:4])
	j = int(content[3][28]) 
	l = int(content[3][53])
	l_i_re = float(content[3][56:79])
	l_i_im = float(content[3][81:105])
	m = int(content[3][128])
	n_p = int(content[3][153])
	n_g = int(content[3][178])
	n_pg = int(content[3][203])
	omega_re = float(content[3][206:229])
	omega_im = float(content[3][231:254])
	E = float(content[3][256:280])
	E_norm = float(content[3][281:304])
	beta = float(content[3][306:330])
	M = float(content[3][331:355])
	R = float(content[3][356:])

	return n, j, l, l_i_re, l_i_im, m, n_p, n_g, n_pg, omega_re, omega_im, E, E_norm, beta, M, R 
	
	
	

#j, l, freq, beta = relevant_data()
#for i in range(len(l)):
#	print(l[i], freq[i], beta[i])
#eigs_data(2, 2, 0)
#xi_data(2, 2, 0)
#metadata(2, 2, 0)