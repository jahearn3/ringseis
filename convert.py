# Convert 
#
# Author: Joseph A'Hearn
# Created 11/17/2020
#
# This program takes the Uranus and Neptune models
#   provided by Nadine Nettelmann and Ravit Helled
#   and outputs files that can be input into GYRE

import numpy as np 
import UN_model_loader as ml 
import useful as uf 
from scipy.optimize import curve_fit
#import compare_models as comp 

def initialize_arrays(n_data):
	return np.zeros(n_data), np.zeros(n_data), np.zeros(n_data), np.zeros(n_data), np.zeros(n_data)

def load_model(filename):
	if(filename[:5] != 'table'):
		m, P, r, T, rho = ml.Nettelmann_model(filename)
	else:
		m, P, r, T, rho = ml.Helled_model(filename)
	return m, P, r, T, rho

def double_gauss(x, a, b, sigma1, sigma2):
	return (a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2)))

def double_gauss_times_one_minus_r(x, a, b, sigma1, sigma2):
	return (1 - x) * ((a * np.exp(-(x)**2/(2*sigma1**2))) + (b * np.exp(-(x)**2/(2*sigma2**2))))

def double_tanh_single_gauss_one_minus_r(x, a, sigma, b, c):
	return (1 - x) * ((b * (1 - np.tanh(25 * (x - 0.8)))) + (c * (1 - np.tanh(15 * (x - 0.2))))) * (a * np.exp(-(x)**2/(2*sigma**2)))

def fit_Net16(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(25 * (x - 0.803049915212053)))) + (c * (1 - np.tanh(15 * (x - 0.1795737879803282))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (2.1123524917657734 * (1 - x))

def fit_Net13s(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(25 * (x - 0.772460027704274)))) + (c * (1 - np.tanh(15 * (x - 0.175965042838437))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.8470085775875662 * (1 - x))

def fit_Net13f(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(25 * (x - 0.757626895535424)))) + (c * (1 - np.tanh(15 * (x - 0.14878409471968446))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.6383646099853342 * (1 - x))

def fit_Sch(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(25 * (x - 0.7490239925459978)))) + (c * (1 - np.tanh(15 * (x - 0.07178373865247327))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (2.471561135633503 * (1 - x))

def fit_Beth(x, a, sigma, b, c):
	return ((b * (1 - np.tanh(25 * (x - 0.7842305293814895)))) + (c * (1 - np.tanh(15 * (x - 0.1813453920113811))))) * (a * np.exp(-(x)**2/(2*sigma**2))) + (1.7356022560861448 * (1 - x))

# previous derivative calculations
#dlnT_dlnP = np.zeros(n_data)
	#dellnP_dellnrho_ad = np.zeros(n_data)
	#dlnT_dlnP_ad = np.zeros(n_data)
	#neg_dellnrho_dellnT_P = np.zeros(n_data)
	
#	iba = 0 
#	for i in range(1, n_data):
#		#dlnp0_dr   = (np.log(P[i]) - np.log(P[i-1])) / (r[i] - r[i-1])
#		#dlnrho0_dr = (np.log(rho[i]) - np.log(rho[i-1])) / (r[i] - r[i-1])
#		if(T[i] == T[i-1]):
#			dlnT_dlnP[i] = dlnT_dlnP[i-1]
#			#dlnT_dlnP[i] = 0
#		elif(P[i] == P[i-1]):
#			dlnT_dlnP[i] = dlnT_dlnP[i-1]
#			#dlnT_dlnP[i] = 1
#		else:
#			dlnT_dlnP[i] = (np.log(T[i]) - np.log(T[i-1])) / (np.log(P[i]) - np.log(P[i-1]))
#		if(np.absolute(dlnT_dlnP[i]) > 0.6):
#			dlnT_dlnP[i] = dlnT_dlnP[i-1]
#		if(dlnT_dlnP[i] == 0):
#			dlnT_dlnP[i] = 1.0E-10
#
#		# gamma
#		if(P[i] == P[i-1]):
#			dellnP_dellnrho_ad[i] = dellnP_dellnrho_ad[i-1]
#			#dellnP_dellnrho_ad[i] = 0
#		elif(rho[i] == rho[i-1]):
#			dellnP_dellnrho_ad[i] = dellnP_dellnrho_ad[i-1]
#			#dellnP_dellnrho_ad[i] = 1
#		else:	
#			#dellnP_dellnrho_ad[i] = (np.log(P[i]) - np.log(P[i-1])) / (np.log(rho[i]) - np.log(rho[i-1]))
#			dellnP_dellnrho_ad[i] = (rho[i] / P[i]) * (P[i] - P[i-1]) / (rho[i] - rho[i-1])
#		if(np.absolute(dellnP_dellnrho_ad[i]) > 5):
#			dellnP_dellnrho_ad[i] = dellnP_dellnrho_ad[i-1]
#		if(dellnP_dellnrho_ad[i] == 0):
#			dellnP_dellnrho_ad[i] = 1.0E-10
#
#		#if(((1 / (dellnP_dellnrho_ad[i])) * dlnp0_dr) < dlnrho0_dr): # I just made up this condition. Need to find out what it should be. 
#		#	N_sq[i] = 1.0E-10
#		#else:
#		#	N_sq[i] = (G * m[i] / ((r[i])**2)) * (((1 / (dellnP_dellnrho_ad[i])) * dlnp0_dr) - dlnrho0_dr) # Equation 3.7 in Unno, p. 11
#		#dlnT_dlnP_ad[i]          = (np.log(T[i]) - np.log(T[i-1])) / (np.log(P[i]) - np.log(P[i-1]))
#		#dlnT_dlnP_ad[i]          = dlnT_dlnP[i]
#		
#		if(rho[i] == rho[i-1]):
#			neg_dellnrho_dellnT_P[i] = neg_dellnrho_dellnT_P[i-1]
#			#neg_dellnrho_dellnT_P[i] = 0
#		elif(T[i] == T[i-1]):
#			neg_dellnrho_dellnT_P[i] = neg_dellnrho_dellnT_P[i-1]
#			#neg_dellnrho_dellnT_P[i] = 1
#		else:
#			neg_dellnrho_dellnT_P[i] = -(np.log(rho[i]) - np.log(rho[i-1])) / (np.log(T[i]) - np.log(T[i-1]))
#		if(neg_dellnrho_dellnT_P[i] < -4):
#			neg_dellnrho_dellnT_P[i] = neg_dellnrho_dellnT_P[i-1]
#		if(neg_dellnrho_dellnT_P[i] == 0):
#			neg_dellnrho_dellnT_P[i] = -1.0E-10
#
#		if((rho[i] < 7) and (rho[i-1] > 7)): # Exclusively in the core is the density greater than 7 g/cm^3
#			icb = i # index of (just outside) the core boundary 
#		if((rho[i] < 1) and (rho[i-1] > 1)) or ((rho[i] < 2) and (rho[i-1] > 2)):
#			iba = i # index of (just outside) the bottom of the atmosphere
#	
#	dlnT_dlnP[0]             = dlnT_dlnP[1]
#	dellnP_dellnrho_ad[0]    = dellnP_dellnrho_ad[1]
#	#N_sq[0]                  = 1.0E-10
#	neg_dellnrho_dellnT_P[0] = neg_dellnrho_dellnT_P[1]
#
#	min_delta = np.absolute(np.amin(neg_dellnrho_dellnT_P))
#	for i in range(n_data):
#		if(neg_dellnrho_dellnT_P[i] < -1.0E-10):
#			neg_dellnrho_dellnT_P[i] /= min_delta # This should make delta be between 0 and -1
#
#
#		if(i < icb): # in the core 
#			dlnT_dlnP[i] = dlnT_dlnP[icb] 
#
#		elif(i == icb): # at the core boundary
#			dellnP_dellnrho_ad[i] = (dellnP_dellnrho_ad[i-1] + dellnP_dellnrho_ad[i+1]) / 2
#
#		if(i == iba): # at the bottom of the atmosphere
#			dellnP_dellnrho_ad[i] = (dellnP_dellnrho_ad[i-1] + dellnP_dellnrho_ad[i+1]) / 2
#
#		if(dellnP_dellnrho_ad[i] > 3):
#			if(r[i] > r[-1] / 2):
#				dellnP_dellnrho_ad[i] = 1.67
#			elif(r[i] < r[-1] / 5):
#				dellnP_dellnrho_ad[i] = 3.0
#			else:
#				dellnP_dellnrho_ad[i] = 2.5
#
#	# median values
#	#gamma_med = np.median(dellnP_dellnrho_ad)
#	#nabla_med = np.median(dlnT_dlnP_ad)
#	#delta_med = np.median(neg_dellnrho_dellnT_P)
#	#print(model)
#	#print('median values:')
#	#print('gamma: ' + str(gamma_med))
#	#print('nabla: ' + str(nabla_med))
#	#print('delta: ' + str(delta_med))
#
#	polynomial = 8 # since these don't have to be monotonic 
#
#	# polynomial fit to smooth out gamma
#	polynomial_coeff = np.polyfit(r, dellnP_dellnrho_ad, polynomial)
#	predict = np.poly1d(polynomial_coeff)
#	gamma_fit = predict(r)
#
#	# polynomial fit to smooth out nabla
#	polynomial_coeff = np.polyfit(r, dlnT_dlnP, polynomial)
#	predict = np.poly1d(polynomial_coeff)
#	nabla_fit = predict(r)
#
#	# polynomial fit to smooth out delta
#	polynomial_coeff = np.polyfit(r, neg_dellnrho_dellnT_P, polynomial)
#	predict = np.poly1d(polynomial_coeff)
#	delta_fit = predict(r)

def log_derivative_with_gutters(derivative, numerator, denominator, min_derivative, max_derivative):
	for i in range(1, len(derivative)):
		if(numerator[i] == numerator[i-1]):
			derivative[i] = derivative[i-1]
		elif(denominator[i] == denominator[i-1]):
			derivative[i] = derivative[i-1]
		else:
			derivative[i] = (np.log(numerator[i]) - np.log(numerator[i-1])) / (np.log(denominator[i]) - np.log(denominator[i-1]))
		if(np.absolute(derivative[i]) > max_derivative):
			derivative[i] = derivative[i-1]
		if(derivative[i] == 0):
			derivative[i] = min_derivative
	derivative[0] = derivative[1]
	return derivative 

def polynomial_fit_to_smooth_out_derivative(derivative, r, polynomial):
	polynomial_coeff = np.polyfit(r, derivative, polynomial)
	predict = np.poly1d(polynomial_coeff)
	derivative_fit = predict(r)
	return derivative_fit


def calculate_derivatives(r, rho, P, T, min_nabla=1.0E-10, max_nabla=0.6, min_gamma=1.0E-10, max_gamma=5.0, min_delta=1.0E-10, max_delta=4, polynomial=8):
	n_data = len(T)
	#N_sq = Brunt-Vaisala frequency squared
	gamma = np.zeros(n_data)
	delta = np.zeros(n_data)
	nabla = np.zeros(n_data)

	nabla = log_derivative_with_gutters(nabla, T, P, min_nabla, max_nabla)
	gamma = log_derivative_with_gutters(gamma, P, rho, min_gamma, max_gamma)
	delta = log_derivative_with_gutters(delta, rho, T, min_delta, max_delta)

	delta = -delta

	# polynomial fit to smooth out gamma
	#gamma_fit = polynomial_fit_to_smooth_out_derivative(gamma, r, polynomial)
	#nabla_fit = polynomial_fit_to_smooth_out_derivative(nabla, r, polynomial)
	#delta_fit = polynomial_fit_to_smooth_out_derivative(delta, r, polynomial)

	#return gamma_fit, delta_fit, nabla_fit, nabla_fit
	return gamma, delta, nabla, nabla

def nudge_method(derivative, tol):
	for i in range(len(derivative) - 2, -1, -1): # starting with the atmosphere and working inward
		if(np.absolute(derivative[i] - derivative[i+1]) > tol):
			if(derivative[i] < derivative[i+1]):
				derivative[i] = derivative[i+1] - tol
			else:
				derivative[i] = derivative[i+1] + tol
	return derivative 

def adiabatic_exponents(r, rho, P, T):
	logrho = np.log(rho)
	logP = np.log(P)
	logT = np.log(T)

	linear_model = np.polyfit(logrho,logP,1)
	#print('average slope for P rho: ' + str(linear_model[0]))
	#comp.quick_plot(logrho, logP)
	linear_model = np.polyfit(logP,logT,1)
	#print('average slope for P T: ' + str(linear_model[0]))
	#comp.quick_plot(logT, logP)
	linear_model = np.polyfit(logT,logrho,1)
	#print('average slope for T rho: ' + str(linear_model[0]))
	#comp.quick_plot(logrho, logT)

	exp_1 = np.zeros(len(r))
	exp_2 = np.zeros(len(r))
	exp_3 = np.zeros(len(r))

	for i in range(1, len(r)):
		if(logrho[i] - logrho[i-1] == 0):
			exp_1[i] = 1
		else:
			exp_1[i] = (logP[i] - logP[i-1]) / (logrho[i] - logrho[i-1])
		if(logT[i] - logT[i-1] == 0):
			exp_2[i] = 1
		else:
			exp_2[i] = (logP[i] - logP[i-1]) / (logT[i] - logT[i-1])
		if(logrho[i] - logrho[i-1] == 0):
			exp_3[i] = 1
		else:
			exp_3[i] = (logP[i] - logP[i-1]) / (logrho[i] - logrho[i-1])

	exp_1[0] = exp_1[1]
	exp_2[0] = exp_2[1]
	exp_3[0] = exp_3[1]

	# fit the exponent to a polynomial
	#polynomial = 8 
	#polynomial_coeff = np.polyfit(r, exp_1, polynomial)
	#predict = np.poly1d(polynomial_coeff)
	#gamma = predict(r)
#
	#polynomial_coeff = np.polyfit(r, exp_2, polynomial)
	#predict = np.poly1d(polynomial_coeff)
	#exponent_2 = predict(r)
	#nabla = (exponent_2 - 1) / exponent_2
#
	#polynomial_coeff = np.polyfit(r, exp_3, polynomial)
	#predict = np.poly1d(polynomial_coeff)
	#exponent_3 = predict(r)
	#delta = -1 / exponent_3

	nabla = np.zeros(len(r))
	delta = np.zeros(len(r))

	gamma = exp_1
	for i in range(len(r)):
		if(exp_2[i] == 0):
			nabla[i] = nabla[i-1]
		else:
			nabla[i] = (exp_2[i] - 1) / exp_2[i] 
		if(exp_3[i] == 0):
			delta[i] = delta[i-1]
		else:
			delta[i] = -1 / exp_3[i] 

	# soften extremes
	soften_extremes = False
	if(soften_extremes == True):
		A = 4
		B = 2 * A
		C = 99
		for i in range(len(r)):
			if(gamma[i] < 1.33):
				gamma[i] = (gamma[i] + (B * 1.33)) / (B + 1)
			elif(gamma[i] > 2.5):
				gamma[i] = (gamma[i] + (A * 2.5)) / (A + 1)
			if(nabla[i] < 0.3):
				nabla[i] = (nabla[i] + (C * 0.3)) / (C + 1)
			elif(nabla[i] > 0.5):
				nabla[i] = (nabla[i] + (C * 0.5)) / (C + 1)
			if(delta[i] < -1.5):
				delta[i] = (delta[i] - (B * 1.5)) / (B + 1)

	nudge = True
	tol = 0.005 # maximum allowable change in the derivative from one data point to the next 
	if(nudge == True):
		gamma = nudge_method(gamma, tol)
		nabla = nudge_method(nabla, tol)
		delta = nudge_method(delta, tol)


	# fit the exponent to a polynomial
	fit_to_polynomials = False
	if(fit_to_polynomials == True):
		polynomial = 4
		polynomial_coeff = np.polyfit(r, gamma, polynomial)
		predict = np.poly1d(polynomial_coeff)
		gamma_fit = predict(r)
	
		polynomial = 2
		polynomial_coeff = np.polyfit(r, nabla, polynomial)
		predict = np.poly1d(polynomial_coeff)
		nabla_fit = predict(r)
	
		polynomial_coeff = np.polyfit(r, delta, polynomial)
		predict = np.poly1d(polynomial_coeff)
		delta_fit = predict(r)
		return gamma_fit, delta_fit, nabla_fit, nabla_fit
	else:
		return gamma, delta, nabla, nabla
	

def calculate_BVfreq(r, rho, P, m, gamma, G=6.67430E-11):
	N_sq = np.zeros(len(r))
	for i in range(1, len(r)):
		if(r[i] - r[i-1] > 0):
			N_sq[i] = (-G*m[i] / r[i]) * (((np.log(rho[i]) - np.log(rho[i-1]))/(np.log(r[i] - np.log(r[i-1])))) - ((1 / gamma[i]) * (((np.log(P[i]) - np.log(P[i-1])))/(np.log(r[i] - np.log(r[i-1]))))))
		else:
			N_sq[i] = N_sq[i-1]
	for i in range(len(r)):
		if(N_sq[i] < 1.0E-10):
			N_sq[i] = 1.0E-10
	return N_sq

def inertia_check(rho_kgm, r_m, M, R):
	#print(np.amax(rho_kgm))
	#print(np.amax(r_m))
	#print(np.amax(M))
	#print(np.amax(R))
	# INERTIA CHECK
	# calculate the constant of inertia based on the density profile
	# Jupiter C = I/(MR^2) = 0.254
	# Saturn  C = I/(MR^2) = 0.210
	# Uranus  C = I/(MR^2) = 0.225
	integrand = (r_m**4) * rho_kgm
	#print(np.amax(r_m))
	#print(np.amax(rho_kgm))
	#comp.quick_plot(r_m, integrand)
	area = np.trapz(integrand, x=r_m)
	#print(area)
	I = 4 * np.pi * area
	#print(I)
	C = (2 /3) * I / (M * (R**2))
	print('Uranus constant of inertia: 0.225')
	print('Umodel constant of inertia: ' + str(C))

	#rho_mean = np.mean(rho_kgm)
	#integrand = r_m**4
	#area = np.trapz(integrand, x=r_m)
	#I = 4 * np.pi * rho_mean * area
	#print(I)
	#C = I / (M * (R**2))
	#print('Uranus constant of inertia: 0.225')
	#print('Umodel constant of inertia: ' + str(C))
#
	#I = 4 * np.pi * rho_mean * (np.amax(r_m)**5) / 5
	#print(I)
	#C = I / (M * (R**2))
	#print('Uranus constant of inertia: 0.225')
	#print('Umodel constant of inertia: ' + str(C))
#
	#print(M)
	#print(R)
	#I_sphere = 0.4 * M * (R**2)
	#print(I / I_sphere)
	#print(I_sphere / I)

	return C

def fit_density_profile_with_moment_of_inertia(r_m, rho, R, M, model, tol=0.1, C_p=0.225):
	error = 1.0 # to ensure at least one density fit
	n_iter = 0 # to count how many iterations it takes to converge
	while(error > tol):
		# density fitting function must be monotonic

		#parameters, covariance = curve_fit(double_tanh_single_gauss_one_minus_r, r_m / R, rho)  
		#rho_r = double_tanh_single_gauss_one_minus_r(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
		if(n_iter == 0): # running through the fitting function with constants unique to each model  #filenames = [, , , , ]
			if(model == 'U_TBL_Ne16.dat'):
				parameters, covariance = curve_fit(fit_Net16, r_m / R, rho)  
				rho_r = fit_Net16(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(model == 'U1_Ne13.dat'):
				parameters, covariance = curve_fit(fit_Net13s, r_m / R, rho)  
				rho_r = fit_Net13s(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(model == 'U2_Ne13.dat'):
				parameters, covariance = curve_fit(fit_Net13f, r_m / R, rho)  
				rho_r = fit_Net13f(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(model == 'U3_Sh19.dat'):
				parameters, covariance = curve_fit(fit_Sch, r_m / R, rho)  
				rho_r = fit_Sch(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
			elif(model == 'Uicy_Be17.dat'):
				parameters, covariance = curve_fit(fit_Beth, r_m / R, rho)  
				rho_r = fit_Beth(r_m / R, parameters[0], parameters[1], parameters[2], parameters[3])
			for i in range(len(r_m)):
				if(rho_r[i] < 5.0E-04): # to ensure that no density is negative (or less than the "edge" of the atmosphere should be)
					rho_r[i] = 5.0E-04 # corresponding to a density that might be near the "edge" of the atmosphere, based on some interior models
		else: # modify the density profile without running the fitting function again
			break
			#if(C > C_p): # if C > C_p, that means the density needs to become more centralized
			#	rho_r += 0.1 * (1 - (r_m / R))
			#else: # if C < C_p, that means the density needs to become less centralized, but we don't want rho to become non-monotonic
			#	rho_r += 0.1 * (r_m / R)


		
		C = inertia_check(rho_r * 1000, r_m, M, R)
		error = np.absolute(C - C_p) / C_p
		n_iter += 1
		if(n_iter > 15):
			print('did not converge after ' + str(n_iter) + ' iterations.')
			break 
	#if(n_iter <= 15):
		#print(str(n_iter) + ' iterations to converge.')
	#comp.quick_plot(r_m / R, rho_r)

	return rho_r

def calculate_mass_by_finite_difference(n_data, r_m, rho_r):
	m_r = np.zeros(n_data)
	m_r[0] = rho_r[0] * 1000 * (4/3) * np.pi * (r_m[0]**3)
	interior_mass = m_r[0]
	for i in range(1, n_data):
		shell_mass = rho_r[i]* 1000 * (4/3) * np.pi * (r_m[i]**3 - r_m[i-1]**3)
		m_r[i] = interior_mass + shell_mass
		interior_mass += shell_mass
	return m_r, interior_mass

def modify_model(m, P, r, T, rho):
	n_data = len(m)
	M_Earth = 5.97237E+27 	# Mass of the Earth in g (Luzum et al. 2009) Bibcode:  2011CeMDA.110..293L
	R_Earth_cm = 6.3781366E+08 # Radius of the Earth in cm (Selected Astronomical Constants, 2011, https://web.archive.org/web/20130826043456/http://asa.usno.navy.mil/SecK/2011/Astronomical_Constants_2011.txt)
	R_Earth = 6.3781366E+06 # Radius of the Earth in m (Selected Astronomical Constants, 2011, https://web.archive.org/web/20130826043456/http://asa.usno.navy.mil/SecK/2011/Astronomical_Constants_2011.txt)
	P_conversion = 1.0E+10 # conversion factor from GPa to dyne/cm^2
	G = 6.67430E-11 # gravitational constant in m^3 kg^-1 s^-2 using the 2018 CODATA-recommended value

	M = float(m[-1]) * M_Earth * 1.0E-03 # conversion to kg
	R = r[-1] * R_Earth # conversion to m 
	r_m = r * R_Earth

	# calculations 
	N_sq = np.zeros(n_data)
	#N_sq, nabla, gamma_ad, nabla_ad, delta = initialize_arrays(n_data)
	
	# MODIFICATION OF DENSITY AND MASS ========================================================================================
	# fit to smooth out density
	# modify density profile until the moment of inertia is close enough
	rho_r = fit_density_profile_with_moment_of_inertia(r_m, rho, R, M, model)
	

	# alter the density profile if the constant of inertia is too small or too large
	# left unfinished

	# MASS
	# mass continuity and conservation of mass
	
	# attempt with integrals
	#m_r = np.zeros(n_data) # interior mass
	#interior_mass = 0
	#for i in range(n_data): 
	#	integrand = 4.0 * np.pi * (r_m[:i]**2) * rho_kgm[:i]
	#	m_r[i] = np.trapz(integrand, x=r_m[:i])
	#	interior_mass += m_r[i]
	#print('integrated mass after first density adjustment: ' + str(interior_mass) + ' kg')
	#print('percent error: ' + str(100 * np.absolute(interior_mass - M) / M) + '%')

	# conservation of mass
	#print('model mass: ' + str(M) + ' kg')
	# integrate rho to calcuate total mass
	#mass = np.absolute(np.trapz(rho_r, x=r))
	#print('integrated mass after first density adjustment: ' + str(mass * M_Earth * 1.0E-03) + ' kg')
	#print('percent error: ' + str(100 * np.absolute(mass - m[-1]) / m[-1]) + '%')

	# finite difference method
	m_r, interior_mass = calculate_mass_by_finite_difference(n_data, r_m, rho_r)
	
	#comp.quick_plot(r_m / R, m_r)

	# multiplying density function by a factor to reduce mass error
	#rho_fit *= np.absolute(m_tot / mass)
	#rho_kgm *= np.absolute(M / m_r)
	#rho_r *= np.absolute(M / m_r)

	#comp.quick_plot(r, rho_r)
	
	# recalculate mass after latest density adjustement
	#mass = np.absolute(np.trapz(rho_r, x=r))
	#interior_mass = 0
	#for i in range(n_data): 
	#	integrand = 4.0 * np.pi * (r_m[:i]**2) * rho_kgm[:i]
	#	m_r[i] = np.trapz(integrand, x=r_m[:i])
	#	interior_mass += m_r[i]
	#print('integrated mass after second density adjustment: ' + str(mass * M_Earth * 1.0E-03) + ' kg')
	#print('integrated mass after second density adjustment: ' + str(interior_mass) + ' kg')
	#print('percent error: ' + str(100 * np.absolute(mass - m[-1]) / m[-1]) + '%')
	#print('percent error: ' + str(100 * np.absolute(interior_mass - M) / M) + '%')

	# recalculate m, P, T according to new density profile
	# m
	#m_fit = np.absolute(rho_fit * (4 / 3) * np.pi * (r**3))
	#print('total mass after first recalculation: ' + str(m_fit[-1] * M_Earth * 1.0E-03))
	#m_fit *= np.absolute(mass / m_fit[-1])
	#print('total mass after second recalculation: ' + str(m_fit[-1] * M_Earth * 1.0E-03))

	# P
	#P_fit = rho_fit * other_terms # this needs to be done using hydrostatic equilibrium
	# dP/dr = -(G*M(r)/r^2)*rho
	P_r = np.zeros(n_data)
	for i in range(n_data - 1, -1, -1): # starting from the edge of the atmosphere and moving inward 
		integrand = (m_r[i:] / (r_m[i:]**2)) * rho_r[i:] * 1000
		P_r[i] = G * np.trapz(integrand, x=r_m[i:])
	
	P_r[-1] = 1.0E+04

	print('central pressure: ' + str(P_r[0] / P_conversion) + ' TPa') # should be in the range of 50-70 TPa = 50,000-70,000 GPa
	#comp.quick_plot(r_m / R, P_r)
	# T
	#T_fit = P_fit * other_terms # not a priority

	return m_r, P_r, r, T_r, rho_r

def find_boundaries(r, rho, R):
	for j in range(1, len(r)):
		if((rho[j-1] / rho[j]) > 1.5):
			if((r[j] / R) < 0.5): # core boundary
				cb = ((r[j] / R) + (r[j-1] / R)) / 2
				idxcb = j
			else: # atmosphere boundary
				ab = ((r[j] / R) + (r[j-1] / R)) / 2
				idxab = j
	return idxab, idxcb 

def insert_double_point(idxab, idxcb, m, P, r, T, rho):
	m = np.insert(m, idxab, m[idxab])
	P = np.insert(P, idxab, P[idxab])
	r = np.insert(r, idxab, r[idxab])
	T = np.insert(T, idxab, T[idxab])
	rho = np.insert(rho, idxab, rho[idxab+1])
	m = np.insert(m, idxcb, m[idxcb])
	P = np.insert(P, idxcb, P[idxcb])
	r = np.insert(r, idxcb, r[idxcb])
	T = np.insert(T, idxcb, T[idxcb])
	rho = np.insert(rho, idxcb, rho[idxcb+1])
	return m, P, r, T, rho

def remove_superfluous_points(): # this function is unfinished 
	r_prev = 0
	j = 0
	for i in range(n_data):
		if(i > 0) and (0.0 < (r[i] - r[i-1]) < 0.002): # skipping indices where radii are too close to each other but not equal to avoid spurious mesh points
			if(r_prev == 0):
				r_prev = r[i-1] # the previous r that was written to the .gyre file
				continue
			elif(0.0 < (r[i] - r_prev) < 0.002):
				r_prev = r[i-1]
				continue
		j +=1
	return

def check_boundaries(m, P, r, T, rho):
	idxab, idxcb = find_boundaries(r, rho, r[-1])
	#print('atmosphere boundary from the outside in: ')
	#print(idxab+1)
	#print(idxab)
	#print(idxab-1)
	#print('core boundary from the outside in: ')
	#print(idxcb+1)
	#print(idxcb)
	#print(idxcb-1)
	m, P, r, T, rho = insert_double_point(idxab, idxcb, m, P, r, T, rho)

	return m, P, r, T, rho

def gyremodelsetup(m, P, r, T, rho, model, modify=False):
	R_Earth_m = 6.3781366E+06
	M_Earth = 5.97237E+24 	# Mass of the Earth in kg (Luzum et al. 2009) Bibcode:  2011CeMDA.110..293L
	if(modify == True):
		m_r, P_r, r_m, T_r, rho_r = modify_model(m, P, r, T, rho)
	else:
		m_r, P_r, r_m, T_r, rho_r = check_boundaries(m, P, r, T, rho)
	M = np.amax(m_r)
	R = np.amax(r_m)
	C = inertia_check(rho_r * 1000, r_m * R_Earth_m, M * M_Earth, R * R_Earth_m)
	# calculate derivatives
	gamma, delta, nabla, nabla_ad = adiabatic_exponents(r_m, rho_r, P_r, T_r)
	#gamma, delta, nabla, nabla_ad = calculate_derivatives(r_m, rho_r, P_r, T_r)

	# make sure units are correct
	#print('r: ' + str(np.amax(r) * R_Earth_cm))
	#print('m_r: ' + str(np.amax(m_r) * 1000))
	#print('P_r: ' + str(np.amax(P_r) * P_conversion / 1.0E+09))
	#N_sq = 0
	calculate_BV = False
	if(calculate_BV == True):
		N_sq = calculate_BVfreq(r_m, rho_r, P_r, m_r, gamma)
	else:
		N_sq = np.zeros(len(r_m)) 
	#for i in range(len(r_m)):
	#	N_sq[i] = 1.0E-10
	write_to_file(model, len(r_m), m_r, r_m, P_r, T_r, rho_r, nabla, N_sq, gamma, nabla_ad, delta)



def write_to_file(model, n_data, m_r, r, P_r, T, rho_r, nabla, N_sq, gamma, nabla_ad, delta):
	R_Earth_cm = 6.3781366E+08
	M_Earth = 5.97237E+27 	# Mass of the Earth in g (Luzum et al. 2009) Bibcode:  2011CeMDA.110..293L
	P_conversion = 1.0E+10
	R = r[-1] 
	#print('first and last mass points: ')
	#print(m_r[0], m_r[-1])

	# write to file
	with open(model[0:2] + 'model.gyre', 'a') as myfile:
		#					number of grid points			Mass 									radius 						luminosity					version number * 100
		#myfile.write('  ' + str(n_data) + '   ' + str("{:.9e}".format(m[0] * M_Earth)) + '   ' + str(r[0] * R_Earth) + '   ' + str(1.00000000e+00) + '   ' + str(101) + '\n')
		myfile.write('  ' + str(n_data) + '   ' + str("{:.9e}".format(m_r[-1] * M_Earth)) + '   ' + str("{:.9e}".format(r[-1] * R_Earth_cm)) + '   1.00000000e+00   101\n')
		# 19 columns
		# 1 grid point index
		# 2 r (cm)
		# 3 interior mass (g)
		# 4 luminosity (erg/s)
		# 5 total pressure (dyn/cm^2)
		# 6 temperature (K)
		# 7 density (g/cm^3)
		# 8 d ln P / d ln rho 
		# 9 Brunt-Vaisala frequency squared (s^-2)
		#10 (del ln P / del ln rho)_ad
		#11 (d ln T / d ln P)_ad
		#12 -(del ln rho / del ln T)_P
		#13 opacity kappa (cm^2/g)
		#14 kappa(del ln kappa / del ln T)_rho (cm^2/g)
		#15 kappa(del ln kappa / del ln rho)_T (cm^2/g)
		#16 nuclear energy generation/loss rate (erg/(s g))
		#17 epsilon_nuc (del ln epsilon_nuc / del ln T)_rho (erg/(s g))
		#18 epsilon_nuc (del ln epsilon_nuc / del ln rho)_T (erg/(s g))
		#19 rotational angular velocity (rad/s)

		for i in range(n_data):
			if(i < 9):
				myfile.write('    ' + str(i + 1) + '  ')
			elif(i < 99):
				myfile.write('   ' + str(i + 1) + '  ')
			elif(i < 999):
				myfile.write('  ' + str(i + 1) + '  ')
			elif(i < 9999):
				myfile.write(' ' + str(i + 1) + '  ')
			elif(i < 99999):
				myfile.write(str(i + 1) + '  ')
			#myfile.write(str("{:.9e}".format(r[i] * R_Earth)) + '  ') # radius 
			myfile.write(str("{:.9e}".format(r[i] * R_Earth_cm)) + '  ') # radius 

			#myfile.write(str("{:.9e}".format(m_fit[i] * M_Earth)) + '  ') # interior mass 
			#myfile.write(str("{:.9e}".format(m[i] * M_Earth)) + '  ') # interior mass 
			myfile.write(str("{:.9e}".format(m_r[i] * M_Earth)) + '  ') # interior mass
			#myfile.write(str("{:.9e}".format(m_r[i] * 1000)) + '  ') # interior mass

			myfile.write('1.00000000e+00   ') # luminosity
			#myfile.write(str("{:.9e}".format(P[i] * P_conversion)) + '  ') # total pressure
			myfile.write(str("{:.9e}".format(P_r[i] * P_conversion)) + '  ') # total pressure
			#myfile.write(str("{:.9e}".format(P_r[i] * P_conversion / 1.0E+09)) + '  ') # total pressure

			myfile.write(str("{:.9e}".format(T[i])) + '  ') # temperature

			#myfile.write(str("{:.9e}".format(rho[i])) + '  ') # density
			#myfile.write(str("{:.9e}".format(rho_fit[i])) + '  ') # density
			myfile.write(str("{:.9e}".format(rho_r[i])) + '  ') # density

			#myfile.write(str("{:.9e}".format(dlnT_dlnP[i])) + '  ') # d ln T / d ln P 
			#myfile.write(str("{:.9e}".format(nabla_med)) + '  ') # d ln T / d ln P 
			#myfile.write('0.30000000e+00   ') # d ln T / d ln P 
			#myfile.write(str("{:.9e}".format(nabla_fit[i])) + '  ') # d ln T / d ln P 
			myfile.write(str("{:.9e}".format(nabla[i])) + '  ') # d ln T / d ln P 

			myfile.write('1.00000000e-10   ') # Brunt-Vaisala frequency squared

			#myfile.write(str("{:.9e}".format(dellnP_dellnrho_ad[i])) + '  ') # (del ln P / del ln rho)_ad
			#myfile.write(str("{:.9e}".format(gamma_med)) + '  ') # (del ln P / del ln rho)_ad
			#myfile.write('3.00000000e+00   ') # (del ln P / del ln rho)_ad
			#myfile.write('1.66666667e+00   ') # (del ln P / del ln rho)_ad (non-relativistic)
			#myfile.write(str("{:.9e}".format(gamma_fit[i])) + '  ')
			myfile.write(str("{:.9e}".format(gamma[i])) + '  ')

			#myfile.write(str("{:.9e}".format(dlnT_dlnP_ad[i])) + '  ') # (d ln T / d ln P)_ad
			#myfile.write(str("{:.9e}".format(nabla_med)) + '  ') # (d ln T / d ln P)_ad
			#myfile.write('0.30000000e+00   ') # (d ln T / d ln P)_ad
			myfile.write(str("{:.9e}".format(nabla_ad[i])) + '  ') # (d ln T / d ln P)_ad

			#myfile.write(str("{:.9e}".format(neg_dellnrho_dellnT_P[i])) + '  ') # -(del ln rho / del ln T)_P
			myfile.write(str("{:.9e}".format(delta[i])) + '  ') # -(del ln rho / del ln T)_P
			#myfile.write('0.00000000e+00   ') # delta (-(del ln rho / del ln T)_P)

			myfile.write('1.00000000e+00   ') # opacity
			myfile.write('0.00000000e+00   ') # kappa(del ln kappa / del ln T)_rho (cm^2/g)
			myfile.write('0.00000000e+00   ') # kappa(del ln kappa / del ln rho)_T (cm^2/g)
			myfile.write('1.00000000e+00   ') # nuclear energy generation/loss rate (erg/(s g))
			myfile.write('0.00000000e+00   ') # epsilon_nuc (del ln epsilon_nuc / del ln T)_rho (erg/(s g))
			myfile.write('0.00000000e+00   ') # epsilon_nuc (del ln epsilon_nuc / del ln rho)_T (erg/(s g))
			#myfile.write('1.01500000e-04 \n') # rotational angular velocity (rad/s) for solid body rotation of Uranus (value corresponds to 17.2 hr)
			myfile.write('0.00000000e+00 \n') # rotational angular velocity (rad/s)

def generate_homogeneous_sphere_model(n_data=10000, M=8.6810E+25, R=2.5559E+07, R_mean=2.5362E+07):
	rho_mean = (3 * M / (4 * np.pi * (R_mean**3))) / 1000 # g/cm^3
	r = np.linspace(0, R, n_data)
	rho = np.zeros(n_data)
	for i in range(n_data):
		rho[i] = rho_mean
	C = inertia_check(rho * 1000, r, M, R)


	#P = np.zeros(n_data)
	#T = np.zeros(n_data)
	#m = calculate_mass_by_finite_difference(n_data, r, rho)


def main():
	#filenames = ['N1_Ne13.dat', 'N2b_Ne13.dat', 'N3_Sh19.dat', 'U_TBL_Ne16.dat', 'U1_Ne13.dat', 'U2_Ne13.dat', 'U3_Sh19.dat', 'Uicy_Be17.dat']
	#filenames = ['U_TBL_Ne16.dat', 'U1_Ne13.dat', 'U2_Ne13.dat', 'U3_Sh19.dat', 'Uicy_Be17.dat']
	#filenames = ['Uicy_Be17.dat']
	#filenames = ['U3_Sh19.dat'] # This is the only Uranus model that worked in GYRE on the first attempt
	#filenames = ['U1_Ne13.dat']
	#filenames = ['U_TBL_Ne16.dat']
	filenames = ['N3_Sh19.dat']
	for model in filenames:
		print(model)
		m, P, r, T, rho = load_model(model)
		# flip models around because MESA format goes from the core outward
		m = np.flip(m, 0)
		P = np.flip(P, 0)
		r = np.flip(r, 0)
		T = np.flip(T, 0)
		rho = np.flip(rho, 0)

		#delta_r = np.zeros(len(r))
		#R = np.amax(r)
		#for i in range(1, len(r)):
		#	delta_r[i] = (r[i] - r[i-1])
		#comp.quick_plot(r/R, delta_r)

		gyremodelsetup(m, P, r, T, rho, model)


#main()
#generate_homogeneous_sphere_model()