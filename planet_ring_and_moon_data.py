# Planet, ring, and moon data
#
# Author: Joseph A'Hearn
# Created 11/03/2021
#
# This program stores data on planets, rings, and moons

import numpy as np 

def import_planet(planet):
	if(planet == 'Jupiter'):
		R = 7.1492E+07
		M = 1.8982E+27
		Omega = 1.7585E-04 # 9.925 hr
		J2 =  0.01469562 
		J4 = -0.00059131 
		J6 =  0.00002078 
		# J2, J4, J6 from
		# Jacobson, R.A.: 2013, JUP310 orbit solution
	elif(planet == 'Saturn'):
		R = 6.0268E+07
		M = 5.683E+26
		#Omega = 1.65269E-04 # rad/s 10.56056 hr from Mankovich et al 2019
		Omega = 1.65269E-04 * 0.8 # for plot comparing resonance locations at different rotation rates
		J2 = 0.01629071
		J4 = -0.00093583
		J6 = 0.00008614
		# J2, J4, J6 from
		# Jacobson, R. A., Antreasian, P. G., Bordi, J. J., Criddle, K. E., Ionasescu,R., Jones, J. B., Mackenzie, R. A., Pelletier, F. J., Owen Jr., W. M., Roth, D. C., and Stauch, J. R.: 
		# 2006 ``The gravity field of the Saturnian system from satellite observations and spacecraft tracking data'', Astronomical Journal 132, 6.
	elif(planet == 'Uranus'):
		R = 2.5559E+07
		M = 8.6810E+25
		Omega = 1.0124E-04 # 17.24 hr (slow rate); I previously had 1.015E-04 listed here
		#Omega = 1.0527E-04 # 16.58 hr (fast rate) Helled et al. 2010
		J2 = 0.0035107 # default
		#J2 = 0.00351069 # A&M thick
		#J2 = 0.00351068 # A&M medium and thin
		J4 = -0.0000342 # default
		#J4 = -0.0000351949 # A&M thick
		#J4 = -0.0000341705 # A&M medium
		#J4 = -0.0000330171 # A&M thin
		J6 = 0.0
		# default J2, J4, J6 from
		# Jacobson, R. A.: 2014 ``The Orbits of the Uranian Satellites and Rings, the Gravity Field of the Uranian System, and the Orientation of the Pole of Uranus'',Astronomical Journal 148,76-88.
	elif(planet == 'Neptune'):
		R = 2.4764E+07 
		M = 1.02413E+26
		#Omega = 1.0834E-04 # 16.11 hr (fast rate)
		Omega = 9.996E-05 # 17.46 hr (slow rate) Helled et al. 2010, Mankovich model
		#J2 = 0.0034091 # Brozovic et al. 2020
		J2 = 0.00340843 # Jacobson 2009, Mankovich model
		J4 = -0.00003340 # Jacobson 2009, Mankovich model
		J6 = 0.0
		# J2, J4, J6 from
		# Jacobson, R. A.: 2009 ``The Orbits of the Neptunian Satellites and the Orientation of the Pole of Neptune'', Astronomical Journal 137,4322.
	return R, M, Omega, J2, J4, J6 

def import_rings(planet):
	if(planet == 'Jupiter'):
		r_halo = 1.0725E+08
		r_main = 1.2575E+08
		r_Amgoss = 1.555E+08
		r_Thgoss = 1.775E+08
		r_rings = [r_halo, r_main, r_Amgoss, r_Thgoss]
		ring_names = ['halo', 'main', 'Amalthea gossamer', 'Thebe gossamer']
	elif(planet == 'Saturn'):
		#D_inner_edge = 6.69E+07
		#D_outer_edge = 7.451E+07
		C_inner_edge = 7.4658E+07
		B_C_boundary = 9.2000E+07
		B_outer_edge = 1.1758E+08
		A_inner_edge = 1.2217E+08
		A_outer_edge = 1.36775E+08
		F_ring = 1.4018E+08
		r_rings = [C_inner_edge, B_C_boundary, B_outer_edge, A_inner_edge, A_outer_edge, F_ring]
		ring_names = ['C inner edge', 'B/C boundary', 'B outer edge', 'A inner edge', 'A outer edge', 'F']
	elif(planet == 'Uranus'): # from Nicholson 2018
		r_6 = 4.18376E+07
		r_5 = 4.22353E+07
		r_4 = 4.25715E+07
		r_alpha = 4.47190E+07
		r_beta = 4.56614E+07
		r_eta = 4.71763E+07
		r_gamma = 4.76268E+07
		r_delta = 4.83006E+07
		r_lambda = 5.00249E+07
		r_epsilon = 5.11496E+07
		r_rings = [r_6, r_5, r_4, r_alpha, r_beta, r_eta, r_gamma, r_delta, r_lambda, r_epsilon]
		ring_names = ['6', '5', '4', r'$\alpha$', r'$\beta$', r'$\eta$', r'$\gamma$', r'$\delta$', r'$\lambda$', r'$\epsilon$']
	elif(planet == 'Neptune'):
		r_Le_Verrier = 5.32E+07
		r_Arago = 5.72E+07
		r_Adams = 6.2932E+07
		r_rings = [r_Le_Verrier, r_Arago, r_Adams]
		ring_names = ['Le Verrier ring', 'Arago ring', 'Adams ring']
	return r_rings, ring_names

def import_broad_rings(planet):
	if(planet == 'Uranus'):
		r_zeta_inner_Voyager = 3.7E+07
		r_zeta_inner_Keck = 3.785E+07
		r_zeta_outer_Voyager = 3.95E+07
		r_zeta_outer_Keck = 4.135E+07
		r_inner_edge = [r_zeta_inner_Voyager, r_zeta_inner_Keck]
		r_outer_edge = [r_zeta_outer_Voyager, r_zeta_outer_Keck] 
		ring_names = [r'$\zeta$ (Voyager)', r'$\zeta$ (Keck)'] # 1986U2R/
	elif(planet == 'Neptune'):
		r_inner_Galle = 4.09E+07
		r_outer_Galle = 4.29E+07
		r_inner_Lassell = 5.32E+07
		r_outer_Lassell = 5.72E+07
		r_inner_edge = [r_inner_Galle, r_inner_Lassell]
		r_outer_edge = [r_outer_Galle, r_outer_Lassell]
		ring_names = ['Galle ring', 'Lassell ring']
	return r_inner_edge, r_outer_edge, ring_names

def import_moons(planet):
	if(planet == 'Uranus'): # from Jacobson 1998
		r_Cordelia  = 4.9751722E+07
		r_Ophelia   = 5.3763390E+07
		r_Bianca    = 5.9165550E+07
		r_Cressida  = 6.1766730E+07
		r_Desdemona = 6.2658364E+07
		r_Juliet    = 6.4358222E+07
		r_Portia    = 6.6097265E+07
		r_Rosalind  = 6.9926795E+07
		r_moons = [r_Cordelia, r_Ophelia, r_Bianca, r_Cressida, r_Desdemona, r_Juliet, r_Portia, r_Rosalind]
		moon_names = ['Cordelia', 'Ophelia', 'Bianca', 'Cressida', 'Desdemona', 'Juliet', 'Portia', 'Rosalind']
	elif(planet == 'Neptune'):
		r_Naiad    = 4.8224E+07
		r_Thalassa = 5.0074E+07
		r_Despina  = 5.2526E+07
		r_Galatea  = 6.1953E+07
		r_Larissa  = 7.3548E+07
		r_moons = [r_Naiad, r_Thalassa, r_Despina, r_Galatea, r_Larissa]
		moon_names = ['Naiad', 'Thalassa', 'Despina', 'Galatea', 'Larissa']
	elif(planet == 'Saturn'):
		r_Pan        = 1.33584E+08
		r_Daphnis    = 1.36505E+08
		r_Atlas      = 1.37670E+08
		r_Prometheus = 1.39380E+08
		r_Pandora    = 1.41720E+08
		r_moons = [r_Pan, r_Daphnis, r_Atlas, r_Prometheus, r_Pandora]
		moon_names = ['Pan', 'Daphnis', 'Atlas', 'Prometheus', 'Pandora']
	elif(planet == 'Jupiter'):
		r_Metis = 1.28852E+08
		r_Adrastea = 1.290E+08
		r_Amalthea = 1.81366E+08
		r_Thebe = 2.22452E+08
		r_moons = [r_Metis, r_Adrastea, r_Amalthea, r_Thebe]
		moon_names = ['Metis', 'Adrastea', 'Amalthea', 'Thebe']
	else:
		r_moons = []
		moon_names = []
	return r_moons, moon_names 

def import_radscan():
	data = np.loadtxt('C2685219_radscan_042621.tab', skiprows=16)
	r = data[:,0]
	i_over_f = data[:,1]
	return r, i_over_f

def lit_compare(r_lmn, ell, m_max, r_min, r_max, rowspan=10, colspan=10, black=False):
	if (black == True):
		fig = plt.figure(facecolor='black', figsize=(16,8))
	else:
		fig = plt.figure(figsize=(16,8))
	ax = plt.subplot2grid((10,10),(0, 0), rowspan=rowspan, colspan=colspan)
	#fig, ax = pa.title_and_axes(fig, ax, 'Discrepancies', r'$m$', 'Percent', m_min - 1, 9)
	fig, ax = pa.title_and_axes(fig, ax, 'Discrepancies', r'$\ell$', 'Percent', 1)
	Omega_pat_Marley = np.zeros((m_max + 1, m_max + 1)) # the first two rows (or columns?) will be zeros for m=0,1 
	Omega_pat_Mankov = np.zeros((m_max + 1, m_max + 1)) # (16, 14) to include the predictions for largest l and largest m 

	Omega_pat_Marley[2][2] = convert_to_pattern_frequency(286.1) # These numbers are from Marley 2014
	Omega_pat_Marley[2][1] = convert_to_pattern_frequency(150.1)
	Omega_pat_Marley[3][3] = convert_to_pattern_frequency(296.7)
	Omega_pat_Marley[3][2] = convert_to_pattern_frequency(211.6)
	Omega_pat_Marley[3][1] = convert_to_pattern_frequency(118.2)
	Omega_pat_Marley[4][4] = convert_to_pattern_frequency(310.7)
	Omega_pat_Marley[4][3] = convert_to_pattern_frequency(249.3)
	Omega_pat_Marley[4][2] = convert_to_pattern_frequency(183.9)
	Omega_pat_Marley[4][1] = convert_to_pattern_frequency(100.1)
	Omega_pat_Marley[5][5] = convert_to_pattern_frequency(324.2)
	Omega_pat_Marley[5][4] = convert_to_pattern_frequency(276.3)
	Omega_pat_Marley[5][3] = convert_to_pattern_frequency(223.3)
	Omega_pat_Marley[5][2] = convert_to_pattern_frequency(163.9)
	Omega_pat_Marley[6][6] = convert_to_pattern_frequency(335.6)
	Omega_pat_Marley[6][5] = convert_to_pattern_frequency(297.5)
	Omega_pat_Marley[6][4] = convert_to_pattern_frequency(254.2)
	Omega_pat_Marley[6][3] = convert_to_pattern_frequency(206.8)
	Omega_pat_Marley[6][2] = convert_to_pattern_frequency(150.0)
	Omega_pat_Marley[7][7] = convert_to_pattern_frequency(348.6)
	Omega_pat_Marley[7][6] = convert_to_pattern_frequency(313.5)
	Omega_pat_Marley[7][5] = convert_to_pattern_frequency(277.8)
	Omega_pat_Marley[7][4] = convert_to_pattern_frequency(238.0)
	Omega_pat_Marley[7][3] = convert_to_pattern_frequency(192.9)
	Omega_pat_Marley[8][8] = convert_to_pattern_frequency(355.0)
	Omega_pat_Marley[8][7] = convert_to_pattern_frequency(327.3)
	Omega_pat_Marley[8][6] = convert_to_pattern_frequency(297.5)
	Omega_pat_Marley[8][5] = convert_to_pattern_frequency(263.8)
	Omega_pat_Marley[8][4] = convert_to_pattern_frequency(225.7)

	Omega_pat_Mankov[2][2]   = convert_from_deg_day_to_rad_s(np.mean([1769.2, 2169.3])) # These first couple correspond to observed detections
	Omega_pat_Mankov[3][3]   = convert_from_deg_day_to_rad_s(np.mean([1730.3, 1736.7]))

	Omega_pat_Mankov[4][4]   = convert_from_deg_day_to_rad_s(np.mean([1657.87, 1673.41])) # These numbers are from Mankovich et al 2019
	Omega_pat_Mankov[5][5]   = convert_from_deg_day_to_rad_s(np.mean([1592.08, 1596.05])) # The first segment corresponds to model predictions whose waves are observed in the rings (Table 1)
	Omega_pat_Mankov[6][6]   = convert_from_deg_day_to_rad_s(np.mean([1537.10, 1539.51]))
	Omega_pat_Mankov[7][7]   = convert_from_deg_day_to_rad_s(np.mean([1491.73, 1493.72]))
	Omega_pat_Mankov[9][7]   = convert_from_deg_day_to_rad_s(np.mean([1655.86, 1657.35]))
	Omega_pat_Mankov[8][8]   = convert_from_deg_day_to_rad_s(np.mean([1453.93, 1455.23]))
	Omega_pat_Mankov[9][9]   = convert_from_deg_day_to_rad_s(np.mean([1421.83, 1422.55]))
	Omega_pat_Mankov[13][9]  = convert_from_deg_day_to_rad_s(np.mean([1626.48, 1627.46]))
	Omega_pat_Mankov[10][10] = convert_from_deg_day_to_rad_s(np.mean([1394.03, 1394.71]))
	Omega_pat_Mankov[13][11] = convert_from_deg_day_to_rad_s(np.mean([1451.53, 1453.07]))
	Omega_pat_Mankov[5][4]   = convert_from_deg_day_to_rad_s(np.mean([1871.22, 1875.42]))
	Omega_pat_Mankov[10][7]  = convert_from_deg_day_to_rad_s(np.mean([1723.99, 1725.28]))
	Omega_pat_Mankov[11][8]  = convert_from_deg_day_to_rad_s(np.mean([1644.89, 1645.81]))
	Omega_pat_Mankov[14][9]  = convert_from_deg_day_to_rad_s(np.mean([1667.72, 1668.85]))

	Omega_pat_Mankov[11][11] = convert_from_deg_day_to_rad_s(np.mean([1368.5, 1371.5])) # This second segment corresponds to model predictions without associated wave detections (Table 2)
	Omega_pat_Mankov[12][12] = convert_from_deg_day_to_rad_s(np.mean([1346.9, 1349.7]))
	Omega_pat_Mankov[13][13] = convert_from_deg_day_to_rad_s(np.mean([1327.7, 1330.1]))
	Omega_pat_Mankov[8][6]   = convert_from_deg_day_to_rad_s(np.mean([1742.1, 1747.6])) 
	Omega_pat_Mankov[10][8]  = convert_from_deg_day_to_rad_s(np.mean([1586.9, 1591.0]))
	Omega_pat_Mankov[11][9]  = convert_from_deg_day_to_rad_s(np.mean([1532.9, 1536.4]))
	Omega_pat_Mankov[12][10] = convert_from_deg_day_to_rad_s(np.mean([1488.4, 1491.6]))
	Omega_pat_Mankov[14][12] = convert_from_deg_day_to_rad_s(np.mean([1419.3, 1421.8]))
	Omega_pat_Mankov[12][8]  = convert_from_deg_day_to_rad_s(np.mean([1695.6, 1699.5]))
	Omega_pat_Mankov[14][10] = convert_from_deg_day_to_rad_s(np.mean([1568.6, 1571.6]))
	Omega_pat_Mankov[15][11] = convert_from_deg_day_to_rad_s(np.mean([1521.4, 1524.1]))
	Omega_pat_Mankov[6][5]   = convert_from_deg_day_to_rad_s(np.mean([1737.0, 1743.6]))
	Omega_pat_Mankov[7][6]   = convert_from_deg_day_to_rad_s(np.mean([1646.4, 1652.1]))
	Omega_pat_Mankov[8][7]   = convert_from_deg_day_to_rad_s(np.mean([1578.0, 1582.8]))
	Omega_pat_Mankov[9][8]   = convert_from_deg_day_to_rad_s(np.mean([1523.9, 1528.1]))
	Omega_pat_Mankov[10][9]  = convert_from_deg_day_to_rad_s(np.mean([1479.8, 1483.5]))
	Omega_pat_Mankov[11][10] = convert_from_deg_day_to_rad_s(np.mean([1443.0, 1446.3]))
	Omega_pat_Mankov[12][9]  = convert_from_deg_day_to_rad_s(np.mean([1581.1, 1584.6]))
	Omega_pat_Mankov[13][10] = convert_from_deg_day_to_rad_s(np.mean([1530.1, 1533.1]))
	Omega_pat_Mankov[15][10] = convert_from_deg_day_to_rad_s(np.mean([1604.7, 1607.6]))

	Omega_pat_Mankov[5][3]   = convert_from_deg_day_to_rad_s(np.mean([2314.3, 2324.3])) # This third segment corresponds to more model predictions (Table 3)
	Omega_pat_Mankov[6][4]   = convert_from_deg_day_to_rad_s(np.mean([2034.3, 2042.2]))
	Omega_pat_Mankov[7][5]   = convert_from_deg_day_to_rad_s(np.mean([1861.1, 1867.6]))
	Omega_pat_Mankov[9][5]   = convert_from_deg_day_to_rad_s(np.mean([2063.1, 2069.4]))
	Omega_pat_Mankov[10][6]  = convert_from_deg_day_to_rad_s(np.mean([1901.5, 1906.8]))
	Omega_pat_Mankov[11][7]  = convert_from_deg_day_to_rad_s(np.mean([1784.6, 1789.0]))
	Omega_pat_Mankov[2][1]   = convert_from_deg_day_to_rad_s(np.mean([3359.2, 3400.7]))
	Omega_pat_Mankov[3][2]   = convert_from_deg_day_to_rad_s(np.mean([2423.7, 2433.4]))
	Omega_pat_Mankov[4][3]   = convert_from_deg_day_to_rad_s(np.mean([2061.8, 2070.8]))
	Omega_pat_Mankov[7][4]   = convert_from_deg_day_to_rad_s(np.mean([2177.3, 2185.2]))
	Omega_pat_Mankov[8][5]   = convert_from_deg_day_to_rad_s(np.mean([1968.1, 1974.5]))
	Omega_pat_Mankov[9][6]   = convert_from_deg_day_to_rad_s(np.mean([1826.1, 1831.4]))
	Omega_pat_Mankov[11][6]  = convert_from_deg_day_to_rad_s(np.mean([1970.5, 1975.6]))
	Omega_pat_Mankov[12][7]  = convert_from_deg_day_to_rad_s(np.mean([1841.5, 1845.9]))
	Omega_pat_Mankov[13][8]  = convert_from_deg_day_to_rad_s(np.mean([1743.8, 1747.5]))

	for o in range(5): # o for order
		for l in range(len(ell)):
			m = l - o 
			if(r_min < r_lmn[l][m] < r_max):
				if(Omega_pat_Marley[l][m] != 0):
					pct_error_Marley = 100 * np.absolute(Omega_pat[l][m] - Omega_pat_Marley[l][m]) / Omega_pat_Marley[l][m]
					ax.scatter(l, pct_error_Marley, color=c[o], marker='s')
				if(Omega_pat_Mankov[l][m] != 0):
					if(l < 11):
						pct_error_Mankov = 100 * np.absolute(Omega_pat[l][m] - Omega_pat_Mankov[l][m]) / Omega_pat_Mankov[l][m]
						ax.scatter(l, pct_error_Mankov, color=c[o], marker='v')
	ymin, ymax = ax.get_ylim()
	spacing = (ymax - ymin) / 20
	ax.text(9.5, ymax - (1 * spacing), r'$\ell = m$'    , color='r'   , fontsize=16)
	ax.text(9.5, ymax - (2 * spacing), r'$\ell - m = 1$', color='c'   , fontsize=16)
	ax.text(9.5, ymax - (3 * spacing), r'$\ell - m = 2$', color='g'   , fontsize=16)
	ax.text(9.5, ymax - (4 * spacing), r'$\ell - m = 3$', color='m'   , fontsize=16)
	ax.text(9.5, ymax - (5 * spacing), r'$\ell - m = 4$', color='b'   , fontsize=16)
	xmarker = 7.35
	xtext = 7.5
	yjump = spacing / 10
	ax.scatter(xmarker, ymax - (1 * spacing) + yjump, color='k', marker='s')
	ax.text(xtext, ymax - (1 * spacing), 'Marley (2014)', color='k', fontsize=16)		
	ax.scatter(xmarker, ymax - (2 * spacing) + yjump, color='k', marker='v')
	ax.text(xtext, ymax - (2 * spacing), 'Mankovich et al. (2019)', color='k', fontsize=16)
	pa.save_and_clear_plot(fig, ax, filename=planet + '_system_discrepancies', black=black)

def marley88():
	r_lmn = np.zeros((6, 6, 3))
	Omega_pat = np.zeros((6, 6, 3))
	dim_switch = True
	m_max = 15
	ell = [0, 1, 2, 3, 4, 5]
	# Marley et al. 1988 Model 6
	r_lmn[2][2][0] = 4.83E+07
	r_lmn[3][3][0] = 3.96E+07
	r_lmn[4][4][0] = 3.81E+07
	r_lmn[5][5][0] = 3.79E+07
	# Marley et al. 1988 Model 5
	r_lmn[2][2][1] = 4.77E+07
	r_lmn[3][3][1] = 3.95E+07
	r_lmn[4][4][1] = 3.80E+07
	r_lmn[5][5][1] = 3.79E+07
	# Marley et al. 1988 Model 4
	r_lmn[2][2][2] = 4.68E+07
	r_lmn[3][3][2] = 3.91E+07
	r_lmn[4][4][2] = 3.78E+07
	r_lmn[5][5][2] = 3.78E+07
	return r_lmn, Omega_pat, dim_switch, m_max, ell 

def calculate_frac_of_breakup(planet, G=6.6743E-11):
	R, M, Omega, J2, J4, J6 = import_planet(planet)
	print(Omega / np.sqrt(G * M / (R**3)))

def ballpark_resonance_location(planet, j, G=6.6743E-11):
	R, M, Omega, J2, J4, J6 = import_planet(planet)
	n = j * Omega
	a_res = ((G * M) / (n**2))**(1/3)
	print(a_res)


#calculate_frac_of_breakup('Uranus')
#calculate_frac_of_breakup('Neptune')
#ballpark_resonance_location('Neptune', 2)