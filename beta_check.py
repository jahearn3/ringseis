# beta check
#
# Author: Joseph A'Hearn
# Created 02/18/2021
#
# This program evaluates the integral of the displacement eigenfunction 
#   to check the consistency of the beta value provided by GYRE
#

import numpy as np 
import get_gyre_output as gg 
import scipy.integrate as integrate 
import matplotlib.pyplot as plt 
import plot_assistant as pa 
import matplotlib.ticker as ticker

def plot_xi(x, xi_r, xi_h, rho, j, l, n, eul_rho, black=False):
	# Plot displacement eigenfunctions 
	fig, ax = pa.initialize_plot(black)
	fig, ax = pa.title_and_axes(fig, ax, '', 'r/R', '')
	ax.plot(x, xi_r, label=r'$\xi_r$')
	ax.plot(x, xi_h, label=r'$\xi_h$')
	#if(l > 5):
	#	plt.yscale('log')
	#print('max x: ' + str(np.amax(x)))
	#ax.plot(x, ((2 * xi_r * xi_h) + (xi_h**2)) * rho * (x**2), label=r'$(2 \xi_r \xi_h + \xi_h^2) \rho (\frac{r}{R})^2$')
	#ax.plot(x, ((xi_r**2) + (l * (l + 1) * (xi_h**2))) * rho * (x**2), label=r'$(\xi_r^2 + \ell (\ell + 1) \xi_h^2) \rho (\frac{r}{R})^2$')
	ax.legend(loc=0)
	if(l == 2):
		ax.figure.savefig('xi_plot_' + str(j) + '_' + str(l) +'_' + str(n) + '.pdf')
	else:
		pa.save_and_clear_plot(fig, ax, filename='xi_plot_' + str(j) + '_' + str(l) +'_' + str(n), black=black)
	write_to_file = True
	if(write_to_file == True):
		with open('xi_data_' + str(l) + '.txt', 'a') as myfile:
			for i in range(len(x)):
				myfile.write('%.8f' % (x[i])) 
				myfile.write('\t\t%.8f' % (xi_r[i]))
				myfile.write('\t\t%.8f' % (xi_h[i]))
				myfile.write('\t\t%.8f' % (eul_rho[i]))
				myfile.write('\n')

def plot_xi_for_paper(l, planet='Uranus', rowspan=12, colspan=12):
	files = []
	for i in range(len(l)):
		files.append('xi_data_' + str(l[i]) + '.txt') 
	c = ['b', 'orange', 'g', 'r']
	fig = plt.figure(figsize=(8,12))
	ax1 = plt.subplot2grid((12,12),(0, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	ax2 = plt.subplot2grid((12,12),(4, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	ax3 = plt.subplot2grid((12,12),(8, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	fig, ax1 = pa.title_and_axes(fig, ax1, '', '', r'$\Delta r / R$', 0, 1)
	fig, ax2 = pa.title_and_axes(fig, ax2, '', '', r'$\rho\'(r)$', 0, 1)
	fig, ax3 = pa.title_and_axes(fig, ax3, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho\'(r) r^{\ell + 2}$', 0, 1)
	ax1.tick_params(axis='x', labelsize=18)
	ax1.tick_params(axis='y', labelsize=18)
	ax2.tick_params(axis='x', labelsize=18)
	ax2.tick_params(axis='y', labelsize=18)
	ax3.tick_params(axis='x', labelsize=18)
	ax3.tick_params(axis='y', labelsize=18)
	ax1.xaxis.set_major_locator(ticker.NullLocator())
	ax2.xaxis.set_major_locator(ticker.NullLocator())
	for i in range(len(files)):
		data = np.loadtxt(files[i])
		r    = data[:,0]
		xi_r = data[:,1] 						 
		xi_h = data[:,2] 
		rho_prime = data[:,3] 

		xi_max = np.amax(xi_r)

		ax1.plot(r, xi_r / xi_max, c=c[i])
		#ax1.plot(r, xi_h / xi_max, linestyle='dashed')
		ax2.plot(r, rho_prime, c=c[i])
		ax3.plot(r, rho_prime * (r ** (l[i]+2)), c=c[i])

	x_leg = 0.02
	y_leg = 0.9
	spacing = 0.1
	fs = 22
	ax1.text(x_leg, y_leg,                 r'$\ell=2$',  fontsize=fs, c='b')
	ax1.text(x_leg, y_leg - (1 * spacing), r'$\ell=6$',  fontsize=fs, c='orange')
	ax1.text(x_leg, y_leg - (2 * spacing), r'$\ell=10$', fontsize=fs, c='g')
	ax1.text(x_leg, y_leg - (3 * spacing), r'$\ell=14$', fontsize=fs, c='r')
	ax1.figure.savefig('xi_plot_for_paper.pdf')

def main(j=2, l=2, n=0, pathname='/Users/josephahearn/Downloads/work/eigs/j', rule='trapezoidal', plots=True, g_modes=False):
	x, xi_r, xi_h, rho, beta_gyre, R, eul_rho = gg.xi_data(j, l, n, pathname, g_modes)
	#print(beta_gyre)
	if(plots == True):
		plot_xi(x, xi_r, xi_h, rho, j, l, n, eul_rho)
	r = x * R
	C_numerator_integrand = ((2 * xi_r * xi_h) + (xi_h**2)) * rho * (x**2)
	C_denominator_integrand = ((xi_r**2) + (l * (l + 1) * (xi_h**2))) * rho * (x**2)
	if(rule == 'trapezoidal'):
		# Compute the area using the composite trapezoidal rule.
		area1 = np.trapz(C_numerator_integrand, dx=0.01)
		#print("area = " + str(area1))
		area3 = np.trapz(C_denominator_integrand, dx=0.01)
		#print("area = " + str(area3))
		beta1 = 1 - (area1 / area3)
		#print(beta1)
	elif(rule == 'simpson'):
		# Compute the area using the composite Simpson's rule.
		area2 = integrate.simps(C_numerator_integrand, dx=0.01)
		#print("area = " + str(area2))
		area4 = integrate.simps(C_denominator_integrand, dx=0.01)
		#print("area = " + str(area4))
		#beta2 = 1 - (area2 / area4)
		#print(beta2)
	elif(rule == 'quad'): # I don't think this worked
		x2 = lambda r: rho * r**2
		i = integrate.quad(x2, 0, 4)
		#print(i)
		C_numerator =   integrate.quad(lambda r: C_numerator_integrand,   0, R)
		C_denominator = integrate.quad(lambda r: C_denominator_integrand, 0, R)
		C = C_numerator / C_denominator
		beta = 1 - C
		#print(beta)

	return beta1

#main(7, 2, 0)
#plot_xi_for_paper([2, 6, 10, 14])
