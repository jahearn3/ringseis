# xi and density perturbation
#
# Author: Joseph A'Hearn
# Created 12/21/2021
#
# This program plots xi and the density perturbation from the output of a GYRE simulation 
#   

import numpy as np 
import matplotlib.pyplot as plt 
import matplotlib.ticker as ticker
import compare_models as cm 

def summary_data(path='/Users/josephahearn/Downloads/work/summary.h5'):
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

def xi_data(j, l, n=0, path='/Users/josephahearn/Downloads/work/eigs/j'):
	if(n >= 0):
		sign = '+'
	else:
		sign = '-'
	path += str(j) + '_l' + str(l) + '_n' + sign + str(n) + '.h5'
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

def title_and_axes(fig, ax, title, xlabel, ylabel, xmin=None, xmax=None, ymin=None, ymax=None):
	ax.set_title(title, fontsize=22)
	ax.set_xlabel(xlabel, fontsize=18)
	ax.set_ylabel(ylabel, fontsize=20)
	plt.rc('xtick', labelsize=18)    # fontsize of the tick labels
	plt.rc('ytick', labelsize=18)    # fontsize of the tick labels
	if((xmin != None) and (xmax != None)):
		ax.set_xlim([xmin,xmax])
	if((ymin != None) and (ymax != None)):
		ax.set_ylim([ymin,ymax])
	return fig, ax 

def make_plot(j, l, path, planet='Uranus', rowspan=12, colspan=12):
	#planet = 'Saturn'
	c = ['b', 'orange', 'g', 'r']
	xmax = 1.001
	fig = plt.figure(figsize=(8,12))
	ax1 = plt.subplot2grid((12,12),(0, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	ax2 = plt.subplot2grid((12,12),(4, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	ax3 = plt.subplot2grid((12,12),(8, 1), rowspan=int(rowspan / 3), colspan=colspan - 1, fig=fig)
	#fig, ax1 = title_and_axes(fig, ax1, '', '', r'$\Delta r / R$', 0, xmax)
	fig, ax1 = title_and_axes(fig, ax1, '', '', r'$\xi_r$', 0, xmax)
	#fig, ax2 = title_and_axes(fig, ax2, '', '', r'$\rho\'(r)$ [g/cm$^3$]', 0, xmax)
	fig, ax2 = title_and_axes(fig, ax2, '', '', r'$\rho\'(r)$', 0, xmax)
	fig, ax3 = title_and_axes(fig, ax3, '', r'$r$' + ' [' + r'$R_{' + planet + '}$' + ']', r'$\rho\'(r) r^{\ell + 2}$', 0, xmax)
	ax1.tick_params(axis='x', labelsize=18)
	ax1.tick_params(axis='y', labelsize=18)
	ax2.tick_params(axis='x', labelsize=18)
	ax2.tick_params(axis='y', labelsize=18)
	ax3.tick_params(axis='x', labelsize=18)
	ax3.tick_params(axis='y', labelsize=18)
	ax1.xaxis.set_major_locator(ticker.NullLocator())
	ax2.xaxis.set_major_locator(ticker.NullLocator())
	ax1.spines['top'].set_visible(False)
	ax1.spines['right'].set_visible(False)
	ax2.spines['top'].set_visible(False)
	ax2.spines['right'].set_visible(False)
	ax3.spines['top'].set_visible(False)
	ax3.spines['right'].set_visible(False)
	include_N = True
	if(include_N):
		r1, m, P, T, rho1, Delta, N_sq, Gamma_1, Delta_ad, delta = cm.read_model('Umedium4096model.gyre')
		N = np.sqrt(N_sq)
		Nmax = np.amax(N)
		R = np.amax(r1)
		ax1.plot(r1/R, N/Nmax, linestyle='dotted', c='k')
	for i in range(len(j)):
		r, xi_r, xi_h, rho, beta, R, eul_rho = xi_data(j[i], l[i], path=path + 'eigs/j')
		xi_max = np.amax(xi_r)
		rho_prime = eul_rho * rho
		#print(np.amax(rho))
		ax1.plot(r, xi_r / xi_max, c=c[i])
		#ax1.plot(r, xi_h / xi_max, linestyle='dashed')
		ax2.plot(r, rho_prime, c=c[i])
		ax3.plot(r, rho_prime * (r ** (l[i]+2)), c=c[i])
	x_leg = 0.02
	y_leg = 0.9
	spacing = 0.2
	fs = 30
	ax1.text(x_leg, y_leg,                 r'$\ell=2$',  fontsize=fs, c='b')
	ax1.text(x_leg, y_leg - (1 * spacing), r'$\ell=6$',  fontsize=fs, c='orange')
	ax1.text(x_leg, y_leg - (2 * spacing), r'$\ell=10$', fontsize=fs, c='g')
	ax1.text(x_leg, y_leg - (3 * spacing), r'$\ell=14$', fontsize=fs, c='r')
	ax1.figure.savefig('xi_plot_for_paper.pdf')

def main(path='/Users/josephahearn/Downloads/work/77_Umedium4096_all_modes_found/', filename='summary.h5'):
	pathname = path + filename
	j_, ell, omega, beta = summary_data(pathname)
	l = [2, 6, 10, 14]
	j = []
	for k in range(len(l)):
		for i in range(len(j_)):
			if(l[k] == ell[i]):
				j.append(j_[i])
	make_plot(j, l, path)

main()
#main('/Users/josephahearn/Downloads/work/80_U3_d_all_modes_found/')
#main('/Users/josephahearn/Downloads/work/78_Uthick4096_all_modes_found/')
#main('/Users/josephahearn/Downloads/work/14_Saturn/')

