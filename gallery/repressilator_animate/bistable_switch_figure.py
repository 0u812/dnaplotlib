#!/usr/bin/env python3

import numpy as np
from scipy.integrate import odeint
import dnaplotlib
import matplotlib.pyplot as plt
from matplotlib import gridspec

__author__  = 'Emerson Glassey <eglassey@mit.edu>, Voigt Lab, MIT'
__license__ = 'MIT'
__version__ = '1.0'

# Initialize Parts
blue = (0.38, 0.65, 0.87)
bright_orange = (1.00, 0.75, 0.17)
red = tuple(x/256.+0.1 for x in (214, 39, 40))
plac_color = blue
ptet_color = bright_orange
pci_color = red
plac = {'name':'P_lac', 'start':40, 'end':35, 'type':'Promoter', 'opts': {'color':plac_color}}
rbs1 = {'name':'RBS', 'start':30, 'end':25, 'type':'RBS', 'opts':{'linewidth': 0, 'color':[0.0, 0.0, 0.0]}}
tetr = {
	'name': 'tetR',
	'start': 50,
	'end':   20,
	'fwd':   False,
	'type': 'CDS',
	'opts': {
		'label': 'tetR',
		'fontsize': 8,
		'label_y_offset': 0,
		'label_x_offset': -2,
		'label_style':'italic',
		'color':ptet_color}
	}
term1 = {'name':'Term', 'start':10, 'end':5, 'type':'Terminator'}

pgamma = {'name':'P_gamma', 'start':56, 'end':65, 'type':'Promoter', 'opts': {'color':pci_color}}
rbs2 = {'name':'RBS', 'start':66, 'end':75, 'type':'RBS', 'opts':{'linewidth': 0, 'color':[0.0, 0.0, 0.0]}}
laci = {'name':'lacI', 'start':76, 'end':95, 'type':'CDS', 'opts':{'label': 'lacI', 'fontsize': 8,  'label_y_offset': 0, 'label_x_offset': -2, 'label_style':'italic', 'color':plac_color}}
term2 = {'name':'Term', 'start':96, 'end':110, 'type':'Terminator'}
ptet = {'name':'P_tet', 'start':111, 'end':120, 'type':'Promoter', 'opts': {'color':ptet_color}}
rbs3 = {'name':'RBS', 'start':121, 'end':130, 'type':'RBS', 'opts':{'linewidth': 0, 'color':[0.0, 0.0, 0.0]}}
gamma = {'name':'gamma', 'start':131, 'end':150, 'type':'CDS', 'opts':{'label': 'cI', 'fontsize': 8, 'label_y_offset': 0, 'label_x_offset': -1, 'label_style':'italic', 'color':pci_color}}
term3 = {'name':'Term', 'start':151, 'end':165, 'type':'Terminator'}

lac_repress = {'from_part':laci, 'to_part':plac, 'type':'Repression', 'opts':{'linewidth':1, 'color':plac_color}}
gamma_repress = {'from_part':gamma, 'to_part':pgamma, 'type':'Repression', 'opts':{'linewidth':1, 'color':pci_color}}
tet_repress = {'from_part':tetr, 'to_part':ptet, 'type':'Repression', 'opts':{'linewidth':1, 'color':ptet_color}}

def plot_construct(ax, t, ymtet, ymlac, ymgamma, ytet, ylac, ygamma):
	tind = int(t*10)
	exp_lims = (1.0, 4.0)
	ax.set_title('t = {}'.format(t), fontsize=8)
	# Set color for each of the CDSs
	tetr['opts']['color'] = [rescale(1 - expression(ymlac[tind], exp_lims), (x, 1.0)) for x in ptet_color]
	laci['opts']['color'] = [rescale(1 - expression(ymlac[tind], exp_lims), (x, 1.0)) for x in plac_color]
	gamma['opts']['color'] = [rescale(1 - expression(ymlac[tind], exp_lims), (x, 1.0)) for x in pci_color]
	# Set transparency for each of the regulatory lines
	lac_repress['opts']['color'] = [*plac_color,
								rescale(repression(ylac[tind], 2.0, 8), (0.2, 1.0))]
	gamma_repress['opts']['color'] = [*pci_color,
								rescale(repression(ygamma[tind], 2.0, 8), (0.2, 1.0))]
	tet_repress['opts']['color'] = [*ptet_color,
								rescale(repression(ytet[tind], 2.0, 8), (0.2, 1.0))]
	# Set width for each of the regulatory lines
	lac_repress['opts']['linewidth'] = rescale(repression(ylac[tind], 2.0, 8), (0.5, 2.0))
	gamma_repress['opts']['linewidth'] = rescale(repression(ygamma[tind], 2.0, 8), (0.5, 2.0))
	tet_repress['opts']['linewidth'] = rescale(repression(ytet[tind], 2.0, 8), (0.5, 2.0))
	dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1, pgamma, rbs2, laci, term2, ptet, rbs3, gamma, term3]],
				[[lac_repress, gamma_repress, tet_repress]])
	ax.set_ylim([-10, 31])

if __name__ == '__main__':
	plt.close()
	plt.figure(figsize=(3.5, 1.5))
	gs = gridspec.GridSpec(1, 1, height_ratios=[1])

	# Plot of repressilator circuit
	ax = plt.subplot(gs[0])
	dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1, pgamma, rbs2, laci, term2, ptet, rbs3, gamma, term3]],
				[[lac_repress, gamma_repress, tet_repress]])
	ax.set_ylim([-10, 31])

	# Update subplot spacing
	plt.subplots_adjust(hspace=0.4, left=0.12, right=0.95, top=0.99, bottom=0.01)

	# Save the figure
	plt.savefig('bistable.pdf', transparent=True)
	# plt.savefig('repressilator_animate.png', dpi=300)
