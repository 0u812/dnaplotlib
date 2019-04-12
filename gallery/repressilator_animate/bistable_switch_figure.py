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
def make_operon(origin, reverse):
	plac = {'name':'P_lac', 'start':origin+75, 'end':origin+65, 'fwd': not reverse, 'type':'Promoter', 'opts': {'color':plac_color}}
	rbs1 = {'name':'RBS', 'start': origin+10, 'end': origin+5, 'fwd': not reverse, 'type':'RBS', 'opts':{'linewidth': 0, 'color':[0.0, 0.0, 0.0]}}
	tetr = {
		'name': 'tetR',
		'start': origin+(50 if reverse else 20),
		'end':   origin+(20 if reverse else 50),
		'fwd':   not reverse,
		'type': 'CDS',
		'opts': {
			'label': 'tetR',
			'fontsize': 8,
			'label_y_offset': 0,
			'label_x_offset': -2,
			'label_style':'italic',
			'color':ptet_color}
		}
	term1 = {'name':'Term', 'start':origin+5, 'end': origin+15, 'fwd': not reverse, 'type':'Terminator'}

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

if __name__ == '__main__':
	plt.close()
	plt.figure(figsize=(3.5, 1.5))
	gs = gridspec.GridSpec(1, 1, height_ratios=[1])

	# Plot of repressilator circuit
	ax = plt.subplot(gs[0])
	# dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1, pgamma, rbs2, laci, term2, ptet, rbs3, gamma, term3]],
	# 			[[lac_repress, gamma_repress, tet_repress]])
	# dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1]])
	dnaplotlib.plot_sbol_designs([ax], [[term1, tetr, rbs1, plac]])
	ax.set_ylim([-15, 31])

	# Update subplot spacing
	plt.subplots_adjust(hspace=0.4, left=0.12, right=0.95, top=0.99, bottom=0.01)

	# Save the figure
	plt.savefig('bistable.pdf', transparent=True)
	# plt.savefig('repressilator_animate.png', dpi=300)
