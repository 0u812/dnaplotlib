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
lac_color = blue
tet_color = bright_orange
ci_color = red
def make_operon(origin, reverse, name, p_color, cds_color):
	p = {'name':'P_'+name, 'start':origin+(75 if reverse else 65), 'end':origin+(65 if reverse else 75), 'fwd': not reverse, 'type':'Promoter', 'opts': {'color':p_color}}
	rbs = {'name':'RBS', 'start': origin+(10 if reverse else 5), 'end': origin+(5 if reverse else 10), 'fwd': not reverse, 'type':'RBS', 'opts':{'linewidth': 0, 'color':[0.0, 0.0, 0.0]}}
	cds = {
		'name': name+'R',
		'start': origin+(50 if reverse else 20),
		'end':   origin+(20 if reverse else 50),
		'fwd':   not reverse,
		'type': 'CDS',
		'opts': {
			'label': name+'R',
			'fontsize': 8,
			'label_y_offset': 0,
			'label_x_offset': -2,
			'label_style':'italic',
			'color':cds_color}
		}
	term = {'name':'Term', 'start':origin+(15 if reverse else 5), 'end': origin+(5 if reverse else 15), 'fwd': not reverse, 'type':'Terminator'}
	return p,rbs,cds,term

p1, rbs1, tetr, term1 = make_operon(origin=0, reverse=True,  name='lac', p_color=tet_color, cds_color=lac_color)
p2, rbs2, laci, term2 = make_operon(origin=0, reverse=False, name='tet', p_color=lac_color, cds_color=tet_color)

# lac_repress = {'from_part':laci, 'to_part':plac, 'type':'Repression', 'opts':{'linewidth':1, 'color':plac_color}}
# gamma_repress = {'from_part':gamma, 'to_part':pgamma, 'type':'Repression', 'opts':{'linewidth':1, 'color':pci_color}}
# tet_repress = {'from_part':tetr, 'to_part':ptet, 'type':'Repression', 'opts':{'linewidth':1, 'color':ptet_color}}

if __name__ == '__main__':
	plt.close()
	plt.figure(figsize=(3.5, 1.5))
	gs = gridspec.GridSpec(1, 1, height_ratios=[1])

	# Plot of repressilator circuit
	ax = plt.subplot(gs[0])
	# dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1, pgamma, rbs2, laci, term2, ptet, rbs3, gamma, term3]],
	# 			[[lac_repress, gamma_repress, tet_repress]])
	# dnaplotlib.plot_sbol_designs([ax], [[plac, rbs1, tetr, term1]])
	dnaplotlib.plot_sbol_designs([ax], [[term1, tetr, rbs1, p1, p2, rbs2, laci, term2]])
	ax.set_ylim([-15, 31])

	# Update subplot spacing
	plt.subplots_adjust(hspace=0.4, left=0.12, right=0.95, top=0.99, bottom=0.01)

	# Save the figure
	plt.savefig('bistable.pdf', transparent=True)
	# plt.savefig('repressilator_animate.png', dpi=300)
