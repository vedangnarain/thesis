#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Fri Aug 23 15:47:27 2024

@author: Vedang Narain (vedang.narain@msdtc.ox.ac.uk)

Plot the oxygen conversion.

Parameters and methods sourced from: 
    
[1] D. R. Grimes et al., ‘Estimating oxygen distribution from vasculature in 
three-dimensional tumour tissue’, J. R. Soc. Interface., vol. 13, no. 116, p. 
20160070, Mar. 2014, doi: 10.1098/rsif.2016.0070.

[2] T. D. Lewin, P. K. Maini, E. G. Moros, H. Enderling, and H. M. Byrne, ‘The 
Evolution of Tumour Composition During Fractionated Radiotherapy: Implications 
for Outcome’, Bull Math Biol, vol. 80, no. 5, pp. 1207–1235, May 2018, doi: 
10.1007/s11538-018-0391-9.
"""

# Initialise libraries
# import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np
from pathlib import Path

# Set LaTex-style font
plt.rcParams['mathtext.fontset'] = 'stix'
plt.rcParams['font.family'] = 'STIXGeneral'
plt.rcParams.update({'font.size': 22})

# Import conversion functions
from convert_oxygen import *

# Chaste's conversion is the same
# 3.1e-5*(1.0)/(22.4e-3*((1.0)**3))*3mmHg/1000*1e9
'''
# Enter parameters from paper [1]
K = 22779
rb_pp = 3
hyp_pp = 15
phys_pp = 38
# tumour_boundary_pp = 100

# Plot concentration (nM) vs partial pressure (mmHg)
pp_range = np.linspace(0, phys_pp, 1000)
nM_range = (pp_range/K)/(0.000000001*32)
# nM_range = np.linspace(0, 50000, 100)
# pp_range = nM_range*0.000000001*32*K
plt.figure(figsize=(15, 7.5))
#plt.title('Relationship between concentration and partial pressure (Grimes et al., 2014)')
plt.plot(pp_range, nM_range, c='black')
#plt.axhline(100, alpha=0.5, c='g', ls='--', label='tumour boundary')
#plt.axhline(hypoxic_pp, alpha=0.5, c='orange', ls='--', label='hypoxia')
#plt.axhline(anoxic_pp, alpha=0.5, c='red', ls='--', label='anoxia')
plt.fill_between([0, hyp_pp], 0, convert_mmHg_to_nM(hyp_pp), facecolor='blue', alpha=0.75)
plt.fill_between([0, rb_pp], 0, convert_mmHg_to_nM(rb_pp), facecolor='deepskyblue', alpha=1)
plt.fill_between([0, hyp_pp], convert_mmHg_to_nM(hyp_pp), convert_mmHg_to_nM(phys_pp), facecolor='red', alpha=0.75)
plt.fill_between([hyp_pp, phys_pp], 0, convert_mmHg_to_nM(phys_pp), facecolor='red', alpha=0.75)
# plt.fill_between([hyp_pp, phys_pp], convert_mmHg_to_nM(hyp_pp), convert_mmHg_to_nM(phys_pp ), facecolor='red', alpha=0.75)
# plt.fill_between([0, phys_pp], 0, convert_mmHg_to_nM(phys_pp ), facecolor='red', alpha=0.75)
# plt.fill_between([convert_mmHg_to_nM(hypoxic_pp), convert_mmHg_to_nM(30)], hypoxic_pp, 30, facecolor='red', alpha=0.75)
# plt.fill_between(pp_range, 0, convert_mmHg_to_nM(anoxic_pp), facecolor='blue', alpha=0.75)
# plt.fill_between(pp_range, convert_mmHg_to_nM(anoxic_pp), convert_mmHg_to_nM(hypoxic_pp), facecolor='deepskyblue', alpha=0.75)
# plt.fill_between(pp_range, convert_mmHg_to_nM(hypoxic_pp), convert_mmHg_to_nM(30), facecolor='red', alpha=0.75)
plt.ylim(0, convert_mmHg_to_nM(phys_pp))
plt.xlim(0,phys_pp)
plt.ylabel('concentration of oxygen (nM)')
plt.xlabel('partial pressure of oxygen (mmHg)')
#plt.legend()

anoxia_patch = mpatches.Patch(color='deepskyblue', label='radiobiological hypoxia')
hypoxia_patch = mpatches.Patch(color='blue', label='physiological hypoxia')
normoxia_patch = mpatches.Patch(color='red', label='physoxia')

plt.legend(handles=[normoxia_patch, hypoxia_patch, anoxia_patch, ])
file_path = Path('~/Desktop/Final Figures/relationship between concentration and partial pressure.png').expanduser()
plt.savefig(file_path, bbox_inches='tight', dpi=500)
plt.show()
'''
# Save image
# file_path = Path('~/Desktop/Final Figures/relationship between concentration and partial pressure.svg').expanduser()
# plt.savefig(file_path, bbox_inches='tight', dpi=500)
