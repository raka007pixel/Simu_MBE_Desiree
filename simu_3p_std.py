#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# |m|i|n|g|c|h|a|o|j|i| @ 2019-03-10 13:43|
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from desiree import *
from matplotlib import cm
from mpl_toolkits import mplot3d
from scipy.special import comb, perm
from itertools import combinations, permutations

import random
import scipy.stats as st
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

mpl.rc('pdf', fonttype=42)
mpl.rc('axes.spines', left=True, top=True, right=True, bottom=True)
mpl.rcParams.update({'font.family': 'serif', 'font.serif': 'Times', 'mathtext.fontset': 'cm', 'font.size': 12})

# ker selection
idk = 0

# check time limit with tac convert [1] as well as neutrals' position if on imd; if save on disk
exp_mode, save2fig = True, True

# gravity center selection, values given by xc, yc and r in order. unit: mm
img_slt_circ = [0, 0, 5.5]

num_ions = 50000

# calc av ion speed after MN depending on drift tube voltage
anion_ke_dt_cent = anion_ke - dt_cent_volt * anion_charge
cation_ke_dt_cent = cation_ke - dt_cent_volt * cation_charge

anion_speed_dt_cent = func.ke2speed(anion_ke_dt_cent, anion_mass)
cation_speed_dt_cent = func.ke2speed(cation_ke_dt_cent, cation_mass)

av_ion_speed = np.mean([anion_speed_dt_cent, cation_speed_dt_cent])

av_L = av_ion_speed / ms_ctr_to_imd
av2_L = av_ion_speed * av_ion_speed / ms_ctr_to_imd

# max energy shared by He
e1_max = sum([neut_amu[1], neut_amu[2]]) / sum([neut_amu[0], neut_amu[1], neut_amu[2]])

# array to save img info, tof info
img_array = [[[] for n in range(2)] for l in range(n_neut)]
dist_array = [[] for l in range(n_neut)]
tof_array = [[] for l in range(n_neut)]
dt_array = [[] for l in range(n_neut)]

p1_ker, p2_ker, p3_ker = [], [], []
area_imd = []

for n in range(num_ions):
    # event position
    pos_dt_bias = random.uniform(dts_to_imd[dt_first], dts_to_imd[dt_last + 1])

    # ker taken by He:
    e1 = random.uniform(0, e1_max)
    p1_ker.append(e1)

 # doesn't work


# plot
fig, ax = plt.subplots(2, 4, figsize=(16, 8))
fig.suptitle(''.join(['HE: ', str(he_extra_volt / 1.0e3), ' kV, LE: ', str(le_extra_volt / 1.0e3),
                      ' kV, DT ', str(dt_first), '-', str(dt_last), ': ', str(int(dt_cent_volt)), ' V']), fontsize=14)
fig.subplots_adjust(left=0.05, bottom=0.08, right=0.95, top=0.92, wspace=0.25, hspace=0.25)

font_big = {'family': 'serif', 'color': 'k', 'size': 14}


ax[0, 0].set_title('ker sharing ratio')
ax[0, 0].hist2d(p1_ker, p2_ker, bins=(num_bins, num_bins), cmin=img_count_diss, cmap=cm.autumn)
ax[0, 0].hist2d(p3_ker, p2_ker, bins=(num_bins, num_bins), cmin=img_count_diss, cmap=cm.winter)
ax[0, 0].set_xlim(0, 1)
ax[0, 0].set_ylim(0, 1)
ax[0, 0].set_xlabel('a (He), c (H)')
ax[0, 0].set_ylabel('b (D)')


ax[0, 1].set_title('img')
# ax[0, 1].hist2d(img_x_pix, img_y_pix, bins=(num_bins, num_bins), cmin=img_count_diss)
ax[0, 1].hist2d(p1_x_pix, p1_y_pix, bins=(num_bins, num_bins), cmin=img_count_diss)
ax[0, 1].hist2d(p2_x_pix, p2_y_pix, bins=(num_bins, num_bins), cmin=img_count_diss)
ax[0, 1].hist2d(p3_x_pix, p3_y_pix, bins=(num_bins, num_bins), cmin=img_count_diss)
circle1 = plt.Circle((0.5 * img_pixel, 0.5 * img_pixel), radius=0.5 * img_pixel, color='red', linewidth=1,
                     fill=False, alpha=0.8)
ax[0, 1].add_patch(circle1)
if exp_mode is True:
    ax[0, 1].set_xlim(0, img_pixel)
    ax[0, 1].set_ylim(0, img_pixel)
    ax[0, 1].xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[0, 1].xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax[0, 1].yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax[0, 1].yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax[0, 1].set_ylabel('Y (pixel)')
ax[0, 1].set_xlabel('X (pixel)')


ax[0, 2].set_title('arrival time')
ax[0, 2].hist(np.array(tof_array[0]) * 1.0e6, bins=num_bins, histtype='step', label='He', alpha=0.75)
ax[0, 2].hist(np.array(tof_array[1]) * 1.0e6, bins=num_bins, histtype='step', label='D', alpha=0.75)
ax[0, 2].hist(np.array(tof_array[2]) * 1.0e6, bins=num_bins, histtype='step', label='H', alpha=0.75)
ax[0, 2].set_xlabel(''.join(['ar. time (', str(chr(956)), 's)']))
ax[0, 2].legend(loc='best', frameon=False)


ax[0, 3].set_title('area')
# ax[0, 3].hist(area_pix, bins=num_bins, histtype='step', alpha=0.75)
ax[0, 3].hist(area_imd, bins=num_bins, histtype='step', alpha=0.75)
ax[0, 3].set_xlim(0, )
ax[0, 3].set_xlabel(r'area ($\rm{mm^{2}}$)')

ax[1, 0].set_title('Dalitz coords')
plt.text(0.5, 0.5, 'too fancy to plot', fontdict=font_big, ha='center', va='center', transform=ax[1, 0].transAxes)


ax[1, 1].set_title('r')
ax[1, 1].hist(np.array(dist_array[0]) * 1.0e3, bins=num_bins, histtype='step', label='D-He', alpha=0.75)
ax[1, 1].hist(np.array(dist_array[1]) * 1.0e3, bins=num_bins, histtype='step', label='H-He', alpha=0.75)
ax[1, 1].hist(np.array(dist_array[2]) * 1.0e3, bins=num_bins, histtype='step', label='H-D', alpha=0.75)
ax[1, 1].set_xlabel('distance (mm)')
ax[1, 1].legend(loc='best', frameon=False)


ax[1, 2].set_title('delta_t')
ax[1, 2].hist(np.array(dt_array[0]) * 1.0e9, bins=num_bins, histtype='step', label='D-He', alpha=0.75)
ax[1, 2].hist(np.array(dt_array[1]) * 1.0e9, bins=num_bins, histtype='step', label='H-He', alpha=0.75)
ax[1, 2].hist(np.array(dt_array[2]) * 1.0e9, bins=num_bins, histtype='step', label='H-D', alpha=0.75)
ax[1, 2].set_xlabel('ar. time diff (ns)')
ax[1, 2].legend(loc='best', frameon=False)


ax[1, 3].set_title('sum ker')
ax[1, 3].hist(calc_ker_perm, bins=num_bins, histtype='step', color='dodgerblue')
ax[1, 3].set_xlim(0, 50)
plt.text(0.95, 0.95, str(len_calc_ker), fontdict=font_big, ha='right', va='top', transform=ax[1, 3].transAxes)
# plt.text(0.25, 0.25, str(len_ker_slt), fontdict=font_big, ha='right', va='top', transform=ax[1, 3].transAxes)
ax[1, 3].set_xlabel('energy (eV)')


# save file or, not
if save2fig is True:
    fig.savefig(''.join(['Simu_MN3p_DT', str(dt_first), '-', str(dt_last), '_n=', str(idk), '.png']), dpi=150)
else:
    plt.show()
