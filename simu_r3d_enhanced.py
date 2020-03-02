#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# author: MingChao Ji
# email: mingchao.ji@fysik.su.se
# date created: 2020-01-15 13:26:45
# last modified: 2020-03-02 14:20:58
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from desiree import *
from matplotlib import cm
from mpl_toolkits import mplot3d

import random
import numpy as np
import scipy.stats as st
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
mpl.rc('pdf', fonttype=42)
mpl.rc('axes.spines', left=True, top=True, right=True, bottom=True)
mpl.rcParams.update({'font.family': 'serif', 'font.serif': 'Times', 'mathtext.fontset': 'cm', 'font.size': 12})

# ker selection
idk = 0

# check time limit with tac convert [1] as well as neutrals' position if on imd;
# 3d plot; if save on disk
exp_mode, fig3d, save2fig = True, False, True

# calc ker_speed shared by anion and cation
anion_speed_ker = np.sqrt(2 * ker[idk] * const.eV_J / (anion_mass * (anion_mass / cation_mass + 1)))
cation_speed_ker = anion_mass * anion_speed_ker / cation_mass

anion_ke_ker = func.speed2ke(anion_speed_ker, anion_mass)
cation_ke_ker = func.speed2ke(cation_speed_ker, cation_mass)

num_ions = 20000

ion_pos_dt_cent = []
anion_speed_x, anion_speed_y, anion_speed_z = [], [], []
cation_speed_x, cation_speed_y, cation_speed_z = [], [], []

for n in range(num_ions):
    pos_dt_cent = random.uniform(dts_to_imd[dt_first], dts_to_imd[dt_last + 1])

    z_frac = random.gauss(0, 1)
    x_frac = random.gauss(0, 1)
    y_frac = random.gauss(0, 1)

    sum_frac = np.sqrt(z_frac * z_frac + x_frac * x_frac + y_frac * y_frac)

    anion_spd_z = anion_speed_ker * z_frac / sum_frac
    anion_spd_x = anion_speed_ker * x_frac / sum_frac
    anion_spd_y = anion_speed_ker * y_frac / sum_frac

    # print(anion_spd_z, anion_spd_x, anion_spd_y)

    cation_spd_x = -anion_spd_x * anion_mass / cation_mass
    cation_spd_y = -anion_spd_y * anion_mass / cation_mass
    cation_spd_z = -anion_spd_z * anion_mass / cation_mass

    ion_pos_dt_cent.append(pos_dt_cent)

    anion_speed_x.append(anion_spd_x)
    anion_speed_y.append(anion_spd_y)
    anion_speed_z.append(anion_spd_z)

    cation_speed_x.append(cation_spd_x)
    cation_speed_y.append(cation_spd_y)
    cation_speed_z.append(cation_spd_z)

# for MN occurs in dt_center
# cal ion speed before MN
anion_ke_dt_cent = anion_ke - dt_cent_volt * anion_charge
anion_speed_dt_cent = func.ke2speed(anion_ke_dt_cent, anion_mass)

cation_ke_dt_cent = cation_ke - dt_cent_volt * cation_charge
cation_speed_dt_cent = func.ke2speed(cation_ke_dt_cent, cation_mass)

# cal ion speed after MN
anion_speed_dt_cent_fin = anion_speed_dt_cent + anion_speed_z
cation_speed_dt_cent_fin = cation_speed_dt_cent + cation_speed_z

# the slower speed will come to detector later
late_ion_speed = [min(anion_speed_dt_cent_fin[i], cation_speed_dt_cent_fin[i]) for i in range(num_ions)]

# cal time separation, position separation
anion_tof_dt_cent = ion_pos_dt_cent / anion_speed_dt_cent_fin
cation_tof_dt_cent = ion_pos_dt_cent / cation_speed_dt_cent_fin
delta_t_dt_cent = abs(anion_tof_dt_cent - cation_tof_dt_cent)

# unit: m
anion_x = anion_tof_dt_cent * anion_speed_x
anion_y = anion_tof_dt_cent * anion_speed_y

cation_x = cation_tof_dt_cent * cation_speed_x
cation_y = cation_tof_dt_cent * cation_speed_y

# simu_data = np.array([delta_t_dt_cent, anion_x, anion_y, cation_x, cation_y]).T

r_dt_cent = np.sqrt((cation_x - anion_x) ** 2 + (cation_y - anion_y) ** 2)

r3d = np.sqrt(r_dt_cent ** 2 + (late_ion_speed * delta_t_dt_cent) ** 2)

if exp_mode is True:
    imd_radius = img_size * 0.5 * 1.0e-3
    exp_diss = []
    for i in range(num_ions):
        # check camac dt limit and spot position
        anion_r0 = np.sqrt(anion_x[i] * anion_x[i] + anion_y[i] * anion_y[i])
        cation_r0 = np.sqrt(cation_x[i] * cation_x[i] + cation_y[i] * cation_y[i])
        if np.any([anion_r0 >= imd_radius, cation_r0 >= imd_radius, delta_t_dt_cent[i] >= tac_convert[1] * 1.0e-9]):
            exp_diss.append(i)

    # qualified data
    delta_t_dt_cent = [delta_t_dt_cent[j] for j in range(num_ions) if j not in exp_diss]
    anion_x = [anion_x[j] for j in range(num_ions) if j not in exp_diss]
    anion_y = [anion_y[j] for j in range(num_ions) if j not in exp_diss]
    cation_x = [cation_x[j] for j in range(num_ions) if j not in exp_diss]
    cation_y = [cation_y[j] for j in range(num_ions) if j not in exp_diss]

    r_dt_cent = [r_dt_cent[k] for k in range(num_ions) if k not in exp_diss]
    r3d = [r3d[l] for l in range(num_ions) if l not in exp_diss]

# convert anion and cation position to pixel
anion_x_pix = np.array(anion_x) * 1.0e3 / img_resolution
anion_y_pix = np.array(anion_y) * 1.0e3 / img_resolution

cation_x_pix = np.array(cation_x) * 1.0e3 / img_resolution
cation_y_pix = np.array(cation_y) * 1.0e3 / img_resolution

center_x_pix = 0.5 * (anion_x_pix + cation_x_pix)
center_y_pix = 0.5 * (anion_y_pix + cation_y_pix)

both_x_pix = np.concatenate((anion_x_pix, cation_x_pix), axis=0) + 0.5 * img_pixel
both_y_pix = np.concatenate((anion_y_pix, cation_y_pix), axis=0) + 0.5 * img_pixel

fig = plt.figure(figsize=(8, 12))
fig.suptitle(''.join(['HE: ', str(he_extra_volt/1.0e3), ' kV, LE: ', str(le_extra_volt / 1.0e3),
                     ' kV, DT ', str(dt_first), '-', str(dt_last), ': ', str(int(dt_cent_volt)), ' V']), fontsize=14)
fig.subplots_adjust(left=0.1, bottom=0.08, right=0.95, top=0.95, wspace=0.25, hspace=0.25)

font_big = {'family': 'serif', 'color': 'k', 'size': 14}

ax1 = fig.add_subplot(3, 2, 1)
ax2 = fig.add_subplot(3, 2, 2)
ax3 = fig.add_subplot(3, 2, 3)
ax4 = fig.add_subplot(3, 2, 4)
ax5 = fig.add_subplot(3, 2, 5)

ax1.hist(np.array(delta_t_dt_cent) * 1.0e9, bins=num_bins, histtype='step',
         label=''.join(['dt: ', str(dt_first), '-', str(dt_last)]), alpha=0.75)
if exp_mode is True:
    ax1.set_xlim(0, 250)
    ax1.xaxis.set_major_locator(ticker.MultipleLocator(25))
    ax1.xaxis.set_minor_locator(ticker.MultipleLocator(5))
else:
    ax1.set_xlim(0, )
ax1.set_xlabel('t (ns)')
ax1.set_ylim(0, )
ax1.set_ylabel('Yield (arb. unit)')

ax1.legend(loc='best', frameon=False)

ax3.hist(np.array(r_dt_cent) * 1.0e3, bins=num_bins, histtype='step',
         label=''.join(['dt: ', str(dt_first), '-', str(dt_last)]), alpha=0.75)
if exp_mode is True:
    ax3.set_xlim(0, 70)
    ax3.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax3.xaxis.set_minor_locator(ticker.MultipleLocator(2))
else:
    ax3.set_xlim(0, )
ax3.set_xlabel('r (mm)')
ax3.set_ylim(0, )
ax3.set_ylabel('Yield (arb. unit)')

ax3.legend(loc='best', frameon=False)

ax5.hist(np.array(r3d) * 1.0e3, bins=num_bins, histtype='step', label='r3d', alpha=0.75)
if exp_mode is True:
    ax5.set_xlim(0, 70)
    ax5.xaxis.set_major_locator(ticker.MultipleLocator(10))
    ax5.xaxis.set_minor_locator(ticker.MultipleLocator(2))
else:
    ax5.set_xlim(0, )
ax5.set_xlabel('r3d (mm)')
ax5.set_ylim(0, )
ax5.set_ylabel('Yield (arb. unit)')

ax5.legend(loc='best', frameon=False)

ax2.hist2d(anion_x_pix + 0.5 * img_pixel, anion_y_pix + 0.5 * img_pixel, bins=(num_bins, num_bins),
           cmin=img_count_diss)
plt.text(0.95, 0.95, 'from Anion', fontdict=font_big, ha='right', va='top', transform=ax2.transAxes)
circle1 = plt.Circle((0.5 * img_pixel, 0.5 * img_pixel), radius=0.5 * img_pixel, color='red', linewidth=1,
                     fill=False, alpha=0.8)
ax2.add_patch(circle1)
if exp_mode is True:
    ax2.set_xlim(0, img_pixel)
    ax2.set_ylim(0, img_pixel)
    ax2.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax2.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax2.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax2.yaxis.set_minor_locator(ticker.MultipleLocator(10))
# ax2.grid(b=None, which='major', axis='both', color='k', linestyle='-', linewidth=0.5)
# ax2.grid(b=None, which='minor', axis='both', color='k', linestyle='-.', linewidth=0.5)
ax2.set_ylabel('Y (pixel)')
ax2.set_xlabel('X (pixel)')

ax4.hist2d(cation_x_pix + 0.5 * img_pixel, cation_y_pix + 0.5 * img_pixel, bins=(num_bins, num_bins),
           cmin=img_count_diss)
plt.text(0.95, 0.95, 'from Cation', fontdict=font_big, ha='right', va='top', transform=ax4.transAxes)
circle2 = plt.Circle((0.5 * img_pixel, 0.5 * img_pixel), radius=0.5 * img_pixel, color='red', linewidth=1,
                     fill=False, alpha=0.8)
ax4.add_patch(circle2)
if exp_mode is True:
    ax4.set_xlim(0, img_pixel)
    ax4.set_ylim(0, img_pixel)
    ax4.xaxis.set_major_locator(ticker.MultipleLocator(50))
    ax4.xaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax4.yaxis.set_major_locator(ticker.MultipleLocator(50))
    ax4.yaxis.set_minor_locator(ticker.MultipleLocator(10))
ax4.set_xlabel('X (pixel)')
ax4.set_ylabel('Y (pixel)')

if fig3d is False:
    ax6 = fig.add_subplot(3, 2, 6)
    ax6.hist2d(center_x_pix + 0.5 * img_pixel, center_y_pix + 0.5 * img_pixel, bins=(num_bins, num_bins),
               cmin=img_count_diss)
    plt.text(0.95, 0.95, 'Center Position', fontdict=font_big, ha='right', va='top', transform=ax6.transAxes)
    circle3 = plt.Circle((0.5 * img_pixel, 0.5 * img_pixel), radius=0.5 * img_pixel, color='red', linewidth=1,
                         fill=False, alpha=0.8)
    ax6.add_patch(circle3)
    if exp_mode is True:
        ax6.set_xlim(0, img_pixel)
        ax6.set_ylim(0, img_pixel)
        ax6.xaxis.set_major_locator(ticker.MultipleLocator(50))
        ax6.xaxis.set_minor_locator(ticker.MultipleLocator(10))
        ax6.yaxis.set_major_locator(ticker.MultipleLocator(50))
        ax6.yaxis.set_minor_locator(ticker.MultipleLocator(10))
    ax6.set_xlabel('X (pixel)')
    ax6.set_ylabel('Y (pixel)')
else:
    xx, yy = np.mgrid[0:img_pixel:100j, 0:img_pixel:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([both_x_pix, both_y_pix])
    kernel = st.gaussian_kde(values)

    f = np.reshape(kernel(positions).T, xx.shape)

    ax6 = fig.add_subplot(3, 2, 6, projection='3d')
    surf = ax6.plot_surface(xx, yy, f * 1.0e4, cmap=cm.coolwarm, linewidth=0, antialiased=False)
    ax6.set_xlabel('X')
    ax6.set_ylabel('Y')
    # ax6.set_zlabel('Z')
    ax6.set_title('Both')
    # fig.colorbar(surf, shrink=0.5, aspect=5) # add color bar indicating the Z
    ax6.view_init(60, 45)

# save file or, not
if save2fig is True:
    fig.savefig(''.join(['Simu_MN_DT', str(dt_first), '-', str(dt_last), '_Enh_n=', str(idk + 1), '.png']), dpi=150)
else:
    plt.show()
