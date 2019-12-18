#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# |m|i|n|g|c|h|a|o|j|i| @ 2019-03-10 13:43|
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from desiree import *

import random
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

mpl.rc('pdf', fonttype=42)
mpl.rc('axes.spines', left=True, top=True, right=True, bottom=True)
mpl.rcParams.update({'font.family': 'serif', 'font.serif': 'Times', 'mathtext.fontset': 'cm', 'font.size': 14})

# save on disk or not
save2fig = True

ion_pos_dt_cent, ion_pos_pu2, ion_pos_pu1 = [], [], []

anion_speed_axial, anion_speed_radial = [], []
cation_speed_axial, cation_speed_radial = [], []

for i in range(n_ker):
    # minimum com energy added to ker
    anion_speed_ker = np.sqrt(2 * (ker[i] + ecm_min)
                              * const.eV_J / (anion_mass * (anion_mass / cation_mass + 1)))
    cation_speed_ker = anion_mass * anion_speed_ker / cation_mass
    anion_ke_ker = func.speed2ke(anion_speed_ker, anion_mass)
    cation_ke_ker = func.speed2ke(cation_speed_ker, cation_mass)
    # print(anion_ke_ker, cation_ke_ker, anion_ke_ker/cation_ke_ker)

    num_ions = int(100000 * br_ratio[i])

    for n in range(num_ions):
        pos_pu2 = random.uniform(dts_to_imd[0], dts_to_imd[dt_first])
        pos_dt_cent = random.uniform(dts_to_imd[dt_first], dts_to_imd[dt_last + 1])
        pos_pu1 = random.uniform(dts_to_imd[dt_last + 1], dts_to_imd[9])

        z_frac = random.gauss(0, 1)
        x_frac = random.gauss(0, 1)
        y_frac = random.gauss(0, 1)

        sum_frac = np.sqrt(z_frac ** 2 + x_frac ** 2 + y_frac ** 2)

        anion_speed_ax = anion_speed_ker * z_frac / sum_frac
        # anion_speed_ax = random.uniform(-anion_speed_ker, anion_speed_ker)
        anion_speed_rad = np.sqrt(anion_speed_ker ** 2 - anion_speed_ax ** 2)

        cation_speed_ax = - anion_mass * anion_speed_ax / cation_mass
        cation_speed_rad = - anion_mass * anion_speed_rad / cation_mass

        ion_pos_dt_cent.append(pos_dt_cent)
        ion_pos_pu2.append(pos_pu2)
        ion_pos_pu1.append(pos_pu1)

        anion_speed_axial.append(anion_speed_ax)
        anion_speed_radial.append(anion_speed_rad)
        cation_speed_axial.append(cation_speed_ax)
        cation_speed_radial.append(cation_speed_rad)

# for MN occurs in dt_center
# cal ion speed before MN
anion_ke_dt_cent = anion_ke - dt_cent_volt * anion_charge
cation_ke_dt_cent = cation_ke - dt_cent_volt * cation_charge
# print(anion_ke_dt_cent)
anion_speed_dt_cent = func.ke2speed(anion_ke_dt_cent, anion_mass)
cation_speed_dt_cent = func.ke2speed(cation_ke_dt_cent, cation_mass)
# cal ion speed after MN
anion_speed_fin_dt_cent = anion_speed_dt_cent + anion_speed_axial
cation_speed_fin_dt_cent = cation_speed_dt_cent + cation_speed_axial

# the slower speed will come to detector later
late_ion_speed = []
for i in range(len(anion_speed_fin_dt_cent)):
    min_sp = min(anion_speed_fin_dt_cent[i], cation_speed_fin_dt_cent[i])
    late_ion_speed.append(min_sp)

# cal time separation, position separation
anion_tof_dt_cent = ion_pos_dt_cent / anion_speed_fin_dt_cent
cation_tof_dt_cent = ion_pos_dt_cent / cation_speed_fin_dt_cent
time_separation_dt_cent = abs(anion_tof_dt_cent - cation_tof_dt_cent)

anion_distance_radial = anion_tof_dt_cent * anion_speed_radial
cation_distance_radial = cation_tof_dt_cent * cation_speed_radial
distance_separation_dt_cent = abs(anion_distance_radial) + abs(cation_distance_radial)

r_3d = np.sqrt(distance_separation_dt_cent ** 2 + (late_ion_speed * time_separation_dt_cent) ** 2)

# for MN occurs in dt_near
anion_speed_fin_pu2 = anion_speed  # be careful that at such high cm energy ker is not sure
cation_speed_fin_pu2 = cation_speed  # + cation_speed_axial

# cal time separation, position separation
anion_tof_dt_near = ion_pos_pu2 / anion_speed_fin_pu2
cation_tof_dt_near = ion_pos_pu2 / cation_speed_fin_pu2
time_separation_pu2 = abs(anion_tof_dt_near - cation_tof_dt_near)

anion_distance_radial = anion_tof_dt_near * anion_speed_radial
cation_distance_radial = cation_tof_dt_near * cation_speed_radial
distance_separation_pu2 = abs(anion_distance_radial) + abs(cation_distance_radial)

# for MN occurs in dt_far
# cal time separation, position separation
anion_tof_pu1 = ion_pos_pu1 / anion_speed_fin_pu2
cation_tof_pu1 = ion_pos_pu1 / cation_speed_fin_pu2
time_separation_pu1 = abs(anion_tof_pu1 - cation_tof_pu1)

anion_distance_radial = anion_tof_pu1 * anion_speed_radial
cation_distance_radial = cation_tof_pu1 * cation_speed_radial
distance_separation_pu1 = abs(anion_distance_radial) + abs(cation_distance_radial)

# save r3d_data for r3d
r3d_data = np.array([time_separation_dt_cent * 1.0e9, time_separation_pu2 * 1.0e9, time_separation_pu1 * 1.0e9,
                     distance_separation_dt_cent * 1.0e3, r_3d * 1.0e3]).T

r3d_header = '\t'.join(['delta_t (ns)', 'delta_t (ns)', 'delta_t (ns)', 'r (mm)', 'r3d (mm)'])
fmt_r3d = ['%.4f'] * 5

fig1, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 12))
fig1.suptitle(''.join(['HE: ', str(he_extra_volt/1.0e3), ' kV,  LE: ', str(le_extra_volt / 1.0e3),
                       ' kV,  DT ', str(dt_first), '-', str(dt_last), ': ', str(int(dt_cent_volt)), ' V']), fontsize=14)
fig1.subplots_adjust(left=0.15, bottom=0.1, right=0.95, top=0.95, wspace=0.1, hspace=0.25)

ax1.hist(time_separation_dt_cent * 1.0e9, bins=num_bins, histtype='step', label='dt_bias', alpha=0.75)
ax1.hist(time_separation_pu2 * 1.0e9, bins=num_bins, histtype='step', label='pu2', alpha=0.75)
ax1.hist(time_separation_pu1 * 1.0e9, bins=num_bins, histtype='step', label='pu1', alpha=0.75)
ax1.set_xlim(0, 500)
ax1.set_xlabel('t (ns)')
ax1.set_ylim(0, )
ax1.set_ylabel('Yield (arb. unit)')
ax1.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax1.xaxis.set_minor_locator(ticker.MultipleLocator(10))
ax1.legend(loc='best', frameon=False)

ax2.hist(distance_separation_dt_cent * 1.0e3, bins=num_bins, histtype='step', label='dt_bias', alpha=0.75)
ax2.hist(distance_separation_pu2 * 1.0e3, bins=num_bins, histtype='step', label='pu2', alpha=0.75)
ax2.hist(distance_separation_pu1 * 1.0e3, bins=num_bins, histtype='step', label='pu1', alpha=0.75)
ax2.set_xlim(0, 70)
ax2.set_xlabel('r (mm)')
ax2.set_ylim(0, )
ax2.set_ylabel('Yield (arb. unit)')
ax2.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax2.xaxis.set_minor_locator(ticker.MultipleLocator(2))
ax2.legend(loc='best', frameon=False)

ax3.hist(r_3d * 1.0e3, bins=num_bins, histtype='step', label='r3d', alpha=0.75)
ax3.set_xlim(0, 70)
ax3.set_xlabel('r3d (mm)')
ax3.set_ylim(0, )
ax3.set_ylabel('Yield (arb. unit)')
ax3.xaxis.set_major_locator(ticker.MultipleLocator(10))
ax3.xaxis.set_minor_locator(ticker.MultipleLocator(2))
ax3.legend(loc='best', frameon=False)

# save file or, not
if save2fig is True:
    fig1.savefig(''.join(['Simu_MN_DCT_DT', str(dt_first), '-', str(dt_last), '.png']), dpi=150)
    np.savetxt('_simu_r3d.txt', r_3d * 1.0e3, fmt=['%.4f'], delimiter='\t')
else:
    plt.show()
