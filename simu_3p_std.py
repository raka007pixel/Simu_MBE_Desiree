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

    # ker taken by D:
    a = - (neut_amu[1] + neut_amu[2]) * (neut_amu[1] + neut_amu[2])

    b = 2 * (neut_amu[2] * neut_amu[2] * (1 - e1) + neut_amu[1] * neut_amu[2] * (1 - e1)
             + e1 * (neut_amu[0] * neut_amu[1] - neut_amu[0] * neut_amu[2]))

    c = neut_amu[2] * neut_amu[2] * (-1 + 2 * e1) - e1 * e1 * (neut_amu[0] + neut_amu[2]) \
        * (neut_amu[0] + neut_amu[2]) + 2 * e1 * neut_amu[0] * neut_amu[2]

    e2_min = (-b - np.sqrt(b * b - 4 * a * c)) / (2 * a)
    e2_max = (-b + np.sqrt(b * b - 4 * a * c)) / (2 * a)

    e2 = random.uniform(e2_min, e2_max)
    p2_ker.append(e2)
    p3_ker.append(1 - e1 - e2)
    # print(e2_min, e2_max, e1, e2)

    # velocity vectors (see Vitali's thesis, p42)
    _v1x = np.sqrt(2 * e1 * ker[idk] * const.eV_J / neut_kg[0])

    _v2x = (neut_kg[2] * (1 - e1 - e2) - e1 * neut_kg[0] - e2 * neut_kg[1]) / neut_kg[1] * np.sqrt(
        ker[idk] * const.eV_J / (2 * e1 * neut_kg[0]))

    _vy = np.sqrt(
        2 * e2 * neut_kg[1] - (neut_kg[2] - e1 * (neut_kg[0] + neut_kg[2]) - e2 * (neut_kg[1] + neut_kg[2])) ** 2
        / (2 * e1 * neut_kg[0])) * np.sqrt(ker[idk] * const.eV_J) / neut_kg[1]

    velo_vct = [[_v1x, 0, 0],
                [_v2x, _vy, 0],
                [-(neut_kg[0] * _v1x + neut_kg[1] * _v2x) / neut_kg[2], - _vy * neut_kg[1] / neut_kg[2], 0]]
    # print(velo_vct)

    # create transformation matrix
    psi = random.uniform(0, 2 * np.pi)
    phi = random.uniform(0, 2 * np.pi)

    cos_theta = random.uniform(-1, 1)
    rd = -1 if random.uniform(-1, 1) < 0 else 1
    sin_theta = np.sqrt(1 - cos_theta ** 2) * rd
    # print(psi, phi, cos_theta, sin_theta)

    R_z1 = [[np.cos(phi), np.sin(phi), 0],
            [-np.sin(phi), np.cos(phi), 0],
            [0, 0, 1]]

    R_y = [[cos_theta, 0, -sin_theta],
           [0, 1, 0],
           [sin_theta, 0, cos_theta]]

    R_z2 = [[np.cos(psi), np.sin(psi), 0],
            [-np.sin(psi), np.cos(psi), 0],
            [0, 0, 1]]

    # final velocity
    velo_fin = np.matmul(velo_vct, np.matmul(R_z2, np.matmul(R_y, R_z1)))
    # print(velo_fin)

    # tof of three particles
    tof = [pos_dt_bias / (av_ion_speed + velo_fin[i, 2]) for i in range(n_neut)]
    for i in range(n_neut):
        tof_array[i].append(tof[i])
    # print(tof, tof_array)

    dt_array[0].append(tof[1] - tof[0])  # D - He
    dt_array[1].append(tof[2] - tof[0])  # H - He
    dt_array[2].append(tof[2] - tof[1])  # H - D

    # position on imd, unit: m
    pos_imd = [[velo_fin[j, 0] * tof[j], velo_fin[j, 1] * tof[j]] for j in range(n_neut)]
    # print(pos_imd)

    # calc distances
    d1 = np.sqrt((pos_imd[1][0] - pos_imd[0][0])**2 + (pos_imd[1][1] - pos_imd[0][1])**2)
    d2 = np.sqrt((pos_imd[2][0] - pos_imd[0][0])**2 + (pos_imd[2][1] - pos_imd[0][1])**2)
    d3 = np.sqrt((pos_imd[2][0] - pos_imd[1][0])**2 + (pos_imd[2][1] - pos_imd[1][1])**2)

    dist_array[0].append(d1)
    dist_array[1].append(d2)
    dist_array[2].append(d3)

    # calc area
    s = 0.5 * sum([d1, d2, d3])
    area = np.sqrt(s*(s-d1)*(s-d2)*(s-d3))*1.0e6    # unit: mm^2
    area_imd.append(area)

    # img array
    for j in range(n_neut):
        img_array[j][0].append(pos_imd[j][0])
        img_array[j][1].append(pos_imd[j][1])
    # print(pos_imd, img_array)

exp_diss = []
# check if particles on imd and arrival time diff in camac range
if exp_mode is True:
    imd_radius = img_size * 0.5 * 1.0e-3
    for k in range(num_ions):
        # check if particle on imd
        p1_r0 = np.sqrt(img_array[0][0][k] ** 2 + img_array[0][1][k] ** 2)
        p2_r0 = np.sqrt(img_array[1][0][k] ** 2 + img_array[1][1][k] ** 2)
        p3_r0 = np.sqrt(img_array[2][0][k] ** 2 + img_array[2][1][k] ** 2)

        if np.any([p1_r0 >= imd_radius, p2_r0 >= imd_radius, p3_r0 >= imd_radius,
                   abs(dt_array[0][k]) >= tac_convert[1] * 1.0e-9,
                   abs(dt_array[1][k]) >= tac_convert[1] * 1.0e-9,
                   abs(dt_array[2][k]) >= tac_convert[1] * 1.0e-9]):
            exp_diss.append(k)

    for i in range(n_neut):
        tof_array[i] = [tof_array[i][k] for k in range(num_ions) if k not in exp_diss]
        dt_array[i] = [dt_array[i][k] for k in range(num_ions) if k not in exp_diss]

        dist_array[i] = [dist_array[i][k] for k in range(num_ions) if k not in exp_diss]

        img_array[i][0] = [img_array[i][0][k] for k in range(num_ions) if k not in exp_diss]
        img_array[i][1] = [img_array[i][1][k] for k in range(num_ions) if k not in exp_diss]

    area_imd = [area_imd[k] for k in range(num_ions) if k not in exp_diss]

num_diss = len(exp_diss)

p1_x_pix = np.array(img_array[0][0]) * 1.0e3 / img_resolution + 0.5 * img_pixel
p2_x_pix = np.array(img_array[1][0]) * 1.0e3 / img_resolution + 0.5 * img_pixel
p3_x_pix = np.array(img_array[2][0]) * 1.0e3 / img_resolution + 0.5 * img_pixel

p1_y_pix = np.array(img_array[0][1]) * 1.0e3 / img_resolution + 0.5 * img_pixel
p2_y_pix = np.array(img_array[1][1]) * 1.0e3 / img_resolution + 0.5 * img_pixel
p3_y_pix = np.array(img_array[2][1]) * 1.0e3 / img_resolution + 0.5 * img_pixel

img_x_pix = np.concatenate((p1_x_pix, p2_x_pix, p3_x_pix), axis=0)
img_y_pix = np.concatenate((p1_y_pix, p2_y_pix, p3_y_pix), axis=0)

# list to store calculated ker from mass perm
calc_ker_perm = []

# re-build ker from simulation
num_perm_np = int(perm(n_neut, 3))
perm_neuts = list(permutations(neut_kg, 3))
# print(num_perm_np, perm_neuts)

# use permutation of neutral mass to calc ker
for k in range(num_ions - num_diss):
    ker_perm, ker_comp = [], []
    for n in range(num_perm_np):
        m1 = perm_neuts[n][0]
        m2 = perm_neuts[n][1]
        m3 = perm_neuts[n][2]
        mt = sum([m1, m2, m3])

        # calc velocity in z (axial) direction
        v1a = av2_L * (m2 * dt_array[0][k] + m3 * dt_array[1][k]) / (m3 + m2 - m1)
        v2a = av2_L * dt_array[0][k] + v1a
        v3a = av2_L * dt_array[1][k] + v1a
        # print(m1, m2, m3, mt, v1a, v2a, v3a)

        # calc velocity and G-center in x direction
        x0 = sum([m1 * img_array[0][0][k], m2 * img_array[1][0][k], m3 * img_array[2][0][k]]) / mt
        v1x = (img_array[0][0][k] - x0) * av_L
        v2x = (img_array[0][0][k] - x0) * av_L
        v3x = (img_array[0][0][k] - x0) * av_L
        # print(x0, v1x, v2x, v3x)

        # calc velocity and G-center in y direction
        y0 = sum([m1 * img_array[0][1][k], m2 * img_array[1][1][k], m3 * img_array[2][1][k]]) / mt
        v1y = (img_array[0][1][k] - y0) * av_L
        v2y = (img_array[1][1][k] - y0) * av_L
        v3y = (img_array[2][1][k] - y0) * av_L
        # print(y0, v1y, v2y, v3y)

        if exp_mode is True:
            dc = np.sqrt((x0 - img_slt_circ[0]*1.0e-3)**2 + (y0-img_slt_circ[1]*1.0e-3)**2)
            if dc <= img_slt_circ[2]*1.0e-3:
                sum_ker = 0.5 * sum([m1 * (v1x * v1x + v1y * v1y + v1a * v1a),
                                     m2 * (v2x * v2x + v2y * v2y + v2a * v2a),
                                     m3 * (v3x * v3x + v3y * v3y + v3a * v3a)]) / const.eV_J
                # print(sum_ker)
                ker_perm.append(sum_ker)
            else:
                ker_perm.append(img_pixel)    # add a large ker to keep the queue
    ker_comp = [abs(ep - ker[0]) for ep in ker_perm]
    ker_comp.extend([abs(ep - ker[1]) for ep in ker_perm])
    min_comp = np.min(ker_comp)
    min_idx = ker_comp.index(min_comp) % num_perm_np
    # print(ker_perm, min_comp, ker_perm[min_idx])
    calc_ker_perm.append(ker_perm[min_idx])

len_calc_ker = len(calc_ker_perm)

calc_ker_slt = [item for item in calc_ker_perm if 0 <= item <= 5.0]
len_ker_slt = len(calc_ker_slt)

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
