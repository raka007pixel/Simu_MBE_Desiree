#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# author: MingChao Ji
# email: mingchao.ji@fysik.su.se
# date created: 2020-01-15 13:25:52
# last modified: 2020-03-02 14:20:47
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from scipy.special import radian

import const
import func
import numpy as np
# import pandas as pd
import periodictable as pt


# for the reaction H+ + H- --> H + H, unit: eV, ker[2], n = 3 reaction, dominant at low cm energy
ker = [0.94]
br_ratio = [1]  # branching ratio, 0.48 stands for 48 %
n_ker = int(len(ker))

# num of bins for hist, bin width for density
num_bins, bin_width = 500, 0.01

# voltage range around zero center of mass energy
dt_volt_slt = [0, 400]

# ------------------------------------------Ion source settings----------------------------------------------
# low energy platform (le)
le_extra_volt = -18.60e3  # unit: V
anion_mass = (pt.H[1].mass + pt.C[13].mass * 4)* const.amu_kg + const.electron_kg
anion_charge = -1  # unit: e
anion_ke = le_extra_volt * anion_charge  # unit: eV
# fwhm of the gaussian-like KE distribution = anion_ke * anion_ke_spread
anion_ke_spread = 0.002
anion_ke_sig = anion_ke * anion_ke_spread / 2.3548

# high energy platform (he)
he_extra_volt = 6.0e3 # unit: V
cation_mass = pt.O[16].mass * const.amu_kg - const.electron_kg  # unit: kg
# print(pt.H[1].mass, pt.He[4].mass)
cation_charge = 1  # unit: e
cation_ke = he_extra_volt * cation_charge
# fwhm of the gaussian-like KE distribution = cation_ke * cation_ke_spread
cation_ke_spread = 0.002
cation_ke_sig = cation_ke * cation_ke_spread / 2.3548

# calc the reduced mass
reduced_mass = anion_mass * cation_mass / (anion_mass + cation_mass)
# -----------------------------------------------------------------------------------------------------------


# -----------------------------------------Desiree Rings settings--------------------------------------------
# geometry of desiree setup
circ_ring_s, circ_ring_a = 8.68, 8.71  # unit: m, hereafter the same
le_to_lebc, he_to_hebc = 3.25, 6.71
lebc_to_s, lebc_to_a = 7.24, 9.39
hebc_to_s, hebc_to_a = 9.93, 7.80

# settings on drift tube
# give the first and last tube numbers (from 1 to 7), if only one used, write one number
drift_tube_bias = [3, 5]

ms_ctr_to_imd = 1.690  # center, unit: m
# distance from imd to: 10def center, dt1-7 start, dt7 end, 10 def center, interval shared by two neighbouring tubes
dts_to_imd = [1.166, 1.410, 1.4875, 1.5665,
              1.6505, 1.7295, 1.8135, 1.8925, 1.970, 2.214]
# diff = [dist_to_imd[i+1] - dist_to_imd[i] for i in range(len(dist_to_imd)-1)]
dt_first, dt_last = min(drift_tube_bias), max(drift_tube_bias)
dt_ctr_to_imd = (dts_to_imd[dt_first] + dts_to_imd[dt_last + 1]) * 0.5

# beam merge condition, assumed to be the best, 0.1 degree
beam_merge_deg = 0.0
beam_merge_rad = radian(beam_merge_deg, 0, 0)

# to calc the minimum center of mass energy
drift_tube_volt = np.linspace(-2000, 2000, 4001)  # unit: v

anion_ke_tube = anion_ke - drift_tube_volt * anion_charge
cation_ke_tube = cation_ke - drift_tube_volt * cation_charge

# eV_J convert can be omitted on both sides
com_energy = reduced_mass * (anion_ke_tube / anion_mass + cation_ke_tube / cation_mass - 2 * np.sqrt(
    anion_ke_tube * cation_ke_tube / (anion_mass * cation_mass)) * np.cos(beam_merge_rad))
# print(drift_tube_volt, com_energy)

# find the min value of ecm_energy and the corresponding drift_tube_volt.
ecm_min = np.min(com_energy)
ecm_min_idx = com_energy.tolist().index(ecm_min)

# to calc the signal at ~ 0 eV center of mass energy
dt_cent_volt = drift_tube_volt[ecm_min_idx]  # unit: V
# -----------------------------------------------------------------------------------------------------------


# ---------------------------------Detector & Data Acquisition Settings--------------------------------------
# Imaging Detector
img_pixel = 352  # pix, calibrated from exp img
img_size = 75.0  # mm
img_resolution = img_size / img_pixel  # mm/pix
# print(img_resolution)

tick_freq = 2083333

# imd spot size, intensity discrimination
# ref value to start: spot_size_diss = [0, 2000], spot_intensity_diss = [0, 1e5]
spot_size_diss = [0, 8000]
spot_intensity_diss = [0, 50000]

max_img_spots, img_count_diss = 6, 1  # max num of spots, plot discrimination

# calibrated from 1 camac_time 1 spot event
# ref to the y values of spots
imd_strip_center = [353, 331.37, 309.65, 287.94, 266.71, 245.11, 221.44, 201.82, 179.73, 158.07, 137.43,
                    115.93, 93.26, 72.22, 50.33, 30.03]
# print(imd_strip_center[::-1])
imd_strip_hw = 20

# data acquisition system
# convert voltage to time, time = 500ns/10V + 95 ns
tac_convert = [50, 180]

# camac basic settings
camac_strips = 16
camac_time_diss = 20  # unit: ns

# used in the camac time calibration
# same discrimination as camac_time_diss, but without unit.
camac_value_diss = 180

# input time = 16 * n + 25  unit: ns, might not be true, use the measured one for calibration
camac_convert = [16, 25]
# camac time measured in exp during the calibration
camac_time = [18.96, 40.66, 61.17, 80.39, 98.04, 116.92, 135.43, 147.79, 157.17, 167.69, 177.80, 188.00, 197.00,
              207.41, 217.23, 227.08]

# match ratio diss when combining camaco and img files
match_ratio_diss = 0.6

# create a list of offset to correct the frame id
max_id_offset = 5
# create a list of offset numbers for img
id_offset = [[0] if i == 0 else [-i, i] for i in range(max_id_offset + 1)]
frame_id_offset = sum(id_offset, [])
# print(frame_id_offset)

max_cycle_offset = 20
# create a list of offset numbers for storage cycle num
cycle_num_offset = [j for j in range(max_cycle_offset + 1)]
# print(cycle_num_offset)
# -----------------------------------------------------------------------------------------------------------


# ---------------------------------------Simple Calcs & Print2Screen-----------------------------------------
print2scr = False

# simple calcs
anion_speed = func.ke2speed(anion_ke, anion_mass)  # unit: m/s
cation_speed = func.ke2speed(cation_ke, cation_mass)  # unit: m/s

tof_le_to_s = (le_to_lebc + lebc_to_s) * 1.0e6 / anion_speed  # us
tof_le_to_a = (le_to_lebc + lebc_to_a) * 1.0e6 / anion_speed  # us
tof_he_to_s = (he_to_hebc + hebc_to_s) * 1.0e6 / cation_speed  # us
tof_he_to_a = (he_to_hebc + hebc_to_a) * 1.0e6 / cation_speed  # us

rev_time_s = circ_ring_s * 1e6 / cation_speed  # unit: microsec
rev_time_a = circ_ring_a * 1e6 / anion_speed  # unit: microsec

neutral_speed = anion_speed
imd_gate_delay = ms_ctr_to_imd * 1.0e6 / neutral_speed  # unit: us

if print2scr is True:
    print('le ion speed (m/s): {0}'
          '\nhe ion speed (m/s): {1}\n'
          '\ntof le to s-ring (us): {2}'
          '\ntof he to a-ring (us): {3}\n'
          '\ns-ring revolution period (us): {4}'
          '\na-ring revolution period (us): {5}\n'
          '\nid signal delay (us): {6}'
          '\n'.format(str(anion_speed), str(cation_speed), str(tof_le_to_s), str(tof_he_to_a), str(rev_time_s), str(rev_time_a), str(imd_gate_delay)))
#------------------------------------------------------------------------------------------------------------
