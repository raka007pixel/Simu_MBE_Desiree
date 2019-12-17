#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# |m|i|n|g|c|h|a|o|j|i| @ 2019-10-09 16:04|
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# global functions

import const
import numpy as np
import scipy as sp


# linear function
def line(x, k, c):
    return k * x + c


# sine function
def sine(t, amp, w, phase, c):
    return amp * np.sin(w * t + phase) + c


# cosine function
def cosine(t, amp, w, phase, c):
    return amp * np.cos(w * t + phase) + c


# exponential function
def expo(x, a, b, c):
    return a * np.exp(-b * x) + c


# power-law function
def power(x, amp, f, h):
    return amp * x ** f + h


# gaussian function
def gaussian(x, amp, ctr, sig):
    return amp * np.exp(-np.power((x - ctr) / sig, 2.) / 2)


# cal ion speed from kinetic energy, eV
def ke2speed(ke_eV, mass_kg):
    return np.sqrt(2 * ke_eV * const.eV_J / mass_kg)


# cal kinetic energy from ion speed, m/s
def speed2ke(ion_speed, mass_kg):
    return 0.5 * mass_kg * ion_speed ** 2 / const.eV_J


# cal photon energy from wavelength
def lambda2ev(lambda_nm):
    return const.planck_eVs * const.speed_of_light * 1.0e9 / lambda_nm


# cm-1 to ev
def pcm2ev(freq):
    return const.planck_eVs * const.speed_of_light * 1.0e2 * freq


# cal wavelength, nm from photon energy, eV
def ev2lambda(energy_eV):
    return const.planck_eVs * const.speed_of_light * 1.0e9 / energy_eV
