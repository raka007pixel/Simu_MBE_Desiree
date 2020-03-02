#!/usr/local/bin/python3
# -*- coding: utf-8 -*-
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
# author: MingChao Ji
# email: mingchao.ji@fysik.su.se
# date created: 2020-01-15 13:26:29
# last modified: 2020-03-02 14:21:01
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+

from desiree import *

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
mpl.rc('pdf', fonttype=42)
mpl.rc('axes.spines', left=True, top=True, right=True, bottom=True)
mpl.rcParams.update({'font.family': 'serif', 'font.serif': 'Times', 'mathtext.fontset': 'cm', 'font.size': 14})

save2fig = True

# beam speed in the drift tube
anion_speed_tube = func.ke2speed(anion_ke_tube, anion_mass)
cation_speed_tube = func.ke2speed(cation_ke_tube, cation_mass)

speed_diff = abs(anion_speed_tube - cation_speed_tube)


# if len(dt_volt) > 0:
#     # two dicts
#     volt_ecm = dict(zip(drift_tube_volt, com_energy))
#     volt_spd_diff = dict(zip(drift_tube_volt, speed_diff))

#     # find corresponding ecm and speed diff
#     exp_ecm = [volt_ecm[v] for v in dt_volt]
#     exp_spd_diff = [volt_spd_diff[v] for v in dt_volt]

#     # create an array to save the related info
#     exp_array = np.array([dt_volt, exp_ecm, exp_spd_diff]).T
#     svf = ['%d'] + ['%.4f'] * 2

#     np.savetxt('_ecm_speed_diff.txt', exp_array, fmt=svf, delimiter='\t')


# plot
fig, ax = plt.subplots(figsize=(8, 6))
fig.suptitle(''.join([r'$\mathrm{E_{CM}}$ vs DT Volt']), fontsize=16)
fig.subplots_adjust(left=0.15, bottom=0.15, right=0.9, top=0.9, wspace=0.1, hspace=0.1)

ax.plot(drift_tube_volt, com_energy, 'r-', label=r'$\mathrm{E_{CM}}$')
ax.set_xlim(dt_volt_slt[0], dt_volt_slt[1])
ax.xaxis.set_major_locator(ticker.MultipleLocator(50))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(5))
ax.set_xlabel("Drift Tube Voltage (V)")

ax.set_ylim(0, 5.0)
ax.yaxis.set_major_locator(ticker.MultipleLocator(0.5))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(0.1))
ax.set_ylabel(r"$\mathrm{E_{CM}}$ (eV)")
ax.grid(b=None, which='major', axis='both', color='k', linestyle='-', linewidth=0.5)
ax.grid(b=None, which='minor', axis='both', color='k', linestyle='-.', linewidth=0.5)

if save2fig is True:
    fig.savefig("Ecm_vs_DT_Volt.png")
else:
    print("index number: {0}"
          "\ncorresp drift tube volt: {1}"
          "\nminimum com_energy: {2}\n"
          "\nanion tube speed: {3}"
          "\ncation tube speed: {4}\n".format(str(ecm_min_idx), str(dt_cent_volt), str(ecm_min),
                                              str(anion_speed_tube[ecm_min_idx]),
                                              str(cation_speed_tube[ecm_min_idx])))
    plt.show()

