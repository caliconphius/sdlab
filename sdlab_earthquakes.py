# ENGSCI233: Lab - Sampled Data
# sdlab_earthquakes.py

# PURPOSE:
# To INVESTIGATE a dataset using interpolation and integration methods.

# PREPARATION:
# Complete the activities in sdlab_practice.py.

# SUBMISSION:
# - YOU MUST submit a plot of NET MASS CHANGE as a function of time (sdlab_earthquakes.png)
# - You MUST submit this file to complete the lab.

import numpy as np
from scipy.integrate import trapz
from sdlab_functions import *
from numpy.linalg import solve, norm
import matplotlib.pyplot as plt

# EXERCISE: Analysis of Net Mass Changes.
#
# Earthquakes are sometimes associated with oil and gas production (taking mass out of the 
# ground) and injection (putting it back in) operations.
# 
# It has been suggested that injection of water at one particular site, which started midway through  
# 1993, has been responsible for a spate of recent earthquakes there. These are the data that were 
# plotted in the first exercise of sdlab_practice.py. The operator of the field has claimed they 
# cannot be responsible, because injection had been ongoing for almost 10 years before any earthquakes
# occurred.
#
# It has been proposed the earthquakes may be related to NET MASS CHANGES in the field. Therefore,
# it is necessary to understand how this quantity has evolved over time.
#
# Although data from the two production wells (mass extractors) - PW1 and PW2 - are reported regularly,
# data reporting from the injection well, IW1, is more irregular. In addition, the operator only 
# reports MASS RATES, not CUMULATIVE production or injection MASS.
#
# TO solve this problem, you will need to use both INTERPOLATION and INTEGRATION.

# get data from dat files
tm,pm1 = np.genfromtxt('PW1.dat',delimiter=',',skip_header=1).T
tq,pm2 = np.genfromtxt('PW2.dat',delimiter=',',skip_header=1).T
ty,iy = np.genfromtxt('IW1.dat',delimiter=',',skip_header=1).T

# get interpolation for pm2
pm2Matrix = spline_coefficient_matrix(tq)
pm2RHS = spline_rhs(tq,pm2)
akpm2 = solve(pm2Matrix,pm2RHS)

pm2_interpo = spline_interpolate(tm, tq, akpm2)

# get interpolation for iy
iyMatrix = spline_coefficient_matrix(ty)
iyRHS = spline_rhs(ty,iy)
akiy = solve(iyMatrix,iyRHS)

iy_interpo = spline_interpolate(tm, ty, akiy)

# total mass flow rate
mass_inflow_rate = iy_interpo - pm2_interpo - pm1

# integrate mass flow 
seconds_in_a_year = 31536000 # to convert kg/s to kg/year
k_to_G_constant = 10**-6 # convert kg to Gg
net_mass_change = [seconds_in_a_year * trapz(mass_inflow_rate[:i], tm[:i]) * k_to_G_constant for i in range(len(tm))]


# create a 
f, ax1 = plt.subplots(nrows=1,ncols=1)
ax2 = ax1.twinx()

ax1.plot(tm, net_mass_change, 'k-', label='Net mass change')


# put in arrows showing when each earthquake happened and each of their magnitudes (arrow height represents magnitude)
ax2.arrow(2003.5, 3.5, 0., -1, length_includes_head=True, head_width=0.3, head_length=0.4, color = 'b',label = '2003 Magnitude 3.5 earthquake')
ax2.text(2003., 3.2, 'M 3.5', ha= 'right', va = 'center', size=10, color = 'b')
ax2.arrow(2004.5, 4, 0., -1, length_includes_head=True, head_width=0.4, head_length=0.3, color = 'r',label = '2004 Magnitude 4.0 earthquake')
ax2.text(2003.2, 4.1, 'M 4.0', ha= 'center', va = 'bottom', size=10, color = 'r')
ax2.arrow(2005., 4.3, 0., -1, length_includes_head=True, head_width=0.3, head_length=0.4, color = 'g',label = '2005 Magnitude 4.3 earthquake')
ax2.text(2005.5, 3.2, 'M 4.3', ha= 'left', va = 'center', size=10, color = 'g')

# add a legend and axes/plot labels 

ax1.legend(loc=2, fontsize = 'small')
ax2.legend(loc=3, fontsize = 'small')
ax2.set_ylim([0,10])
ax1.set_ylabel('net change in field mass flow [Gg]')
ax2.set_ylabel('Earthquake Magnitude')
ax1.set_xlabel('time [yr]')
ax2.set_title('Total Mass change and major earthquakes in the field over time')


# save figure
plt.savefig('sdlab_earthquakes.png',dpi=300)

# TO DO:
# - In sdlab_functions.py, COMPLETE the functions SPLINE_COEFFICIENT_MATRIX, SPLINE_RHS, and
#   SPLINE_INTERPOLATE.
# - Write a Newton-Cotes integration function OR find and use a built-in Python function.
# - Produce a plot of NET MASS CHANGE as a function of time.
# - ANSWER the questions in sdlab_questions.txt


# HINTS:
# - To add or difference two quantities, they should be measured or interpolated at the same time.
# - You should consider a sensible strategy for the event that an interpolation point lies outside 
#   the range of data (extrapolation).
# - MASS RATE is a derivative (per unit time) and CUMULATIVE MASS is its integral.
# - You will be assessed on the ability of your figure to convey information. Things that help:
#    o Sensible labels (and units) for the axes.
#    o A legend.
#    o Sensible use of lines versus markers.
#    o Juxtaposition of information in a way that highlights correlation.
# - The plotting exercise in PRACTICE TASK ONE is relevant here.


