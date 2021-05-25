# ENGSCI233: Lab - Sampled Data
# sdlab_earthquakes.py

# PURPOSE:
# To INVESTIGATE a dataset using interpolation and integration methods.

# PREPARATION:
# Complete the activities in sdlab_practice.py.

# SUBMISSION:
# - YOU MUST submit a plot of NET MASS CHANGE as a function of time (sdlab_earthquakes.png)
# - You MUST submit this file to complete the lab.


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


