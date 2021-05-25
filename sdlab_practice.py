# ENGSCI233: Lab - Sampled Data
# sdlab_practice.py

# PURPOSE:
# To INTRODUCE the basics of plotting in Python. 
# To TEST your implementation of cubic spline interpolation.

# PREPARATION:
# Notebook sampling.ipynb, ESPECIALLY Section 1.3.1 theory of Cubic splines.

# SUBMISSION:
# DO NOT submit this file. 


# imports
from matplotlib import pyplot as plt    # MATPLOTLIB is THE plotting module for Python
from sdlab_functions import*
from numpy.linalg import norm, solve
tol = 1.e-6

# Boolean variables to run each exercise
run_ex1 = False # boolean to run EXERCISE 1
run_ex2 = True # boolean to run EXERCISE 2
run_ex3 = False # boolean to run EXERCISE 3


# EXERCISE 1: Plotting sampled data (max 20 mins)
# This exercise introduces some basic plotting commands. Everything is already complete.

# TO DO:
# - RUN the code and inspect the generated plot.
# - Add COMMENTS as you understand the code.
# - Note, more context will be given to these data in the ASSESSED task.

if run_ex1:

	# File I/O commands to read in the data
	tm,pm1 = np.genfromtxt('PW1.dat',delimiter=',',skip_header=1).T
	tq,pm2 = np.genfromtxt('PW2.dat',delimiter=',',skip_header=1).T
	ty,iy = np.genfromtxt('IW1.dat',delimiter=',',skip_header=1).T

	# **Whats happening here?**
	f,ax1 = plt.subplots(nrows=1,ncols=1)
	ax2 = ax1.twinx()				# twinned plots are a powerful way to juxtapose data

	# **Whats happening here?**
	ax1.plot(tm, pm1, 'k-', label='PW1')
	ax1.plot(tq, pm2, 'b-', label='PW2')
	ax2.plot(ty, iy, 'r*', markersize = 7)

	# **Whats happening here?**
	ax1.arrow(2003.5, 15, 0., -0.75, length_includes_head=True, head_width=0.2, head_length=0.1, color = 'b')
	ax1.text(2003., 14.65, 'M 3.5', ha= 'right', va = 'center', size=10, color = 'b')
	ax1.arrow(2004.5, 15, 0., -0.75, length_includes_head=True, head_width=0.2, head_length=0.1, color = 'b')
	ax1.text(2004.5, 15.2, 'M 4.0', ha= 'center', va = 'bottom', size=10, color = 'b')
	ax1.arrow(2005., 15, 0., -0.75, length_includes_head=True, head_width=0.2, head_length=0.1, color = 'b')
	ax1.text(2005.5, 14.65, 'M 4.3', ha= 'left', va = 'center', size=10, color = 'b')

	# **Whats happening here?**
	ax2.set_ylim([0,40])
	ax1.legend(loc=2)
	ax1.set_ylim([0,20])
	ax1.set_ylabel('production rate [kg/s]')
	ax2.set_ylabel('injection rate [kg/s]')
	ax2.set_xlabel('time [yr]')
	ax2.set_title('Summary of field operations and seismicity at field X')

	# EITHER show the plot to the screen OR save a version of it to the disk
	save_figure = False
	if not save_figure:
		plt.show()
	else:
		plt.savefig('sdlab_plot.png',dpi=300)

	
# EXERCISE 2: Cubic spline interpolation
# In Section 1.3.1 of the Sampling notebook, we developed the theory around cubic
# spline interpolation. Here, you will implement it.

# TO DO:
# - Review the 'by-hand' construction of the spline coefficient matrix and RHS vector
#   for the in-class exercise.
# - COMPLETE the functions SPLINE_COEFFICIENT_MATRIX, SPLINE_RHS, and SPLINE_INTERPOLATE
#   in sdlab_functions.py
# - How do I know its CORRECT? Test you implementation by PASSING the ASSERTS below.
# - How do I DEBUG? Add rows one at a time and check the printed matrix and RHS match
#   your BY HAND example.
# - You DO NOT need to modify any of the commands below.

if run_ex2:
		# i. define the data points (class example)
	xi = np.array([1., 2., 4.])
	yi = np.array([2., 5., 4.])

		# ii. assemble the coefficient matrix
	A = spline_coefficient_matrix(xi)

		# iii. assemble the RHS vector
	b = spline_rhs(xi,yi)

		# DISPLAY A and b (comment this if its annoying)
	display_matrix_equation(A,b)
		
		# CHECK A and b
	A_soln = np.array([
	[ 1, 0, 0, 0, 0, 0, 0, 0],
	[ 1, 1, 1, 1, 0, 0, 0, 0],
	[ 0, 0, 0, 0, 1, 0, 0, 0],
	[ 0, 0, 0, 0, 1, 2, 4, 8],
	[ 0, 1, 2, 3, 0,-1, 0, 0],
	[ 0, 0, 2, 6, 0, 0,-2, 0],
	[ 0, 0, 2, 0, 0, 0, 0, 0],
	[ 0, 0, 0, 0, 0, 0, 2,12],
	])	 									# correct output of SPLINE_COEFFICIENT_MATRIX
	b_soln = np.array([2,5,5,4,0,0,0,0])	# correct output of SPLINE_RHS
	assert norm(A - A_soln) < tol 			# test agreement
	assert norm(b - b_soln) < tol

		# iv. solve matrix equations to obtain spline polynomial coefficients
	ak = solve(A,b)					# (could have used your LU factor method from last week)
	ak_soln = np.array([2, 43/12, 0, -7/12, 5, 11/6, -1.75, 7/24])	# correct coefficients
	assert norm(ak - ak_soln) < tol 								# test agreement

		# v. interpolate points
	xj = np.array([1.5, 2.5, 3.5])   
	yj = spline_interpolate(xj, xi, ak)
	yj_soln = np.array([119/32, 353/64, 307/64])	# correct output of SPLINE_INTERPOLATE
	assert norm(yj - yj_soln) < tol 				# test agreement



# EXERCISE 3 i.e. FINAL TEST:
# Run the code below to test your spline interpolation. If you pass the asserts,
# a plot will automatically be generated that should match the notebook/lab document.

if run_ex3:
	# interpolate
	xi,yi = get_data()                  			# get notebook data
	A = spline_coefficient_matrix(xi)   			# get spline matrix
	print(A)
	b = spline_rhs(xi,yi)							# get rhs vector
	ak = solve(A,b)									# get polynomial coefficients
	assert(np.sum(abs(ak - ak_check()))<1.e-6)		# test agreement

	# plotting commands to show the cubic splines and data
	f,ax = plt.subplots(1,1)
	label = 'cubic spline interpolation'
	ax.plot(xi,yi,'ko',label='data',zorder=2)
	for i, xi1, xi2 in zip(range(len(xi)-1), xi[:-1], xi[1:]):
		xv = np.linspace(xi1,xi2,101) 
		yv = polyval(ak[4*i:4*(i+1)], xv-xi1)
		ax.plot(xv,yv,'r-',label=label)
		label = None
	ax.legend(loc=2)
	ax.set_xlabel('time')
	ax.set_ylabel('temperature')		
	plt.show()


	






