# ENGSCI233: Lab - Sampled Data
# sdlab_functions.py

# PURPOSE:
# To IMPLEMENT cubic spline interpolation.

# PREPARATION:
# Notebook sampling.ipynb, ESPECIALLY Section 1.3.1 theory of cubic splines.

# SUBMISSION:
# - YOU MUST submit this file to complete the lab. 
# - DO NOT change the file name.

# TO DO:
# - COMPLETE the functions spline_coefficient_matrix(), spline_rhs() and spline_interpolation().
# - COMPLETE the docstrings for each of these functions.
# - TEST each method is working correctly by passing the asserts in sdlab_practice.py.
# - DO NOT modify the other functions.

import numpy as np

# **this function is incomplete**
#					 ----------
def spline_coefficient_matrix(xi):	
	''' 
	Creates a matrix detailing the equation satisfied by the coefficients of the spline interpolating points xi

	inputs:
	-------
	xi : np array
		List of the points which the spline is interpolating
	Outputs:
	-------
	A : np array
		Matrix of coefficients needed to solve the spline equation
	'''
	
	# create an array of zeros with the correct dimensions
	#   **what are the correct dimensions? how to determine this from xi?**
	#   **use np.zeros() to create the array**
	npoints = len(xi)
	dim = 4*(npoints-1)
	A = np.zeros([dim,dim])
	delta_x = []
	[delta_x.append(xi[j+1]-xi[j]) for j in range(npoints-1)]
	
	
	# Loop over the subintervals, add matrix coefficients for equations:
	# - polynomial passes through lefthand point of subinterval
	# - polynomial passes through righthand point of subinterval
	#   **how many subintervals should there be (in terms of length of xi)?**
	#   **how to modify loop index so it jumps along a row in increments of 4?**
	#   **how to define width of the subinterval in terms of indices of xi?**
	#   **what values go into matrix A and how do they relate to subinterval width?**
	
	# loop sets spline as compatible with data points
	for i in range(0,2*(npoints-1),2):
		# set polynomial equal to yi at left hand side
		A[i,2*i] = 1
		# set polynomial equal to yi+1 on right hand side
		for j in range(4):
			A[i+1,2*i+j] = delta_x[int(i/2)]**j
		
	
		
		
	# Loop over neighbouring subintervals, add matrix coefficients for equations:
	# - polynomial gradient continuous at shared point
	# - polynomial second derivative continuous at shared point
	#   **how many shared points should there be (in terms of length of xi)?**
	#   **what values go into matrix A and how do they relate to subinterval width?**

		# loop sets spline as compatible with data points
	for i in range(0,npoints-2):
		# set first derivative compatibility
		
		A[2*(npoints-1) + i,(4*i+1):(4*i+4)] = [j*delta_x[i]**(j-1)	for j in range(1,4)]	
		A[2*(npoints-1) + i, 4*i+5] = -1		

	for i in range(0,npoints-2):
		# set second derivative compatibility
		A[3*(npoints)-4 + i,(4*i+2):(4*i+4)] = [j*(j-1)*delta_x[i]**(j-2) for j in range(2,4)]			
		A[3*(npoints)-4 + i, 4*i+6] = -2	

	
	# For the beginning and end points, add matrix coefficients for equations:
	# - the polynomial second derivative is zero
	
	# 2a_2^(0) = 0
	A[dim-2, 2] = 2

	# 2a_2^(n-1) + 6a_3^(n-1)delta_x[n-1] = 0
	A[dim-1, dim-2] = 2
	A[dim-1, dim-1] = 6*delta_x[-1]
	return A

# **this function is incomplete**
#					 ----------
def spline_rhs(xi, yi):
	''' creates the right-hand side vector in the spline matrix equation


		inputs
		------
		xi : np array
			list of points on the x-axis the spline is interpolating
		yi : np array
			list of points on the y-axis the spline is interpolating

		outputs
		-------
		b : np array
			rhs vector for the spline coefficient equation Ax = b
	'''
	npoints = len(yi)
	# initialise vector
	b = np.zeros(4 * (npoints-1))
	
	# generate a vector with each intermediate point's y value duplicated
	yi_prime = [yi[0]]
	[yi_prime.extend([yi[i],yi[i]]) for i in range(1,npoints-1)]
	yi_prime.append(yi[-1])

	# set the first 2(npoints - 1) elements to these values
	b[:2*(npoints-1)] = np.array(yi_prime)

	return list(b)
	
# **this function is incomplete**
#					 ----------
def spline_interpolate(xj, xi, ak):
	''' interpolates the spline through points xi at the points xj

		xj : np array 
			interpolating points
		xi : np array 
			sampled points
		ak : np array 
			spline coefficients

		output:
		-------
		yj : np.array
			interpolated y values
	
		Notes
		-----
		You may assume that the interpolation points XJ are in ascending order.
		Evaluate polynomial using polyval function DEFINED below.
	'''
	
	# Suggested strategy (you could devise another).
	# 1. Initialise FIRST subinterval (and polynomial) as CURRENT subinterval (and polynomial).
	# 2. FOR each interpolation point.
	# 3. WHILE interpolation point NOT inside CURRENT subinterval, iterate
	#    to NEXT subinterval (and polynomial).
	# 4. Evaluate CURRENT polynomial at interpolation point.
	# 5. RETURN when all interpolation points evaluated.
	
	# raise an error if any point is outside the range of the sampled data
	if (any([points < xi[0] for points in xj]) or any([points > xi[-1] for points in xj])):
		raise ValueError("A Point is out of the sampling range")

	# initialise subinterval and polynomial coefficients as well as return value
	nsamples = len(xi)
	intervals = [xi[i:i+2] for i in range(0,nsamples-1)]
	polynomials = [ak[i:i+4] for i in range(0, 4*(nsamples-1),4)]
	yj = np.array([])
	# iterated index specifying which interval we are in
	interv = 0

	# find the point x1
	for ipoint in xj:
		# change interval until it contains the interpolating point
		while not ((ipoint >= intervals[interv][0]) and (ipoint <= intervals[interv][1])):
			# iterate subinterval
			interv +=1
		yj = np.append(yj,polyval(polynomials[interv],ipoint-xi[interv]))

	return yj
	
# this function is complete
def display_matrix_equation(A,b):
	''' Prints the matrix equation Ax=b to the screen.
	
		Parameters
		----------
		A : np.array
			Matrix.
		b : np.array
			RHS vector.
			
		Notes
		-----
		This will look horrendous for anything more than two subintervals.	
	'''
	
	# problem dimension
	n = A.shape[0]
	
	# warning
	if n > 8:
		print('this will not format well...')
		
	print(' _'+' '*(9*n-1) +'_  _       _   _        _')
	gap = ' '
	for i in range(n):
		if i == n - 1:
			gap = '_'
		str = '|{}'.format(gap)
		str += ('{:+2.1e} '*n)[:-1].format(*A[i,:])
		str += '{}||{}a_{:d}^({:d})'.format(gap,gap,i%4,i//4+1)+'{}|'.format(gap)
		if i == n//2 and i%2 == 0:
			str += '='
		else:
			str += ' '
		if b:
			str += '|{}{:+2.1e}{}|'.format(gap,b[i],gap)
		print(str)
	
# this function is complete
def get_data():
	# returns a data vector used during this lab
	xi = np.array([2.5, 3.5, 4.5, 5.6, 8.6, 9.9, 13.0, 13.5])
	yi = np.array([24.7, 21.5, 21.6, 22.2, 28.2, 26.3, 41.7, 54.8])
	return xi,yi
		
# this function is complete
def ak_check():
	# returns a vector of predetermined values
	out = np.array([2.47e+01, -4.075886048665986e+00,0.,8.758860486659859e-01,2.15e+01,
		-1.448227902668027e+00,2.627658145997958e+00,-1.079430243329928e+00,2.16e+01,
		5.687976593381042e-01,-6.106325839918264e-01,5.358287012458253e-01,2.22e+01,
		1.170464160078432e+00,1.157602130119396e+00,-2.936967278262911e-01,2.82e+01,
		1.862652894849505e-01,-1.485668420317224e+00,1.677900564431842e-01,2.63e+01,
		-2.825777017172887e+00,-8.312872001888050e-01,1.079137281294699e+00,4.17e+01,
		2.313177016138269e+01,9.204689515851896e+00,-6.136459677234598e+00])
	return out
	
# this function is complete
def polyval(a,xi):
	''' Evaluates a polynomial.
		
		Parameters
		----------
		a : np.array
			Vector of polynomial coefficients.
		xi : np.array
			Points at which to evaluate polynomial.
		
		Returns
		-------
		yi : np.array
			Evaluated polynomial.
			
		Notes
		-----
		Polynomial coefficients assumed to be increasing order, i.e.,
		
		yi = Sum_(i=0)^len(a) a[i]*xi**i
		
	'''
	# initialise output at correct length
	yi = 0.*xi
	
	# loop over polynomial coefficients
	for i,ai in enumerate(a):
		yi = yi + ai*xi**i
		
	return yi