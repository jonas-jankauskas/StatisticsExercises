##  @docstring	StatisticsExercises 
# 	@file		ci_diffmeans_vars_known_z.sage
#	@brief		Z-test for the difference of two sample means with known variances
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes Z-distribution confidence intervals for the difference of two sample means with known variances. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#data from Exercise 38
x = [49.7, 51.9, 49.0, 49.6, 51.2, 49.2, 49.8, 49.9, 46.5, 51.7]
y = [48.4, 46.3, 51.1, 49.2, 48.7, 49.4, 49.9, 48.3]

#significance and confidence levels
conf = 0.99
sig = 1-conf

#the variances are given
varx = 2.56
vary = 1.96

#element counts and means
nx = len(x)
ny = len(y)

mx = mean(x)
my = mean(y)

mxy = mx-my

#variance and standard dev. of X-Y
varxy = varx/nx + vary/ny
stdxy = sqrt(varxy)

#Distribution and quantiles
Z = RealDistribution('gaussian', 1)

z_sym = Z.cum_distribution_function_inv(1-sig/2)
z_val = Z.cum_distribution_function_inv(conf)

#Confidence intervals
sym_start = mxy - stdxy * z_sym
sym_end = mxy + stdxy * z_sym

left_open_start = -infinity;
left_open_end = mxy + stdxy * z_val

right_open_start = mxy - stdxy * z_val
right_open_end = +infinity

#Output
print 'n = ', nx
print 'm = ', ny
print 'x = ', x
print 'y = ', y
print

print 'x mean    = ', mx
print 'y mean    = ', my
print 'x-y mean  = ', mxy
print 'x var     = ', varx
print 'y var     = ', vary
print 'x-y var   = ', varxy
print 'x-y stdev = ', stdxy
print

print sig, '-significance Z-quantiles:'
print 'z(1-sig/2) = ', z_sym
print 'z(sig)     = ', z_val
print
 
print 'x-y mean ', conf, '-confidence intervals:'
print '---------------'
print 'symmetric: (', sym_start, ', ', sym_end, ')'
print 'left_open: (', left_open_start, ', ', left_open_end, ')'
print 'right_open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'