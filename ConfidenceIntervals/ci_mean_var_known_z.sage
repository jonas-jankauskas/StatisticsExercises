##  @docstring	StatisticsExercises 
# 	@file		ci_means_var_known_z.sage
#	@brief		Z-test for the sample mean with known variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes sample mean Z-distribution confidence intervals and tests hypothesis about the mean when variance is known. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#data from Exercise 34
#sample size, sample mean, population variance
nx = 15
mx =  34500.0
sigma = 1320.0

#sample variance
sdev = sigma/n(sqrt(nx))

#significance and confidence levels
conf = 0.95;
sig = 1-conf;

#sigma2-distribution quantiles
Z = RealDistribution('gaussian', 1)

z_sym = Z.cum_distribution_function_inv(1-sig/2)
z_val = Z.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = mx - z_sym*sdev
sym_end = mx + z_sym*sdev

left_open_start = -infinity
left_open_end = mx + z_val*sdev

right_open_start = mx - z_val*sdev
right_open_end = +infinity

print 'x mean   = ', mx
print 'x sigma  = ', sigma
print 'sample mean variance = ', sdev

print sig, '-significance Z-quantiles:'
print 'z(sig/2) = ', -z_sym, ' z(1-sig/2) = ', z_sym
print 'z(sig)   = ', -z_val, ' z(conf)    = ', z_val
print
 
print 'sigma-squared', conf, '-confidence intervals:'
print '---------------'
print 'symmetric : (', sym_start, ', ', sym_end, ')'
print 'left open : (', left_open_start, ', ', left_open_end, ')'
print 'right open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'