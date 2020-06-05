##  @docstring	StatisticsExercises 
# 	@file		ci_mean_var_ukn_t.sage
#	@brief		t-test for the sample mean with unknown variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes t-distribution confidence intervals for sample mean with unknown variance. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()
import sage.probability as pr

#data from Exercise 35
#sample size, sample mean, sample variance
nx=15
mx = 17250.0
stdx = 660.0

#confidence and significance levels
#symmetric сl=90% for ex.35, left-open cl=95% for ex.36
conf = 0.95
sig = 1-conf

#degrees of freedom
dof = nx-1

#sample mean variance estimate
sigma = stdx/n(sqrt(nx))

#t-distribution quantiles
T = RealDistribution('t', dof)

t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = mx - t_sym*sigma
sym_end = mx + t_sym*sigma

left_open_start = -infinity
left_open_end = mx + t_val*sigma

right_open_start = mx - t_val*sigma
right_open_end = +infinity

print 'sample size n =', nx
print 'x mean  = ', mx
print 'x stdev = ', stdx
print 'sample mean variance = ', sigma
print

print 'degrees of freedom = ', dof
print sig, '-significance t-quantiles:'
print 't(sig/2) = ', -t_sym, ' t(1-sig/2) = ', t_sym
print 't(sig)   = ', -t_val, ' t(conf)    = ', t_val
print
 
print 'sigma-squared', conf, '-confidence intervals:'
print '---------------'
print 'symmetric : (', sym_start, ', ', sym_end, ')'
print 'left open : (', left_open_start, ', ', left_open_end, ')'
print 'right open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'