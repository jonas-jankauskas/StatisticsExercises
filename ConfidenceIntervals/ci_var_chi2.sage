##  @docstring	StatisticsExercises 
# 	@file		ci_var_chi2.sage
#	@brief		Sample variance chi-squared test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes confidence intervals for sample variance. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#data from Exercise 37
x = [0.74, 0.75, 0.73, 0.75, 0.74, 0.72]

#significance and confidence levels
conf = 0.95;
sig = 1-conf;

#element count
n = len(x)

#degrees of freedom
dof = n-1

#sample variance and estimate for sigma-squared (chi2)
varx = variance(x)
sigma2 = dof*varx

#chi2-distribution quantiles
chi2 = RealDistribution('chisquared', dof)

chi2_sym_left = chi2.cum_distribution_function_inv(sig/2)
chi2_sym_right = chi2.cum_distribution_function_inv(1-sig/2)

chi2_val_left = chi2.cum_distribution_function_inv(sig)
chi2_val_right = chi2.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = sigma2/chi2_sym_right
sym_end = sigma2/chi2_sym_left

left_open_start = 0.0;
left_open_end = sigma2/chi2_val_left

right_open_start = sigma2/chi2_val_right
right_open_end = +infinity

#output
print 'n=', n

print 'x=', x
print

print 'x var    = ', varx
print 'sigma2   = ', sigma2
print

print 'degrees of freedom = ', dof
print sig, '-significance chi2-quantiles:'
print 'chi2(sig/2) = ', chi2_sym_left, ' chi2(1-sig/2) = ', chi2_sym_right
print 'chi2(sig)   = ', chi2_val_left, ' chi2(conf)    = ', chi2_val_right
print
 
print 'sigma-squared', conf, '-confidence intervals:'
print '---------------'
print 'symmetric : (', sym_start, ', ', sym_end, ')'
print 'left open : (', left_open_start, ', ', left_open_end, ')'
print 'right open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'