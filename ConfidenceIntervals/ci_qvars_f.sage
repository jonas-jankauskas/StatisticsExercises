##  @docstring	StatisticsExercises 
# 	@file		ci_qvars_f.sage
#	@brief		Two sample variance F-test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes confidence intervals for quatient of two sample variances. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#data from ex. 40
x = [70, 90, 100, 90, 110, 80]
y = [100, 130, 120, 120, 120, 130, 120]

#significance and confidence levels
conf = 0.9;
sig = 1-conf;

#element counts
nx = len(x)
ny = len(y)

#degrees of freedom
dofx = nx-1
dofy = ny-1

#variances
varx = variance(x)
vary = variance(y)

#F-estimate for quotient of variances
qvars = varx/vary

#F-distribution quantiles
F = RealDistribution('F', [dofx, dofy])

f_sym_left = F.cum_distribution_function_inv(sig/2)
f_sym_right = F.cum_distribution_function_inv(1-sig/2)

f_val_left = F.cum_distribution_function_inv(sig)
f_val_right = F.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = qvars/f_sym_right
sym_end = qvars/f_sym_left

left_open_start = 0.0;
left_open_end = qvars/f_val_left

right_open_start = qvars/f_val_right
right_open_end = +infinity

#output
print 'n=', nx
print 'm=', ny

print 'x=', x
print 'y=', y
print

print 'x var    = ', varx
print 'y var    = ', vary
print 'corr. quotient vars = ', qvars
print

print 'degrees of freedom = ', [dofx, dofy] 
print sig, '-significance f-quantiles:'
print 'f(sig/2) = ', f_sym_left, ' f(1-sig/2) = ', f_sym_right
print 'f(sig)   = ', f_val_left, ' f(conf)    = ', f_val_right
print
 
print 'var(x)/var(y) ', conf, '-confidence intervals:'
print '---------------'
print 'symmetric : (', sym_start, ', ', sym_end, ')'
print 'left open : (', left_open_start, ', ', left_open_end, ')'
print 'right open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'