##  @docstring	StatisticsExercises 
# 	@file		ht_qvars_f.sage
#	@brief		Two sample variance F-test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes confidence intervals and performs variance comparison tests. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 46
#-------------------------------------------------------

#two samples
x = [80, 25, 75, 75, 90, 95, 105, 45, 40, 70, 95, 40]
y = [50, 40, 60, 55, 70, 45, 25]

#significance and confidence levels
conf = 0.9;
sig = 1-conf;

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#element counts
nx = len(x)
ny = len(y)

#degrees of freedom
dofx = nx-1
dofy = ny-1

#variances
varx = variance(x)
vary = variance(y)

#F-distribution
F = RealDistribution('F', [dofx, dofy])

#f-score and p-value
f = varx/vary
p = F.cum_distribution_function(f);

#F-distribution quantiles

f_sym_left = F.cum_distribution_function_inv(sig/2)
f_sym_right = F.cum_distribution_function_inv(1-sig/2)

f_val_left = F.cum_distribution_function_inv(sig)
f_val_right = F.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = f/f_sym_right
sym_end = f/f_sym_left

left_open_start = 0.0;
left_open_end = f/f_val_left

right_open_start = f/f_val_right
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print 'first sample  size n = ', nx
print 'second sample size m = ', ny
print 'first sample       x = ', [el.n(digits=4) for el in x]
print 'second sample      y = ', [el.n(digits=4) for el in y]

print '---------------'
print 'first sample variance  x var = ', varx.n(digits=4)
print 'second sample variance y var = ', vary.n(digits=4)
print 'f-score = ', f.n(digits=4)
print 'p-value = ', p.n(digits=4)

print '---------------'
print 'degrees of freedom = ', [dofx, dofy] 
print (100*sig).n(digits=3), '-significance f-quantiles:'
print 'f(sig/2) = ', f_sym_left.n(digits=4), ' f(1-sig/2) = ', f_sym_right.n(digits=4)
print 'f(sig)   = ', f_val_left.n(digits=4), ' f(conf)    = ', f_val_right.n(digits=4)

print '---------------'
print (100*conf).n(digits=4), '%-confidence intervals for var(x)/var(y):'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: sigma_x^2 = sigma_y^2, against H_A: sigma_x^2 <> sigma_y^2, with ', (100*conf).n(digits=4),'%-confidence'
if (f >= f_sym_left) and (f <= f_sym_right):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: sigma_x^2 = sigma_y^2, against H_A: sigma_x^2 <  sigma_y^2, with ', (100*conf).n(digits=4),'%-confidence'
if (f >= f_val_left):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: sigma_x^2 = sigma_y^2, against H_A: sigma_x^2 >  sigma_y^2, with ', (100*conf).n(digits=4),'%-confidence'
if (f <= f_val_right):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'