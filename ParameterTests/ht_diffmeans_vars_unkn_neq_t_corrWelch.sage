##  @docstring	StatisticsExercises 
# 	@file		ht_diffmeans_vars_unkn_eq_t.sage
#	@brief		Welsch-corrected t-test for the difference of two sample means with unknown variances
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes Welsch-corrected t-distribution confidence intervals and compares means of two samples with unknown and possibly unequal variances. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 48
#-------------------------------------------------------

#two samples
#x = [49.7, 51.9, 49.0, 49.6, 51.2, 49.2, 49.8, 49.9, 46.5, 51.7]
#y = [48.4, 46.3, 51.1, 49.2, 48.7, 49.4, 49.9, 48.3]

x=[80, 25, 75, 75, 90, 95, 105, 45, 40, 70, 95, 40]
y=[50,40, 60, 55, 70, 45, 25]


#population variances unknow and possibly unequal

#test for equal means (difference=0)
muxy = 0;

#significance and confidence levels
conf = 0.99;
sig = 1-conf;

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#element counts
nx = len(x)
ny = len(y)

#means
mx = mean(x)
my = mean(y)
mxy = mx-my

#sample variances
varx = variance(x)
vary = variance(y)

#Welch-corrected variance and deviation
cvar = varx/nx + vary/ny
cdev = sqrt(cvar)

#Welch-corrected degrees of freedom
wcof = cvar^2/((varx/nx)^2/(nx-1)+(vary/ny)^2/(ny-1))
cdof = floor(wcof)


#Student t-distribution
T = RealDistribution('t', cdof)

#t-quantiles
t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf)

#Welch-corrected t-score and p-value
t = (mxy - muxy)/cdev
p = T.cum_distribution_function(t);

#confidence intervals
sym_start = mxy - cdev * t_sym
sym_end = mxy + cdev * t_sym

left_open_start = - infinity;
left_open_end = mxy + cdev * t_val

right_open_start = mxy - cdev * t_val
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '---------------'
print 'first sample size  n = ', nx
print 'second sample size m = ', ny
print 'first sample  x = ', [el.n(digits=4) for el in x]
print 'second sample y = ', [el.n(digits=4) for el in y]

print '---------------'
print 'first sample  x mean       = ', mx.n(digits=4)
print 'second sample y mean       = ', my.n(digits=4)
print 'difference of sample means = ', mxy.n(digits=4)
print 'first sample  x var        = ', varx.n(digits=4)
print 'second sample y var        = ', vary.n(digits=4)

print '---------------'
print 'Welch correction = ', wcof.n(digits=4)
print 'Welch-corr. var  = ', cvar.n(digits=4)
print 'Welch-corr. dev  = ', cdev.n(digits=4)

print '---------------'
print 'Welch-corrected difference (x-y) t-score and p-value:'
print 't-score = ', t.n(digits=4)
print 'p-value = ', p.n(digits=4);

print '---------------'
print 'Welch-corrected degrees of freedom = ', cdof
print (100*sig).n(digits=3), '%-significance t-quantiles:'
print 't(sig/2) = ', -t_sym.n(digits=4), ' t(1-sig/2) = ', t_sym.n(digits=4)
print 't(sig)   = ', -t_val.n(digits=4), ' t(1-sig)   = ', t_val.n(digits=4)

print '---------------' 
print (100*conf).n(digits=4), '%-confidence intervals for difference x-y mean:'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: mu_x = mu_y, against H_A: mu_x <> mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (t >= -t_sym) and (t <= t_sym):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: mu_x = mu_y, against H_A: mu_x <  mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (t >= -t_val):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: mu_x = mu_y, against H_A: mu_x >  mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (t <= t_val):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'