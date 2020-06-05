##  @docstring	StatisticsExercises 
# 	@file		ht_diffmeans_vars_known_z.sage
#	@brief		Z-test for the difference of two sample means with known variances
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes Z-distribution confidence intervals and compares means of two samples with known variances. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 47
#-------------------------------------------------------

#two samples
x = [80, 25, 75, 75, 90, 95, 105, 45, 40, 70, 95, 40]
y = [50, 40, 60, 55, 70, 45, 25]

#significance and confidence levels
conf = 0.95
sig = 1-conf

#test equal means (difference=0)
muxy = 0;

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#variances are known
varx = 676
vary = 256
 

#sample sizes
nx = len(x)
ny = len(y)

#sample means and their differences
mx = mean(x)
my = mean(y)

mxy = mx-my

#variance and standard dev. of X-Y
varxy = (varx/nx + vary/ny).n()
stdxy = sqrt(varxy)

#Gaussian Distribution
Z = RealDistribution('gaussian', 1)

#difference (x-y) z-score and p-value 
z = (mxy - muxy)/stdxy
p = Z.cum_distribution_function(z)

#Z-quantiles
z_sym = Z.cum_distribution_function_inv(1-sig/2)
z_val = Z.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = mxy - stdxy * z_sym
sym_end   = mxy + stdxy * z_sym

left_open_start = -infinity;
left_open_end = (mxy + stdxy * z_val).n()

right_open_start = (mxy - stdxy * z_val)
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------
print '---------------'
print 'first sample size  n = ', nx
print 'second sample size m = ', ny
print 'first sample       x = ', [el.n(digits=4) for el in x]
print 'Second sample      y = ', [el.n(digits=4) for el in y]

print '---------------'
print 'first sample  x mean         = ', mx.n(digits=4)
print 'second sample y mean         = ', my.n(digits=4)
print 'first  population variance x var = ', varx.n(digits=4)
print 'second population variance y var = ', vary.n(digits=4)

print '---------------'
print '(x-y) difference sample mean = ', mxy.n(digits=4)
print 'difference x-y sample var    = ', varxy.n(digits=4)
print '(x-y) stdev                  = ', stdxy.n(digits=4)

print '---------------'
print 'difference (x-y) z-score and p-values:'
print 'z-score = ', z.n(digits=4)
print 'p-value =', p.n(digits=4);

print '---------------'
print (100*sig).n(digits=3), '%-significance Z-quantiles:'
print 'z(sig/2) = ', -z_sym.n(digits=4), ' z(1-sig/2) = ', z_sym.n(digits=4)
print 'z(sig)   = ', -z_val.n(digits=4), ' z(1-sign)  = ', z_val.n(digits=4)

print '---------------' 
print (100*conf).n(digits=4), '%-confidence intervals for x-y mean:'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left_open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right_open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: mu_x = mu_y, against H_A: mu_x <> mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (z >= -z_sym) and (z <= z_sym):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: mu_x = mu_y, against H_A: mu_x <  mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (z >= -z_val):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: mu_x = mu_y, against H_A: mu_x >  mu_y, with ', (100*conf).n(digits=4),'%-confidence'
if (z <= z_val):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'