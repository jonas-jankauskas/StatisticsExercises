##  @docstring	StatisticsExercises 
# 	@file		ht_diffmeans_vars_unkn_eq_t.sage
#	@brief		t-test for the difference of two sample means with unknown equal variances
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes t-distribution confidence intervals and compares means of two samples with unknown equal variances. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 48
#-------------------------------------------------------

#two samples
x = [80, 25, 75, 75, 90, 95, 105, 45, 40, 70, 95, 40]
y = [50, 40, 60, 55, 70, 45, 25]

#significance and confidence levels
conf = 0.99;
sig = 1-conf;

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#element counts
nx = len(x)
ny = len(y)

#degrees of freedom
dof = nx+ny-2

#sample means
mx = mean(x)
my = mean(y)
mxy = mx - my

#test equal means (difference=0)
muxy = 0;

#population variances are unknown, but are assumed equal: sigma_x = sigma_y

#sample variances:
varx = variance(x)
vary = variance(y)

#pooled sample variance
varpool = ((nx-1)*varx + (ny-1)*vary)/dof

#correction factor
corr = 1/nx + 1/ny

#corrected sample variance and standard dev. of x-y

varxy = varpool*corr
stdxy = sqrt(varxy)

#Student t-distribution quantiles
T = RealDistribution('t', dof)

#t-quantiles
t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf) #conf=1-sig

#difference (x-y) t-score and p-value 
t = (mxy - muxy)/stdxy
p = T.cum_distribution_function(t)

#Confidence intervals
sym_start = mxy - stdxy * t_sym
sym_end = mxy + stdxy * t_sym

left_open_start = -infinity;
left_open_end = mxy + stdxy * t_val

right_open_start = mxy - stdxy * t_val
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '---------------'
print 'first sample size  n = ', nx
print 'second sample size m = ', ny
print 'first sample  x=', [el.n(digits=4) for el in x]
print 'second sample y=', [el.n(digits=4) for el in y]

print '---------------'
print 'x mean   = ', mx.n(digits=4)
print 'y mean   = ', my.n(digits=4)
print 'x var    = ', varx.n(digits=4)
print 'y var    = ', vary.n(digits=4)

print '---------------'
print '(x-y) mean         = ', mxy.n(digits=4)
print 'degrees of freedom = ', dof
print 'pooled var.        = ', varpool.n(digits=4)
print 'correcion          =', corr, ' = ', corr.n(digits=4)
print 'corr. pooled var   = ', varxy.n(digits=4)
print 'corr. pooled dev   = ', stdxy.n(digits=4)

print '---------------'
print 'difference (x-y) t-score and p-values:'
print 't-score = ', t.n(digits=4)
print 'p-value = ', p.n(digits=4);

print '---------------'
print 'degrees of freedom = ', dof
print (100*sig).n(digits=3), '%-significance t-quantiles:'
print 't(sig/2) = ', -t_sym.n(digits=4), ' t(1-sig/2) = ', t_sym.n(digits=4)
print 't(sig)   = ', -t_val.n(digits=4), ' t(1-sig)   = ', t_val.n(digits=4)

print '---------------'
print (100*conf).n(digits=4), '%-confidence intervals for the mean of x-y:'
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