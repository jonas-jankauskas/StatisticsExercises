##  @docstring	StatisticsExercises 
# 	@file		ci_diffmeans_vars_unkn_eq_t.sage
#	@brief		t-test for the difference of two sample means with unknown equal variances
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes t-distribution confidence intervals for the difference of two sample means when variances are unknown but equal. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#data from ex.38
x = [49.7, 51.9, 49.0, 49.6, 51.2, 49.2, 49.8, 49.9, 46.5, 51.7]
y = [48.4, 46.3, 51.1, 49.2, 48.7, 49.4, 49.9, 48.3]

#significance and confidence levels
conf = 0.9;
sig = 1-conf;

#element counts
nx = len(x)
ny = len(y)

#degrees of freedom
dof = nx+ny-2

#means
mx = mean(x)
my = mean(y)
mxy = mx-my

#variances unknow, assumed equal
varx = variance(x)
vary = variance(y)

#pooled variance
varpool = ((nx-1)*varx + (ny-1)*vary)/dof

#correction factor
corr = 1/nx + 1/ny

#variance and standard dev. of X-Y
varxy = varpool*corr
stdxy = sqrt(varxy)

#Student t-distribution quantiles
T = RealDistribution('t', dof)

t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf) #conf=1-sig

#Confidence intervals
sym_start = mxy - stdxy * t_sym
sym_end = mxy + stdxy * t_sym

left_open_start = -infinity;
left_open_end = mxy + stdxy * t_val

right_open_start = mxy - stdxy * t_val
right_open_end = +infinity

#output
print 'n=', nx
print 'm=', ny

print 'x=', x
print 'y=', y
print

print 'x mean   = ', mx
print 'y mean   = ', my
print 'x-y mean = ', mxy
print 'x var    = ', varx
print 'y var    = ', vary
print 'corr. pooled var = ', varpool
print 'corr. pooled dev = ', sqrt(varpool)
print

print 'degrees of freedom = ', dof
print sig, '-significance t-quantiles:'
print 't(1-sig/2) = ', t_sym
print 't(sig)     = ', t_val
print
 
print 'x-y mean ', conf, '-confidence intervals:'
print '---------------'
print 'symmetric: (', sym_start, ', ', sym_end, ')'
print 'left open: (', left_open_start, ', ', left_open_end, ')'
print 'right open: (', right_open_start, ', ', right_open_end, ')'
print '---------------'