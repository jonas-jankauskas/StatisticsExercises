##  @docstring	StatisticsExercises 
# 	@file		ht_mean_var_ukn_t.sage
#	@brief		t-test for the sample mean with unknown variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes t-distribution confidence intervals and tests the mean of a sample with unknown variance. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 42
#-------------------------------------------------------

#one sample
x = [3.57 ,3.66 ,3.48 ,3.62 ,3.88 ,3.80 ,3.82 ,3.73 ,3.23 ,3.58 ,3.70 ,3.52 ,3.54 ,3.34 ,3.62 ,3.28 ,3.33 ,3.46 ,3.38 ,3.68 ,3.68 ,3.81]

#test population mean
mu = 3.5

#confidence and significance levels

conf = 0.95
sig = 1 - conf

#-------------------------------------------------------
# data for Exercise 45
#-------------------------------------------------------

#one sample
#x = [1.19, 1.23, 1.18, 1.21, 1.27, 1.14]

#test population mean
#mu = 1.26

#confidence and significance
#conf = 0.9
#sig = 1 - conf

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#sample size, sample mean, sample variance

nx = len(x)
mx = mean(x)
stdv = std(x)

#degrees of freedom
dof = nx-1

#sample variance estimate
sdev = stdv/sqrt(nx)

#t-distribution 
T = RealDistribution('t', dof)

#t-score and p-value
t = (mx-mu)/(sdev)
p = T.cum_distribution_function(t)

#t-quantiles
t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = mx - t_sym*sdev
sym_end = mx + t_sym*sdev

left_open_start = -infinity
left_open_end = mx + t_val*sdev

right_open_start = mx - t_val*sdev
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '---------------'
print 'sample size    = ', nx
print 'first sample x = ', [el.n(digits=4) for el in x]

print '---------------'
print 'conjectured population mean mu = ', mu.n(digits=4);
print 'sample mean                    = ', mx.n(digits=4)
print 'sample stdev                   = ', stdv.n(digits=4)
print 'population stdev estimate      = ', sdev.n(digits=4)

print '---------------'
print 't-score = ', t.n(digits=4)
print 'p-value = ', p.n(digits=4)

print '---------------'
print 'degrees of freedom = ', dof
print (100*sig).n(digits=3), '%-significance t-quantiles:'
print 't(sig/2) = ', -t_sym.n(digits=4), ' t(1-sig/2) = ', t_sym.n(digits=4)
print 't(sig)   = ', -t_val.n(digits=4), ' t(1-sig)   = ', t_val.n(digits=4)

print '---------------'
print (100*conf).n(digits=4), '%-confidence intervals for population mean:'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu <> ', mu.n(digits=4), 'with ', (100*conf).n(digits=4),'%-confidence'
if (mx >= sym_start) and (mx <= sym_end):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu <  ', mu.n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (mx >= right_open_start):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu >  ', mu.n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (mx <= left_open_end):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'