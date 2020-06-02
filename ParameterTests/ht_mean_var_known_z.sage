##  @docstring	StatisticsExercises 
# 	@file		ht_means_var_known_z.sage
#	@brief		Z-test for the sample mean with known variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes Z-distribution confidence intervals and tests the mean of a sample with known variance. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 41
#-------------------------------------------------------

#one sample
x = [3.57 ,3.66 ,3.48 ,3.62 ,3.88 ,3.80 ,3.82 ,3.73 ,3.23 ,3.58 ,3.70 ,3.52 ,3.54 ,3.34 ,3.62 ,3.28 ,3.33 ,3.46 ,3.38 ,3.68 ,3.68 ,3.81]

#test population mean value
mu = 3.5
#known population variance
sigma = 0.2

#significance and confidence levels
conf = 0.99
sig = 1 - conf

#-------------------------------------------------------
# data for Exercise 44
#-------------------------------------------------------

#test population mean value
#mu = 150
#known population variance
#sigma = 12

#sample size and predicted mean
#nx = 100
#mx = mu

#confidence and significance levels
#conf = 0.95
#sig = 1 - conf

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#sample size, sample mean
nx = len(x)

mx =  mean(x)


#sample standard deviation
sdev = sigma/sqrt(nx)

#sigma2-distribution quantiles
Z = RealDistribution('gaussian', 1)

#z-score and p-value
z = (mx-mu)/(sdev)
p = Z.cum_distribution_function(z)

z_sym = Z.cum_distribution_function_inv(1-sig/2)
z_val = Z.cum_distribution_function_inv(conf)

#confidence intervals
sym_start = mx - z_sym*sdev
sym_end = mx + z_sym*sdev

left_open_start = -infinity
left_open_end = mx + z_val*sdev

right_open_start = mx - z_val*sdev
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '---------------'
print 'sample size = ', nx
print 'sample    x = ', [el.n(digits=4) for el in x]

print '---------------'
print 'conjectured population mean mu = ', mu.n(digits=4);
print 'assumed population st.dev.     = ', sigma.n(digits=4)
print 'sample mean    = ', mx.n(digits=4)
print 'sample st.dev. = ', sdev.n(digits=4)

print '---------------'
print 'z-score = ', z.n(digits=4)
print 'p-value = ', p.n(digits=4)

print '---------------'
print (100*sig).n(digits=4), '%-significance Z-quantiles:'
print 'z(sig/2) = ', -z_sym.n(digits=4), ' z(1-sig/2) = ', z_sym.n(digits=4)
print 'z(sig)   = ', -z_val.n(digits=4), ' z(1-sign)  = ', z_val.n(digits=4)

print '---------------'
print (100*conf).n(digits=4), '%-confidence intervals for mean of x:', mu.n(digits=4), ':'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu <>', mu.n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (mx >= sym_start) and (mx <= sym_end):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu < ', mu.n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (mx >= right_open_start):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: mu = ', mu.n(digits=4), 'against H_A: mu > ', mu.n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (mx <= left_open_end):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'