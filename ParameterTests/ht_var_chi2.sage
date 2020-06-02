##  @docstring	StatisticsExercises 
# 	@file		ht_var_chi2.sage
#	@brief		Sample variance chi-squared test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes confidence intervals and tests hypothesis about variance. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Move common functionality into a separate scripts.

reset()

import sage.probability as pr


#-------------------------------------------------------
# data for Exercise 43
#-------------------------------------------------------

#one sample
x = [3.57 ,3.66 ,3.48 ,3.62 ,3.88 ,3.80 ,3.82 ,3.73 ,3.23 ,3.58 ,3.70 ,3.52 ,3.54 ,3.34 ,3.62 ,3.28 ,3.33 ,3.46 ,3.38 ,3.68 ,3.68 ,3.81]

#conjectured population stdev
sigma = 0.2

#significance and confidence levels
conf = 0.95
sig = 1-conf

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#sample size
nx = len(x)

#sample mean
mx = mean(x)

#degrees of freedom
dof = nx-1

#sample variance and de-averaged sample sum-of-squares
varx  = variance(x)
sumsq = dof*varx

#chi2-distribution with n-1 degrees of freedom
Chi2 = RealDistribution('chisquared', dof)

#chi2-score and p-value 

chi2 = sumsq/sigma^2
p = Chi2.cum_distribution_function(chi2)

#chi2-quantiles
chi2_sym_left = Chi2.cum_distribution_function_inv(sig/2)
chi2_sym_right = Chi2.cum_distribution_function_inv(1-sig/2)

chi2_val_left = Chi2.cum_distribution_function_inv(sig)
chi2_val_right = Chi2.cum_distribution_function_inv(conf)

#confidence intervals for sigma^2
sym_start = sumsq/chi2_sym_right
sym_end = sumsq/chi2_sym_left

left_open_start = 0.0;
left_open_end = sumsq/chi2_val_left

right_open_start = sumsq/chi2_val_right
right_open_end = +infinity

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '---------------'
print 'sample size    = ', nx
print 'first sample x = ', [el.n(digits=4) for el in x]

print '---------------'
print 'sample mean = ', mx.n(digits=4)
print 'sample var. = ', varx.n(digits=4)
print 'conj. sigma = ', sigma.n(digits=4)

print '---------------'
print 'degrees of freedom = ', dof
print 'chi2-score = ', chi2.n(digits=5)
print 'p-value    = ', p.n(digits=4)

print '---------------'
print (100*sig).n(digits=3), '%-significance chi2-quantiles:'
print 'chi2(sig/2) = ', chi2_sym_left.n(digits=4), ' chi2(1-sig/2) = ', chi2_sym_right.n(digits=4)
print 'chi2(sig)   = ', chi2_val_left.n(digits=4), ' chi2(1-sig)    = ', chi2_val_right.n(digits=4)

print '---------------'
print (100*conf).n(digits=4), '%-confidence intervals for sigma^2:'
print 'symmetric : (', sym_start.n(digits=4), ', ', sym_end.n(digits=4), ')'
print 'left open : (', left_open_start.n(digits=4), ', ', left_open_end.n(digits=4), ')'
print 'right open: (', right_open_start.n(digits=4), ', ', right_open_end.n(digits=4), ')'

print '---------------'
print 'H_0: sigma^2 =', (sigma^2).n(digits=4), 'against H_A: sigma^2 <> ', (sigma^2).n(digits=4), 'with ', (100*conf).n(digits=4),'%-confidence'
if (sigma >= sym_start) and (sigma <= sym_end):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
	 
print 'H_0: sigma^2 =', (sigma^2).n(digits=4), 'against H_A: sigma^2 <  ', (sigma^2).n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (sigma >= left_open_end):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'

print 'H_0: sigma^2 =', (sigma^2).n(digits=4), 'against H_A: sigma^2 >  ', (sigma^2).n(digits=4), ' with ', (100*conf).n(digits=4),'%-confidence'
if (sigma <= right_open_start):
	 print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'