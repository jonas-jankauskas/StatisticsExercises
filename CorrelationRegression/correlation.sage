##  @docstring	StatisticsExercises 
# 	@file		correlation.sage
#	@brief		Two sample correlation
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes correlation coefficient r = cov(X, Y)/sqrt(cov(X, X)*cov(Y, Y)), its confidence intervals and tests hypothesis rho=0. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 52
#-------------------------------------------------------

#two samples
x = [0.14, 0.17, 0.18, 0.21, 0.23, 0.17, 0.18]
y =[0.16, 0.18, 0.18, 0.12, 0.20, 0.20, 0.21]

#significance and confidence levels
conf = 0.95
sig = 1-conf

#-------------------------------------------------------
# internal settings
#-------------------------------------------------------

#default decimal print precision
prec = 5

#print function for lists of reals
def print_real_list(rlist):
	return str([el.n(digits=prec) for el in rlist])
	
#-------------------------------------------------------
# calculations
#-------------------------------------------------------

def cov(a, b):
	num = min(len(a), len(b))
	return ((sum(a[i]*b[i] for i in range(num)) - sum(a)*sum(b)/num)/(num-1))

def cor(a, b):
	return cov(a, b)/sqrt(variance(a)*variance(b))

#correlation coefficient
r = cor(x, y)

#degrees of freedom
n = min(len(x), len(y))
dof = n-2

#t-score
t = r * sqrt(dof/(1-r^2))

#Student t-distribution
T = RealDistribution('t', dof)

#p-value and t-quantiles
pval = T.cum_distribution_function(t)
t_sym = T.cum_distribution_function_inv(1-sig/2)
t_val = T.cum_distribution_function_inv(conf)

#-------------------------------------------------------
#output
#-------------------------------------------------------

print '--------------------------------------------------------------------'
print 'sample sizes    = ', n
print 'first sample  x = ', print_real_list(x)
print 'second sample y = ', print_real_list(y)
print

print '--------------------------------------------------------------------'
print '     ', 'x',                         '          ', 'y'
print 'mean ', mean(x).n(digits=prec),      '      ',     mean(y).n(digits=prec)
print 'var  ', variance(x).n(digits=prec),  '      ',     variance(y).n(digits=prec)
print 'std  ', std(x).n(digits=prec),       '      ',     std(y).n(digits=prec)
print
print 'cov  ', (cov(x, y)).n(digits=prec)
print 'cor  ', (cor(x,y)).n(digits=prec)
print

print '--------------------------------------------------------------------'
print 't-score = ', t.n(digits=prec), ', p-value = ', pval.n(digits=prec)

print '---------------'
print 'degrees of freedom = ', dof
print (100*sig).n(digits=prec-1), '%-significance t-quantiles:'
print 't(sig/2) = ', -t_sym.n(digits=prec), ' t(1-sig/2) = ', t_sym.n(digits=prec)
print 't(sig)   = ', -t_val.n(digits=prec), ' t(1-sig)   = ', t_val.n(digits=prec)

print '--------------------------------------------------------------------'
print 'H_0: qho = 0 vs. H_A: qho <> 0 with ', (100*conf).n(digits=prec),'%-confidence'
if (pval <= 1-sig/2) and (pval >= sig/2):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print
	 
print 'H_0: qho = 0 vs. H_A: qho > 0 with ', (100*conf).n(digits=prec),'%-confidence'
if pval <= 1-sig:
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print

print 'H_0: qho = 0 vs. H_A: qho < 0 with ', (100*conf).n(digits=prec),'%-confidence'
if pval >= sig:
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print