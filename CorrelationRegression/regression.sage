##  @docstring	StatisticsExercises 
# 	@file		regression.sage
#	@brief		Least-squares linear regression y=beta*x+alpha
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Performs least-square line fit, calculates confidence intervals and tests hypothesis for regression parameters: alpha, beta, variance of residues, predicted value Y(x) at x=x_0, mean of the predicted value Y(x). Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.
				
reset()

import sage.probability as pr


#-------------------------------------------------------
# data for Exercises 53 and 54
#-------------------------------------------------------

x = [8, 9, 9, 12, 13, 14, 14, 16, 16, 16]
y = [-55, -55, -53, -50, -50, -35, -37, -35, -27, -24]

x_0 = 6
#x_0 = 20

#significance and confidence levels
conf = 0.9
sig = 1-conf

#conf = 0.9 for ex. 53b) and conf = 0.99 for Ex.53 c) and Ex. 54

#-------------------------------------------------------
# control data from Bsp.7.1, S. 73 Statistikskriptum
#-------------------------------------------------------
#x = [4, 9, 10, 14, 4, 7, 12, 22, 1, 17]
#y = [16, 43, 50, 58, 22, 29, 45, 76, 6, 69]

#x_0 = 15
#x_0 = mean(x)

#significance and confidence levels
#conf = 0.95
#sig = 1-conf


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

#sample sizes
num = min(len(x), len(y))

#degrees of freedom
dof = num-2

#covariation
def cov(x, y):
	return ((sum(x[i]*y[i] for i in range(num)) - sum(x)*sum(y)/num)/(num-1))

#correlation
def cor(x, y):
	return cov(x, y)/sqrt(variance(x)*variance(y))
	
#beta-score
beta = cov(x, y)/cov(x, x)

#residual variance and std. dev.
sig2res = (num-1)*(cov(y, y) - beta*cov(x, y))/(num-2)
stdres = sqrt(sig2res)

#beta variance and std. dev.
sig2beta = sig2res/((num-1)*cov(x, x))
stdbeta = sqrt(sig2beta)

#alpha-score, its variance and std. dev.
alpha = mean(y) - beta*mean(x)
sig2alpha = sig2res*(1/num + mean(x)^2/((num-1)*cov(x, x)))
stdalpha = sqrt(sig2alpha)


#mean value of Y at x, its variance and std. dev.
my_x = beta*x_0 + alpha
sig2my_x = sig2res*(1/num + (x_0-mean(x))^2/((num-1)*cov(x, x)))
stdmy_x = sqrt(sig2my_x)

#Y at x, its variance and stddev
sig2y_x = sig2my_x + sig2res
stdy_x = sqrt(sig2y_x)

#default null-values
alpha_0 = alpha
beta_0 = 0
my_x_0 = my_x
sigma2_0 = sig2res

#t-scores for alpha, beta, and mean of y at x:
t_alpha = (alpha-alpha_0)/stdalpha
t_beta = (beta - beta_0)/stdbeta
t_my_x = (my_x - my_x_0)/stdmy_x

#chi2-score for variance of residuals:
chi2_res = dof*sig2res/sigma2_0

#distributions
T = RealDistribution('t', dof)
Chi2 = RealDistribution('chisquared', dof)

#p-scores for alpha, beta, and mean of y at x, and variance of the residuals:
p_alpha = T.cum_distribution_function(t_alpha)
p_beta = T.cum_distribution_function(t_beta)
p_my_x = T.cum_distribution_function(t_my_x)
p_res = Chi2.cum_distribution_function(chi2_res)

	
def t_intervals(mu, stdev, sig):

	#t-quantiles
	t_sym = T.cum_distribution_function_inv(1-sig/2)
	t_val = T.cum_distribution_function_inv(1-sig)

	#confidence interval endpoints
	sym_start = mu - t_sym*stdev
	sym_end = mu + t_sym*stdev

	left_open_start = -infinity
	left_open_end = mu + t_val*stdev
	

	right_open_start = mu - t_val*stdev
	right_open_end = +infinity
	
	#confidence intervals
	symmetric = (sym_start, sym_end)
	left_open = (left_open_start, left_open_end)
	right_open = (right_open_start, right_open_end)
	
	return (symmetric, left_open, right_open)
	
def chi2_intervals(s2, dof, sig):
	
	#chi2-quantiles
	sym_left = Chi2.cum_distribution_function_inv(sig/2)
	sym_right = Chi2.cum_distribution_function_inv(1-sig/2)

	val_left = Chi2.cum_distribution_function_inv(sig)
	val_right = Chi2.cum_distribution_function_inv(1-sig)

	#confidence interval endpoints for sigma^2
	sum_sq = dof*s2
	sym_start = sum_sq/sym_right
	sym_end = sum_sq/sym_left

	left_open_start = 0.0;
	left_open_end = sum_sq/val_left

	right_open_start = sum_sq/val_right
	right_open_end = +infinity
	
	#confidence intervals
	symmetric = (sym_start, sym_end)
	left_open = (left_open_start, left_open_end)
	right_open = (right_open_start, right_open_end)
	
	return (symmetric, left_open, right_open)
	
#confidence intervals for alpha, beta, mu Y at x, and sigma2
ivals_alpha = t_intervals(alpha, stdalpha, sig)
ivals_beta = t_intervals(beta, stdbeta, sig)
ivals_my_x = t_intervals(my_x, stdmy_x, sig)
ivals_y_x = t_intervals(my_x, stdy_x, sig)
ivals_sig2res = chi2_intervals(sig2res, dof, sig) 
	
def print_quantiles(distr, sig):
	
	name = str(distr)
	
	q_sym_left  = distr.cum_distribution_function_inv(sig/2)
	q_sym_right = distr.cum_distribution_function_inv(1-sig/2)
	q_open_left  = distr.cum_distribution_function_inv(sig)
	q_open_right = distr.cum_distribution_function_inv(1-sig)
	
	print (100*sig).n(digits=prec-1), '%-significance ' + name + '-quantiles:'
	print name + '(sig/2) = ', q_sym_left.n(digits=prec), ', ' + name + '(1-sig/2) = ', q_sym_right.n(digits=prec)
	print name + '(sig)   = ', q_open_left.n(digits=prec),', ' + name + '(1-sig)   = ',q_open_right.n(digits=prec)
	
	return

	
	
def print_intervals(ivals, conf, name_msg):
	
	print (100*conf).n(digits=prec), '%-confidence intervals for ' + name_msg + ':'
	print 'symmetric : (', ivals[0][0].n(digits=prec), ', ', ivals[0][1].n(digits=prec), ')'
	print 'left open : (', ivals[1][0].n(digits=prec), ', ', ivals[1][1].n(digits=prec), ')'
	print 'right open: (', ivals[2][0].n(digits=prec), ', ', ivals[2][1].n(digits=prec), ')'
	print
	
	return

	
def decide_hypothesis(param_name, param_value, pval, sig):
	
	#two sided alternative
	msg = 'H_0: ' + param_name + ' = ' + str(param_value.n(digits=prec)) + ' vs.H_A: ' + param_name + ' <> ' + str(param_value.n(digits=prec))+ ' | '
	if (pval >= sig/2) or (pval <= 1-sig/2):
		msg += 'H_0 accepted, H_A rejected!'
	else:
		msg += 'H_0 rejected, H_A accepted!'
	print msg
	
	#left-open alternative
	msg = 'H_0: ' + param_name + ' = ' + str(param_value.n(digits=prec)) + ' vs.H_A: ' + param_name + '  < ' + str(param_value.n(digits=prec))+ ' | '
	if (pval >= sig):
		msg += 'H_0 accepted, H_A rejected!'
	else:
		msg += 'H_0 rejected, H_A accepted!'
	print msg
	
	#right-open alternative
	msg = 'H_0: ' + param_name + ' = ' + str(param_value.n(digits=prec)) + ' vs.H_A: ' + param_name + '  > ' + str(param_value.n(digits=prec))+ ' | '
	if (pval <= 1-sig):
		msg += 'H_0 accepted, H_A rejected!'
	else:
		msg += 'H_0 rejected, H_A accepted!'
	print msg
	
	print
	
	return


#-------------------------------------------------------
#output
#-------------------------------------------------------

print '--------------------------------------------------------------------'
print 'sample sizes    =', num
print 'first sample  x =', print_real_list(x)
print 'second sample y =', print_real_list(y)
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
print 'linear regression y = beta*x + alpha =', beta.n(digits=prec), '*x +', alpha.n(digits=prec)
print '--------------------------------------------------------------------'
print 'mean of Y(' +str(x_0.n(digits=prec)) + ') =', my_x.n(digits=4)
print '--------------------------------------------------------------------' 
print '              ', 'var     ',               '      ', 'std'
print 'residue       ', sig2res.n(digits=prec),   '      ', stdres.n(digits=prec)
print 'alpha         ', sig2alpha.n(digits=prec), '      ', stdalpha.n(digits=prec)
print 'beta          ', sig2beta.n(digits=prec),  '    ', stdbeta.n(digits=prec)
print 'mean Y(' +str(x_0.n(digits=prec)) +')',   sig2my_x.n(digits=prec), '      ', stdmy_x.n(digits=prec)
print 'y(' +str(x_0.n(digits=prec)) +')     ',   sig2y_x.n(digits=prec), '      ', stdy_x.n(digits=prec)

print '--------------------------------------------------------------------'
print 'degrees of freedom = ', dof
print
print '                    ', 't-score',              '     ', 'p-value'
print 'alpha               ', t_alpha.n(digits=prec), '     ', p_alpha.n(digits=prec)
print 'beta                ', t_beta.n(digits=prec),  '     ', p_beta.n(digits=prec)
print 'mean Y(' +str(x_0.n(digits=prec)) +')      ', t_my_x.n(digits=prec),  '     ', p_alpha.n(digits=prec)
print
print '                    ', 'chi2-score',           '  ', 'p-value'
print 'sigma2              ', chi2_res.n(digits=prec),'      ', p_res.n(digits=prec)
print
print_quantiles(T, sig)
print
print_quantiles(Chi2, sig)
print
print '--------------------------------------------------------------------'
print 'confidence intervals:\n'

print_intervals(ivals_alpha, conf, 'alpha')
print_intervals(ivals_beta, conf, 'beta')
print_intervals(ivals_my_x, conf, 'mean of Y at x='+str(x_0))
print_intervals(ivals_y_x, conf, 'y('+str(x_0)+')')
print_intervals(ivals_sig2res, conf, 'sigma^2')

print '--------------------------------------------------------------------'
decide_hypothesis('alpha', alpha_0, p_alpha, sig)
decide_hypothesis('beta', beta_0, p_beta, sig)
decide_hypothesis('Y('+str(x_0.n(digits=prec))+')', my_x_0, p_my_x, sig)
decide_hypothesis('sigma2', sigma2_0, p_res, sig)