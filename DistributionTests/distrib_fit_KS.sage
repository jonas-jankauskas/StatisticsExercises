##  @docstring	StatisticsExercises 
# 	@file		distrib_fit_KS.sage
#	@brief		KS distribution test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Performs Kolmogorov-Smirnov test against normal distribution. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr
import sage.plot.histogram as hg

#sage.probability does not support KS distribution
#we import it from scipy.stats
import numpy as np
import scipy.stats as st

#-------------------------------------------------------
# data for Exercise 51
#-------------------------------------------------------

#sample
x = [2.45, 2.44, 2.45, 2.46, 2.47, 2.46, 2.45, 2.47, 2.43, 2.44, 2.44, 2.46, 2.43, 2.45, 2.46, 2.45, 2.47, 2.44, 2.46, 2.45]

#test Z-distribution parameters
mu = 2.45
sigma = 0.012

#confidence and significance levels
conf = 0.99
sig = 1-conf


#-------------------------------------------------------
# internal settings
#-------------------------------------------------------

#default decimal print precision
prec = 4

#print function for lists of reals
def print_real_list(rlist):
	return str([el.n(digits=prec) for el in rlist])

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#sample size
nx = len(x);

#sample mean and st.dev.
#mx = mean(x)
#s = std(x)

#standardization
z_scores = [n((el - mu)/sigma) for el in x]
z_scores.sort()

#Z-distribution
Z = RealDistribution('gaussian', 1)

#distinct z-scores
zvals = list(set(z_scores))
zvals.sort()

#Gaussian probabilities of z-scores
p_scores = [Z.cum_distribution_function(z) for z in zvals]

#observed frequencies and their cumulative sums
count = [z_scores.count(z) for z in zvals]
freqs  = [n(cnt)/nx for cnt in count]
cumsum  = [n(sum(count[:i+1]))/nx for i in range(len(count))]

#distance between empirical and Gaussian CDF.
dists_left = [abs(cumsum[i]-p_scores[i+1]) for i in range(len(zvals)-1)]
dists_right = [abs(cumsum[i]-p_scores[i]) for i in range(len(zvals))]
dist = max(max(dists_left), max(dists_right))

#KS-estimate
ksval = n(sqrt(nx)*dist)

#KS-distribution
KS = st.kstwobign()

#KS p-value and and conf-% quantile
pval = KS.cdf(ksval)
cval = KS.ppf(conf)

#-------------------------------------------------------
#output
#-------------------------------------------------------

#pdf plot
gauss_pdf = Z.plot(xmin=zvals[0]-0.1, xmax=zvals[-1]+0.1, thickness=3, color='red');
hist_pdf = histogram(z_scores, bins=len(zvals), normed=True, color='grey')
plot(hist_pdf+gauss_pdf).show()

#cdf plot
gauss_cdf = plot(Z.cum_distribution_function, xmin=zvals[0]-0.1, xmax=zvals[-1]+0.1, thickness=3, color='blue');
hist_cdf = histogram(z_scores, bins=zvals, normed=True, cumulative=True, color='grey')
plot(hist_cdf+gauss_cdf).show()

print '---------------'
print 'size     = ', nx
print 'sample   = ', print_real_list(x)

print '---------------'
print 'conj. pop. mean     mu    = ', mu.n(digits=prec)
print 'conj. pop. st. dev. sigma = ', sigma.n(digits=prec)


print 'z-values   = ', print_real_list(zvals)
print 'p-scores   = ', print_real_list(p_scores)
print 'obs. freqs = ', print_real_list(freqs)
print 'cum. sum   = ', print_real_list(cumsum)

print '---------------'
print 'left dist  = ', print_real_list(dists_left)
print 'right dist = ', print_real_list(dists_right)
print 'max dist   = ', dist.n(digits=prec)
print 'cum. sum   = ', print_real_list(cumsum)
print 'max. dist  = ', dist.n(digits=prec)
print 'KS-score   = ', ksval.n(digits=prec)
print 'p-value    = ', n(pval, digits=prec)

print '---------------'
print (100*sig).n(digits=prec-1), '%-significance KS-quantile:'
print 'ks(1-sig) = ', n(cval, digits=prec)

print '---------------'
print 'H_0: sample is Z(mu, sigma)-distributed vs.\nH_A: sample is not Z(mu, sigma)-distributed:'
if (pval <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'