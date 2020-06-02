##  @docstring	StatisticsExercises 
# 	@file		distrib_fit_chi2.sage
#	@brief		chi-2 distribution test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Performs chi-2 test against normal distribution. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr
import sage.plot.histogram as hg


#-------------------------------------------------------
# data for Exercise 49
#-------------------------------------------------------

#sample
x = [440, 433, 429, 438, 432, 445, 435, 441, 438, 442, 420, 407, 428, 425, 444, 449, 440, 443, 461, 425]

#conjectured  population mean
sigma = 11.0

#confidence and significance levels
conf = 0.95
sig = 1-conf


#-------------------------------------------------------
# internal settings
#-------------------------------------------------------

#default decimal print precision
prec = 4

#print function for lists of reals
def print_real_list(rlist):
	return str([el.n(digits=prec) for el in rlist])

#recomended minimal average bin count
mincount = 4

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#sample size
nx = len(x);

#sample mean and st.dev.
mx = mean(x)
s = std(x)

#number of estimated params (mean) = 1
num_params = 1

#standardization
z_scores = [n((val - mx)/sigma) for val in x]
z_scores.sort()

#number of bin and bin width
num_bin  = len(z_scores) // mincount
bin_width = (max(z_scores)-min(z_scores))/num_bin

#evenly distributed bin endpoints in the range
bin = [-oo]+[min(z_scores)+ j*bin_width for j in range(1, num_bin)]+[+oo];
bin.sort()

#bin probabilities
Z = RealDistribution('gaussian', 1)
p_scores = [Z.cum_distribution_function(z) for z in bin]
p_values = [p_scores[j+1] - p_scores[j] for j in range(len(p_scores)-1)]

#expected bin counts
e_values = [nx * p for p in p_values];

#observed bin counts
o_values = [len([z for z in z_scores if ((z > bin[j]) and (z <= bin[j+1]))]) for j in range(len(bin)-1)]

#chi2-score
chi2 = sum([(o_values[j] - e_values[j])^2/e_values[j] for j in range(len(bin)-1)])

#Chi2 distribution
dof = num_bin - num_params - 1
Chi2 = RealDistribution('chisquared', dof)

#Chi2 conf-% quantile and p-score at chi2 value
pval = Chi2.cum_distribution_function(chi2);
cval = Chi2.cum_distribution_function_inv(conf);

#-------------------------------------------------------
#output
#-------------------------------------------------------

#pdf plot
gauss_pdf = Z.plot(xmin=z_scores[0]-0.1, xmax=z_scores[-1]+0.1, thickness=3, color='red');
hist_pdf = histogram(z_scores, bins=len(bin)-1, normed=True, color='grey')
plot(hist_pdf+gauss_pdf).show()

#cdf plot
gauss_cdf = plot(Z.cum_distribution_function, xmin=z_scores[0]-0.1, xmax=z_scores[-1]+0.1, thickness=3, color='blue');
hist_cdf = histogram(z_scores, bins=len(bin)-1, normed=True, cumulative=True, color='grey')
plot(hist_cdf+gauss_cdf).show()


print '---------------'
print 'size     = ', nx
print 'sample   = ', print_real_list(x)
print 'z-scores = ', print_real_list(z_scores)

print '---------------'
print 'sample mean           mx = ', mx.n(digits=prec)
#print 'sample st. dev.         = ', s.n(digits=prec)
print 'population st.dev. sigma = ', sigma.n(digits=prec)


print '---------------'
print 'no. of bin  = ', num_bin
print 'bin         = ', print_real_list(bin)
print 'obs. counts = ', o_values
print 'bin probs   = ', print_real_list(p_values)
print 'exp. count  = ', print_real_list(e_values)

print '---------------'
print 'degrees of freedom = ', dof
print 'chi2-score = ', chi2.n(digits=prec)
print 'p-value    = ', pval.n(digits=prec)

print '---------------'
print (100*sig).n(digits=prec-1), '%-significance chi2-quantile:'
print 'chi2(1-sig) = ', cval

print '---------------'
print 'H_0: sample is Z-distributed against H_A: sample is not Z-distributed:'
if (pval <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'