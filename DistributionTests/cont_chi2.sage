##  @docstring	StatisticsExercises 
# 	@file		cont_chi2.sage
#	@brief		chi-2 contigency test
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Tests independence of row and column factors. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr
import sage.plot.histogram as hg

#-------------------------------------------------------
# data for Exercise 50
#-------------------------------------------------------

#sample table
observations = [[14, 33], [32, 21]]

#confidence and significance levels
conf = 0.9
sig = 1-conf

#-------------------------------------------------------
# internal settings
#-------------------------------------------------------

#default decimal print precision
prec = 4

#print function for lists of reals
def print_real_list(rlist):
	print [el.n(digits=prec) for el in rlist]

#recomended minimal average bin count
mincount = 4

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

obs = matrix(observations)

row_sums = [sum(it) for it in obs.rows()]
col_sums = [sum(it) for it in obs.columns()]

numrows = obs.nrows() 
numcols = obs.ncols()
total = sum(row_sums)

a = matrix(numrows, 1, row_sums)
b = matrix(1, numcols, col_sums)
c = matrix(1, 1, total)

obs_table = block_matrix(2, 2, [obs, a, b, c])

obs_probs = obs_table/total
obs_probs.subdivide(obs_table.subdivisions())

exp_table = matrix(numrows, numcols, [[rs*cs/total for rs in row_sums] for cs in col_sums])

exp_vals  = block_matrix(2, 2, [exp_table, a, b, c])
exp_vals.subdivide(obs_table.subdivisions())

exp_probs = exp_vals/total
exp_probs.subdivide(obs_table.subdivisions())

#chi2-score
chi2 = sum([(obs_table[i][j] - exp_table[i][j])^2/exp_table[i][j] for i in range(numrows) for j in range(numcols)])

#Chi2 distribution
dof = (numrows-1)*(numcols-1)
Chi2 = RealDistribution('chisquared', dof)

#Chi2 conf-% quantile and p-score at chi2 value
pval = Chi2.cum_distribution_function(chi2);
cval = Chi2.cum_distribution_function_inv(conf);


#-------------------------------------------------------
#output
#-------------------------------------------------------

print 'observations (& totals)'
print obs_table
print

print '---------------'
print 'probabilities (observed):'
print
print obs_probs.n(digits=prec-1)
print

print '---------------'
print 'probabilities (expected)'
print
print exp_probs.n(digits=prec-1)
print

print '---------------'
print 'expectations:'
print
print exp_vals.n(digits=prec-1)
print

print '---------------'
print 'degrees of freedom = ', dof
print 'chi2-score = ', chi2.n(digits=prec)
print 'p-value    = ', pval.n(digits=prec-1)
print

print '---------------'
print (100*sig).n(digits=prec-1), '%-significance chi2-quantile:'
print 'chi2(1-sig) = ', cval.n(digits=prec)
print

print '---------------'
print 'H_0: row and column events are independant vs. \nH_A: row and column events are dependant:'
if (pval <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'