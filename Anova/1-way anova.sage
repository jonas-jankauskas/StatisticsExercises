##  @docstring	StatisticsExercises 
# 	@file		1-way anova.sage
#	@brief		One-factor analysis of variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes mean sums of squares accross rows, F-tests the equality of means. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 55
#-------------------------------------------------------

#sample table
table = [
[17.5, 16.9, 15.8, 18.6],
[16.4, 19.2, 17.7, 15.4],
[20.3, 15.7, 17.8, 18.9],
]

#confidence and significance levels
conf = 0.95
sig = 1-conf

#-------------------------------------------------------
# data for Exercise 56
#-------------------------------------------------------

#sample table
#table = [
#[72, 81, 63, 59],
#[83, 94, 91, 86],
#[68, 57, 73, 61],
#[55, 73, 77, 66],
#]

#confidence and significance levels
#conf = 0.99
#sig = 1-conf


# control example from Statistikskriptum Bsp 6.1. (S.63)
#table = [
#[48.3, 49.6, 48.4, 48.1],
#[56.1, 56.3, 56.8],
#[52.1, 51.1, 51.6, 52.1, 51.1],
#]

#control confidence and significance levels
#conf = 0.95
#sig = 1-conf

#uncomment to perform ANOVA on columns instead of rows
#transposed sample table
#table = [[row[i] for row in table] for i in range(len(table))]

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

#no. of rows, no. of entries in each row, total no. of entries
num_rows = len(table)
row_len = [len(row) for row in table]
total = sum(num for num in row_len)

#sample row means and total mean
row_mean = [mean(row) for row in table]
table_mean = n(sum(sum(row) for row in table))/total

#between-rows sums-of-squares, degrees-of-freedom and mean-sum-of-squares
SSQ_b = sum(row_len[i]*row_mean[i]^2 for i in range(num_rows)) - total*table_mean^2
dof_b = num_rows - 1
MSSQ_b = n(SSQ_b)/dof_b

#within-rows sum of squares, degrees-of-freedom and mean-sum-of-squares
row_sum_sq = [sum(num^2 for num in row) for row in table]
SSQ_t = sum(row_sum_sq) - total*table_mean^2
SSQ_w = SSQ_t - SSQ_b
dof_w = total - num_rows
MSSQ_w = n(SSQ_w)/dof_w

#f-estimate
f= MSSQ_b/MSSQ_w

#F-distribution
F = RealDistribution('F', [dof_b, dof_w])

#p-value and conf-% level quantile
pval = F.cum_distribution_function(f)
cval = F.cum_distribution_function_inv(conf)

#-------------------------------------------------------
#output
#-------------------------------------------------------

#data table augmented with row means and sums-of-squares

#if row lengths are different, fill-in with the 0
max_row_len = max(row_len)
for i in range(num_rows):
	table[i] += (max_row_len - len(table[i]))*[NaN]

A = matrix(QQ, table)
B= A.augment(vector(row_mean), subdivide=True)
C = B.augment(vector(row_sum_sq), subdivide=True)

print
print 'data table with row means and row sums of squares'
print '-----------------------------------------------------------------'
print C.n(digits=prec)
print

print '/////////////////////////////////////////////////////////////////'
print

print 'entries:', total
print 'table mean  = ', table_mean.n(digits=prec)
print 'type          dof  SSQ     MSSQ    F-scores  p-values  f(1-sig)'
print '-----------------------------------------------------------------'         
print 'between-rows :', dof_b, ' ', SSQ_b.n(digits=prec), '', MSSQ_b.n(digits=prec), '', f.n(digits=prec), ' ', pval.n(digits=prec), ' ', cval.n(digits=prec) 
print 'within-rows  :', dof_w, ' ', SSQ_w.n(digits=prec), '', MSSQ_w.n(digits=prec)
print 'accross table:', total-1,'', SSQ_t.n(digits=prec), '  '
print

print '-----------------------------------------------------------------'
print 'H_0: row means are are equal vs.H_A: row means are not all equal'
if (pval <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print