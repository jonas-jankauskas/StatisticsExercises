##  @docstring	StatisticsExercises 
# 	@file		2-way anova.sage
#	@brief		Two-factor analysis of variance
#	@author		Jonas Jankauskas
#	@date		June 6, 2022
#	@version	1.0
#	@note		SAGE version 8.1
#	@details	Computes row and column mean sums of squares, F-tests the equality of means accross rows and columns. Prints out detailed data report and intermediate variables.
#   @todo		Automatize the handling of parameter null-values. Standardize the variable names and result output accross scripts. Separate common functionality.

reset()

import sage.probability as pr

#-------------------------------------------------------
# data for Exercise 56
#-------------------------------------------------------

#sample table
table = [
[72, 81, 63, 59],
[83, 94, 91, 86],
[68, 57, 73, 61],
[55, 73, 77, 66],
]

#confidence and significance levels
conf = 0.99
sig = 1-conf


#-------------------------------------------------------
# control data from Bsp. 6.2., S.66 Statistik Skriptum
#-------------------------------------------------------
#table = [
#[56.7, 45.7, 48.3, 54.6, 37.7],
#[64.5, 53.4, 54.3, 57.5, 52.3],
#[56.7, 50.6, 49.5, 56.5, 44.7],
#]

# control confidence and significance levels
#conf = 0.95
#sig = 1-conf

# control answers:
#F_{0.95}(2,8) = 4.46, F_{0.95}(4,8) = 3.84
#f_row = 13.1 > 4.46: H_0: row means are equal rejected
#f_col = 16.2 > 3.84: H_0: col means are equal rejected

#-------------------------------------------------------
# internal settings
#-------------------------------------------------------

#default decimal printprecision
prec = 5

#print function for lists of reals
def print_real_list(rlist):
	return str([el.n(digits=prec) for el in rlist])

#-------------------------------------------------------
# calculations
#-------------------------------------------------------

#no. of rows, no. of entries in each row, total no. of entries
num_rows = len(table)
num_cols = len(table[0])
total = num_rows*num_cols

#trnasposed table (columns)
transposed_table = [[table[i][j] for i in range(num_rows)] for j in range(num_cols)]

#row, column and table means
row_mean = [mean(row) for row in table]
col_mean = [mean(col) for col in transposed_table]
table_mean = n(sum(sum(row) for row in table))/total

#between-rows sums-of-squares, degrees-of-freedom and mean-sum-of-squares
SSQ_br = sum(num_cols*r_mean^2 for r_mean in row_mean) - total*table_mean^2
dof_br = num_rows - 1
MSSQ_br = n(SSQ_br)/dof_br

#between-columns sums-of-squares, degrees-of-freedom and mean-sum-of-squares
SSQ_bc = sum(num_rows*c_mean^2 for c_mean in col_mean) - total*table_mean^2
dof_bc = num_cols - 1
MSSQ_bc = n(SSQ_bc)/dof_bc

#table sum of squares, degrees-of-freedom and mean-sum-of-squares
sum_sq = sum(sum(num^2 for num in row) for row in table)
SSQ_tab = sum_sq - total*table_mean^2
SSQ_err = SSQ_tab - SSQ_br - SSQ_bc
dof_err = dof_br*dof_bc
MSSQ_err = n(SSQ_err)/dof_err

#f-estimates
f_rows= MSSQ_br/MSSQ_err
f_cols= MSSQ_bc/MSSQ_err

#F-distribution
F_rows = RealDistribution('F', [dof_br, dof_err])
F_cols = RealDistribution('F', [dof_bc, dof_err])

#p-values and conf-% quantiles
pval_rows = F_rows.cum_distribution_function(f_rows)
cval_rows = F_rows.cum_distribution_function_inv(conf)

pval_cols = F_cols.cum_distribution_function(f_cols)
cval_cols = F_cols.cum_distribution_function_inv(conf)

#-------------------------------------------------------
#output
#-------------------------------------------------------

#data table augmented with row means and sums-of-squares

A = matrix(table)
B= transpose(matrix(row_mean))
C = matrix(col_mean)
D = matrix(1, 1, table_mean)

M = block_matrix(2, 2, [A, B, C, D], subdivide=True)

print
print 'data table with row, column and total means:'
print '--------------------------------------------------------------------'
print M.n(digits=prec)
print ''

print '////////////////////////////////////////////////////////////////////'
print ''

print 'entries : ', total
print 'type            dof   SSQ     MSSQ     F-scores  p-values  f(1-sig)'
print '--------------------------------------------------------------------'
print 'between-rows    ', dof_br, ' ', SSQ_br.n(digits=prec),  ' ', MSSQ_br.n(digits=prec), ' ', f_rows.n(digits=prec), '  ', pval_rows.n(digits=prec), ' ', cval_rows.n(digits=prec)
print 'between-columns ', dof_bc, ' ', SSQ_bc.n(digits=prec),  ' ', MSSQ_bc.n(digits=prec), ' ', f_cols.n(digits=prec), '  ', pval_rows.n(digits=prec), ' ', cval_cols.n(digits=prec)
print 'accross table   ', dof_err,' ', SSQ_err.n(digits=prec), ' ', MSSQ_err.n(digits=prec)
print 'total           ', num_rows*num_cols-1, '', SSQ_tab.n(digits=prec)
print

print '--------------------------------------------------------------------'
print 'H_0: row means are equal vs.H_A: row means are not all equal'
if (pval_rows <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print
	
print 'H_0: column means are equal vs.H_A: column means are not all equal'
if (pval_cols <= conf):
	print 'H_0 accepted, H_A rejected!'
else:
	print 'H_0 rejected, H_A accepted!'
print