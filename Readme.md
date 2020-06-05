#SAGE scripts for Statistics Exercises

These are my personal simple SAGE scripts that calculate step-by-step solutions for Übungen in Statistik that I have tought in Sommer Semeter 2017-2020. These scripts perform confidence interval calculations and classical statistical hypothesis tests and output the values of intermediate expressions and variables (like sample variance, degrees of freedom, correction coefficients) to the terminal. This helps a lot when checking the homework solutions of students, to find where the first error in the solution occurs.

## What statistical computations it does?

### Parameter tests and confidence intervals

Z-test and confidence interval for the sample mean when variance is known
t-test and confidence interval for the sample mean when variance is unknown
chi^2-test and confidence interval for the sample variance
Z-test and confidence intervals for the difference of two sample means when variances are known
t-test and confidence intervals for the difference of two sample means when variances are unknown but equal
t-test and confidence intervals with Welch correction factor for the difference of two sample means when variances are unknown 
F-test and confidence intervals for the quotient of two sample variances

### Goodness-of-fit to normal distribution tests

chi^2-test
Kolmogorov-Smirnov-test

### Contingency tests

chi^2-independence test for 2x2 table

### Correlation and regression

r-coefficient and the t-test for the correlation
Simple least-squares linear regression y=beta*x + alpha, confidence intervals and hypothesis testing for: alpha, beta, variance of the residues, mean value and variance of E(Y|x) and Y(x) 

### Analysis of Variance:

1-factor ANOVA
2-factor ANOVA

## How this works

Edit the .sage script file, copy-paste the numerical sample data in the list, as

```
x= [1.01, -2.02, 3.03, -5.05]
```

and/or conjectured values of distribution parameters (population mean, variance etc.). Then launch the script in the Sage terminal and copy-paste text output to your file.

## Coded with
SAGE v.8.1. Internaly, it uses [sage.probability](https://doc.sagemath.org/html/en/reference/probability/sage/probability/probability_distribution.html), which is built on [GSL statistics library](https://www.gnu.org/software/gsl/doc/html/statistics.html). For KS-test, it imports [scipy.stats](https://docs.scipy.org/doc/scipy/reference/stats.html)

## Known bugs
Currently, it won't work with Sage 9.1 and higher Python-3 based versions.

## To do:

Separate common functionality into stand-alone scripts
Standardize variable/parameter names accross scripts
Automatize handling of unknown parameters in hypothesis testing
Standardize text output by using Python Tabulate (or similar) package
Make it work with v.9.1 and higher SAGE versions that uses Python 3 (i.e., rewrite all 'print' statements as 'print()')
Think about output to .tex files with rejection region/critical value grapics and automatic generation of solution PDFs.
