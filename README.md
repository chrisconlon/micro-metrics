## Empirical Methods for Applied Micro

This is a PhD level course in applied econometrics and computational economics, targeted at students conducting applied research (as opposed to econometricians).

In addition to traditional econometric approaches, this course draws connections to recent literature on machine learning.

BC students should also checkout the [syllabus](syllabus.md). 

Much of the material was (graciously) forked from Chris Conlon's [micro-metrics](https://github.com/chrisconlon/micro-metrics) repo. 

It is based on some combination of 
- [Microeconometrics by Cameron and Trivedi](https://www.amazon.com/Microeconometrics-Methods-Applications-Colin-Cameron/dp/0521848059)
- [Elements of Statistical Learning by Hastie, Tibshirani, and Friedman](https://statweb.stanford.edu/~tibs/ElemStatLearn/)
- Other lectures borrowed/stolen from various sources

The goal is to provide an overview of a number of topics in Microeconometrics including:

1. Nonparametrics and Identification
	- k-NN, Kernels, Nadaraya-Watson
	- Bootstrap and Cross Validation
2. Model Selection and Penalized Regression
	- Ridge, Lasso, LAR, BIC, AIC
3. Treatment Effects and Selection
	- Potential Outcomes, LATE, Diff in Diff, RDD, MTE
4. Binary Discrete Choice (including endogeneity)
	- MLE, Special Regressors, Control Functions
5. Multinomial Discrete Choice
	- Logit, Nested Logit, Mixed Logit
6. Dynamic Discrete Choice
	- Rust Models (NFXP), Hotz+Miller (CCP)
7. Duration Models
8. Bayesian Methods
	- Empirical Bayes, MCMC, James-Stein
