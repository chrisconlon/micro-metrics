## Empirical Methods for Applied Micro

This is a PhD level course in applied econometrics and computational economics, targeted at students conducting applied research (as opposed to econometricians).

In addition to traditional econometric approaches, this course draws connections to recent literature on machine learning.

BC students should also checkout the [syllabus](syllabus.md). 

Much of the material was (gratefully) forked from Chris Conlon's [micro-metrics](https://github.com/chrisconlon/micro-metrics) repo. 

It is based on some combination of 
- [Microeconometrics by Cameron and Trivedi](https://www.amazon.com/Microeconometrics-Methods-Applications-Colin-Cameron/dp/0521848059)
- [Elements of Statistical Learning by Hastie, Tibshirani, and Friedman](https://statweb.stanford.edu/~tibs/ElemStatLearn/)
- Other lectures borrowed/stolen from various sources
	- Bruce Hansen's [online text](https://www.ssc.wisc.edu/~bhansen/econometrics/)
	
The goal is to provide an overview of a number of topics in Microeconometrics including:

1. Nonparametrics and Identification
	- Density estimatation, k-NN, Kernels, Nadaraya-Watson
	- Bootstrap and Cross Validation
2. Model Selection and Penalized Regression
	- Ridge, Lasso, LAR, BIC, AIC
3. Treatment Effects and Selection
	- Potential Outcomes, LATE, Diff in Diff, RDD, MTE
4. Binary Discrete Choice (including endogeneity)
	- MLE, Special Regressors, Control Functions
5. Computational 
	- Root finding, Optimization
	- Differentiation, Integration
6. Multinomial Discrete Choice
	- Logit, Nested Logit, Mixed Logit
7. Dynamic Discrete Choice
	- Rust Models (NFXP), Hotz+Miller (CCP)
8. Partial Identificaton

Repo material likely not covered: 
- Duration Models
- Bayesian Methods
	- Empirical Bayes, MCMC, James-Stein
