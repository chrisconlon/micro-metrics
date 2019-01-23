# Problem Set 1 - Auctions
## Prof. Richard L. Sweeney

In the Nonparametrics lecture, we discussed a prominent application of nonparametric density estimation: auctions. In these problems, we are typically interested in recovering the distribution of values for some object. Once we have these, we can project onto covariates of interest, or simulate an alternative mechanism. 

Guerre, Perrigne and Voung (ECMA 2000) show that the distibution of (unobserved) bidder values can be recovered nonparametrically if we know how the auction works and can estimate the pdf and cdf of bids. 

In this problem, we are going to simulate $N$ auctions, and then try to recover the (known) distribution of valuations ($v$) using GPV. 

[Hickman (IJIO, 2010)](https://www.sciencedirect.com/science/article/pii/S0167718709001076) provides code to simulate bid functions from several distributions in Matlab. This code, available on [his website](http://home.uchicago.edu/hickmanbr/research.html), has been placed in this directory. 

1. Simulate bids from 100 IPV first price auctions with 3 bidders. 
    - To get the bid function, call the `FPauction.m` function. 
    - Assume the distribution of values is **lognormal**, and that the mean and standard deviation of the underlying normal are $\mu = 1$ and $\sigma = .5$ 
    - Note that in the `FPauction.m` these parameters are $a$ and $b$. Also, to speed things up assume a max bid of `vbar` $= \exp(2.5)$

2. Estimate the density of bids using an assumed normal distribution, a Gaussian kernel, and an Epanechnikov kernel. For the kernels, use a plug-in estimate for the bandwidth. 
    - **Bonus:** Try using least-squares cross-validation to pick the bandwidth!

3. Plot the estimated density functions, overlaid on top of a histogram of the bids. Which one fits best?

4. Use GPV and the Epanechnikov kernel to recover the valuation implied for each bid. 
$$ \hat v = b + \frac{\hat G_B (b)}{(n-1)\hat g_B(b)},$$ where $n=3$ is the number of bidders.

5. Finally estimate the distribution of $v$ using another Epanechnikov kernel. 

6. Plot this estimated density against the known density that generated the data. 

7. Use this distribution to construct an estimate of the impact of adding an additional bidder on seller revenues. Steps:
    - Approximate the recovered distribution of values $\hat v$ using a log normal. 
    - Use these estimated log normal parameters to recompute the bid function for actions with 3 and 4 bidders. 
    - Simulate $M = 100$ sets of $N=100$ auction pairs (with 3 and 4 bidders), drawing $v$ from the recovered set of values $\hat v$ with replacement.
    - For each $m \in M$, compute $\theta_m$ equals the expected difference in the winning bid going from 3 to 4 bidders. 
    - Report 95% confidence intervals on $\theta$
