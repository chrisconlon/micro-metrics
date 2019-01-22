# Problem set 1 - Auctions
## Prof. Richard L. Sweeney

In the Nonparametrics lecture, we discussed that one prominent application of nonparametric density estimation is auctions. In these problems, we are typically interested in recovering the distribution of values for some object. Once we have these, we can project onto covariates of interest, or simulate an alternative mechanism. 

Guerre, Perrigne and Voung (ECMA 2000) show that the distibution of (unobserved) bidder values can be recovered nonparametrically if we know how the auction works and can estimated the pdf and cdf of bids. 

In this problem, we are going to simulate N auctions, and then try to recover the (known) distribution of valuations using GPV. It turns out that except for a couple simple cases, such as when $v$ is uniformly distributed, even generating the resulting equilibrium bid function is non-trivial. 

Fortunately, [Hickman (IJIO, 2010)](https://www.sciencedirect.com/science/article/pii/S0167718709001076) provides code to simulate bid functions from several distributions in Matlab. This code, available on [his website](http://home.uchicago.edu/hickmanbr/research.html), has been placed in this directory. 

