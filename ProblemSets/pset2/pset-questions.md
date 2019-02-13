
## Problem set 2
### Prof. Richard L. Sweeney

### Problem 1: Electricity Prices

Deregulated electricity markets solve a series of interellated, high-frequency auctions to ensure that supply equals demand at all points in time and space. The dataset `nodelmp.csv` contains one year worth "locational marginal prices" (lmp's) for a node in Brighton. [If you're interested, more info [here](https://www.iso-ne.com/participate/support/faq/lmp)].

In this problem you will explore the relationship between these prices and temperature (var: `temp`).

1) Make a binner scatter plot of the relationship between the maximum price each day and the temperature. Does the relationship look linear? 

2) Predict the relationship using a polynomial in temperature, for polynomials varying from degree $p=1$ to $p=10$. Which fits best *in sample*?

3) Repeat the above exercise, using k-fold cross validation to select between polynomials of degree 1 through 10. Use $k=10$. Which $p$ fits best now?

4) Plot the predicted lmp for $p=1,2,p^{cv},10$, where $p^{cv}$ is your preferred degree from above. 

5) Predict lmp using a natural cubic spline in temperature. Use k-fold cross-validation to find the optimal number of knots, ranging from 1 to 10.

6) Finally predict lmp using a lowess/ loess regression (using the default bandwidth/ span). 

Plot your preferred polynomial prediction ($p^{cv}$), your prefered spline prediction, and the lowess prediction. 

7) Now plot the residuals for each of these curves against temperature. Which parts of the distribution does each do better/ worse on?

### Problem 2: Hedonic Regression

Harrison & Rubinfeld ([JEEM, 1978](https://www-sciencedirect-com.proxy.bc.edu/science/article/pii/0095069678900062)) conducted one of the first hedonic property value studies of willingness to pay for clean air. Their sample consisted of 506 Census tracts in Massachussets. Their primary interest was the relationship between [NOx](https://www3.epa.gov/region1/airquality/nox.html) concentration levels and median home values. 

A slightly modified version of the HR data (which has been used extensively due to it's inclusion in R's `MASS` package), has been provided in the file `boston_cl.csv`


 There are 506 observations, representing Boston census tracts.
 Variables in order:   
- CRIM     per capita crime rate by town
- ZN       proportion of residential land zoned for lots over 25,000 sq.ft.
- INDUS    proportion of non-retail business acres per town
- CHAS     Charles River dummy variable (= 1 if tract bounds river; 0 otherwise)
- NOX      nitric oxides concentration (parts per 10 million)
- RM       average number of rooms per dwelling
- AGE      proportion of owner-occupied units built prior to 1940
- DIS      weighted distances to five Boston employment centres
- RAD      index of accessibility to radial highways
- TAX      full-value property-tax rate per \$10,000
- PTRATIO  pupil-teacher ratio by town
- BLACK        $1000(Bk - 0.63)^2$ where Bk is the proportion of African-Americans by town
- LSTAT    % lower status of the population
- MEDV     Median value of owner-occupied homes in $1000's (topcoded at 50K)

Set aside 20 percent of the data as a test sample. For each of the models specified below, estimate them using only the remaining 80 percent of the data (training data). 

1) How correllated are these variables?

2) Estimate the original HR model using the training data. Project the the log(Median House Price) onto all of the other variables. Everything should enter linearly, except for `NOx` and `RM`, which should only enter quadratically. 

3) Now estimate the model using LASSO. Use k=10 fold cross validation to select lambda. Then, select the model with the largest lambda (penalty) such that the MSE is within one standard error of the minimum MSE.

3) Do the same thing for ridge regression. Then, select the model with the largest lambda (penalty) such that the MSE is within one standard error of the minimum MSE. 

4) HR's decision to have only `NOx` and `RM` enter quadractically seems sort of arbitrary. Expand the data to contain the square term of all variables. Then run Lasso on this expanded data set. Which coefficients survive now? 

5) Report the internal MSE and test data MSE for HR's original model; Lasso and Ridge on the original covariates; and Lasso on the full set of second order terms. Which one fits best? 
