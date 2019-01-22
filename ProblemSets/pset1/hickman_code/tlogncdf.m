function y=tlogncdf(x,l,u,mu,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%TLOGNCDF(x,l,u,mu,sigma) calculates a truncated lognormal cdf with parameters 
%mu and sigma, lower bound l and upper bound u.  The minimal number
%of arguments is 3; if mu/sigma are left unspecified, then their default values 
%of 0 and 1 will be assigned.  Components of x that lie outside the interval
%[l,u] will be assigned functional values of zero or one.
% 
%See also TNORMCDF, TRAYLCDF, TWBLCDF, TLOGNPDF, TNORMPDF, TRAYLPDF, 
%TWBLPDF

%First, make sure that bounds have been properly specified:
if nargin < 3
    error('YOU MUST AT LEAST SPECIFY AN INPUT AND LOWER/UPPER BOUNDS!');
end
%Next, return a warning message if the lower bound is set to be less than
%the natural lower bound of 0:
if l<0
    display('WARNING: a negative lower bound on the Lognormal is trivially satisfied');
    display('Lower bound has been redefined as "0".');
    l=0;
end

if nargin < 4
    mu = 0;
end
if nargin < 5
    sigma = 1;
end

if l==0 & isinf(u)          % if the bounds supplied are the natural bounds
    y=logncdf(x,mu,sigma)
elseif l==0 & isfinite(u)   % if the right tail only is truncated
    y=logncdf(x,mu,sigma)*(1/logncdf(u,mu,sigma));
elseif l>0 & isinf(u)       % if the left tail only is truncated
    y=(logncdf(x,mu,sigma) - logncdf(l,mu,sigma))*(1/(1-logncdf(l,mu,sigma)));
elseif l>0 & isfinite(u)    %if both tails are truncated
    y=(logncdf(x,mu,sigma) - logncdf(l,mu,sigma))*(1/(logncdf(u,mu,sigma)-logncdf(l,mu,sigma)));
end

y(find(x<=l))=0; % for functional values outside the truncated domain
y(find(x>=u))=1;
