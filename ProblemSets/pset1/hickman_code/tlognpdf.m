function y=tlognpdf(x,l,u,mu,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%TLOGNPDF(x,l,u,mu,sigma) calculates a truncated lognormal pdf with parameters 
%mu and sigma, lower bound l and upper bound u.  The minimal number
%of arguments is 3; if mu/sigma are left unspecified, then their default values 
%of 0 and 1 will be assigned.  Components of x that lie outside the interval
%[l,u] will be assigned functional values of zero.
% 
%See also TLOGNCDF, TNORMCDF, TRAYLCDF, TWBLCDF, TNORMPDF, TRAYLPDF, 
%TWBLPDF 

%First, make sure that bounds have been properly specified:
if nargin < 3
    error('YOU MUST AT LEAST SPECIFY AN INPUT AND LOWER/UPPER BOUNDS!');
end
%Next, return a warning message if the lower bound is set to be less than
%the natural lower bound of 0:
if l<0
    display('WARNING: a negative lower bound on the Lognormal is trivially satisfied!');
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
    y=lognpdf(x,mu,sigma)
    y(find(x<=l))=0;
    y(find(x>=u))=0;
elseif l==0 & isfinite(u)   % if the right tail only is truncated
    y=lognpdf(x,mu,sigma)*(1/logncdf(u,mu,sigma));
    y(find(x<=l))=0;
    y(find(x>u))=0;
elseif l>0 & isinf(u)>0     % if the left tail only is truncated
    y=lognpdf(x,mu,sigma)*(1/(1-logncdf(l,mu,sigma)));
    y(find(x<l))=0;
    y(find(x>=u))=0;
elseif l>0 & isfinite(u)    % if both tails are truncated
    y=lognpdf(x,mu,sigma)*(1/(logncdf(u,mu,sigma)-logncdf(l,mu,sigma)));
    y(find(x<l))=0;
    y(find(x>u))=0;
end

