function y=tnormcdf(x,l,u,mu,sigma)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%TNORMCDF(x,l,u,mu,sigma) calculates a truncated normal cdf with mean mu, 
%standard deviation sigma, lower bound l and upper bound u.  The minimal number
%of arguments is 3; if mu & sigma are left unspecified, then their default values 
%of 0 and 1 will be assigned.  Components of x that lie outside the interval
%[l,u] will be assigned functional values of zero or one.
% 
%See also TLOGNCDF, TRAYLCDF, TWBLCDF, TLOGNPDF, TNORMPDF, TRAYLPDF, 
%TWBLPDF 

%First, make sure that bounds have been properly specified:
if nargin < 3
    error('YOU MUST AT LEAST SPECIFY AN INPUT AND LOWER/UPPER BOUNDS!');
    error('For one-sided trunctation, enter +inf/-inf as one of the bounds.');
end

if l>=u
    error('Lower bound must be strictly less than the upper bound!');
end

if nargin < 4
    mu = 0;
end
if nargin < 5
    sigma = 1;
end

if isinf(l) & isinf(u)          % if the bounds supplied are the natural bounds
    y=normcdf(x,mu,sigma)
elseif isinf(l) & isfinite(u)   % if the right tail only is truncated
    y=normcdf(x,mu,sigma)*(1/normcdf(u,mu,sigma));
elseif isfinite(l) & isinf(u)   % if the left tail only is truncated
    y=(normcdf(x,mu,sigma) - normcdf(l,mu,sigma))*(1/(1-normcdf(l,mu,sigma)));
elseif isfinite(l) & isfinite(u)% if both tails are truncated
    y=(normcdf(x,mu,sigma) - normcdf(l,mu,sigma))*(1/(normcdf(u,mu,sigma)-normcdf(l,mu,sigma)));
end


y(find(x<=l))=0;                % for functional values outside the truncated domain.
y(find(x>=u))=1;
