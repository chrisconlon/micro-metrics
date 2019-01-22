function y=twblpdf(x,l,u,a,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%TWBLPDF(x,l,u,a,b) calculates a truncated Weibull pdf with scale parameter 
%a, shape parameter b, lower bound l and upper bound u.  The minimal number
%of arguments is 3; if a/b are left unspecified, then their default values 
%of 1 will be assigned.  Components of x that lie outside the interval
%[l,u] will be assigned functional values of zero.
%
%When b=1, this function computes a truncated exponential pdf with
%parameter a.
% 
%See also TLOGNCDF, TNORMCDF, TRAYLCDF, TWBLCDF, TLOGNPDF, TNORMPDF, 
%TRAYLPDF

%First, make sure that bounds have been properly specified:
if nargin < 3
    error('YOU MUST AT LEAST SPECIFY AN INPUT AND LOWER/UPPER BOUNDS!');
end
%Next, return a warning message if the lower bound is set to be less than
%the natural lower bound of 0:
if l<0
    display('WARNING: a negative lower bound on the Weibull is trivially satisfied!');
    display('Lower bound has been redefined as "0".');
    l=0;
end

if nargin < 4
    a = 1;
end
if nargin < 5
    b = 1;
end

if l==0 & isinf(u)          % if the bounds supplied are the natural bounds
    y=wblpdf(x,a,b)
elseif l==0 & isfinite(u)   % if the right tail only is truncated
    y=wblpdf(x,a,b)*(1/wblcdf(u,a,b));
elseif l>0 & isinf(u)       % if the left tail only is truncated
    y=wblpdf(x,a,b)*(1/(1-wblcdf(l,a,b)));
elseif l>0 & isfinite(u)    % if both tails are truncated
    y=wblpdf(x,a,b)*(1/(wblcdf(u,a,b)-wblcdf(l,a,b)));
end

y(find(x<l))=0;
y(find(x>u))=0;
