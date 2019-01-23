function y=traylpdf(x,l,u,b)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%TRAYLPDF(x,l,u,b) calculates a truncated Rayleigh pdf with parameter 
%b, lower bound l and upper bound u.  The minimal number
%of arguments is 3; if b is left unspecified, then the default value 
%of 1 will be assigned.  Components of x that lie outside the interval
%[l,u] will be assigned functional values of zero.
% 
%See also TLOGNCDF, TNORMCDF, TRAYLCDF, TWBLCDF, TLOGNPDF, TNORMPDF, 
%TWBLPDF 

%First, make sure that bounds have been properly specified:
if nargin < 3
    error('YOU MUST AT LEAST SPECIFY AN INPUT AND LOWER/UPPER BOUNDS!');
end
%Next, return a warning message if the lower bound is set to be less than
%the natural lower bound of 0:
if l<0
    display('WARNING: a negative lower bound on the Rayleigh is trivially satisfied!');
    display('Lower bound has been redefined as "0".');
    l=0;
end

if nargin < 4
    b = 1;
end

if l==0 & isinf(u)          % if the bounds supplied are the natural bounds
    y=raylpdf(x,b)
    y(find(x<=l))=0;
    y(find(x>=u))=0;
elseif l==0 & isfinite(u)   % if the right tail only is truncated
    y=raylpdf(x,b)*(1/raylcdf(u,b));
    y(find(x<=l))=0;
    y(find(x>u))=0;
elseif l>0 & isinf(u)       % if the left tail only is truncated
    y=raylpdf(x,b)*(1/(1-raylcdf(l,b)));
    y(find(x<l))=0;
    y(find(x>=u))=0;
elseif l>0 & isfinite(u)    % if both tails are truncated
    y=raylpdf(x,b)*(1/(raylcdf(u,b)-raylcdf(l,b)));
    y(find(x<l))=0;
    y(find(x>u))=0;
end

