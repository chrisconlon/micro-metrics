function [pv,bid,opponent_dens,opponent_dist] = FPauction(N,vbar,dist,a,b,bidstop)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 04/2009                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%FPauction(N,vbar,'dist',a,b,bidstop) calculates the equilibrium 
%bid function of a Dutch or first-price auction for the case of independent 
%private values distributed on the interval [0,vbar] with N bidders.  
%
%%%%%%%%%%%%%%%%%%%%%%%   INPUTS:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The mandatory inputs are N, vbar, delta and 'dist'. The other two are
%optional.  Although This program explicitly computes the equilibrium of a
%model where private values belong to a compact set, the unbounded case can
%be computed by choosing a large vbar where the difference between the CDF 
%F(vbar) and 1 is near machine precision.
%
%Dist, a and b determine the distribution of private values.  The five 
%options for 'dist' are 'Lognormal', 'Normal', 'Rayleigh', 'Uniform', 
%'Weibull', and 'Power' (note that 'dist' is a string argument).  Arguments
%a and b are distributional parameters.   
%
%For the lognormal case, a is the mean of the underlying normal and b is
%the standard deviation.  If a/b are left unspecified, then b is set
%equal to vbar/6 and a is chosen so that the mean is at vbar/2.  For the 
%normal case, a is the mean and b is the standard deviation.  If a
%is left unspecified, then it is set equal to vbar/2.  If b is left
%unspecified, it is set to vbar/6.  For the Rayleigh case, a is the only 
%distribution parameter and b is ignored.  If a is left unspecified, then 
%it is chosen so that the mean is vbar/2.  For the Weibull case, a is the 
%scale parameter and b is the shape parameter.  When b=1, the Weibull 
%becomes an Exponential distribution with parameter a^(-1): 
%a^(-1)*exp(-x/a).  When a and b are left unspecified, b is set equal to 1 
%and a is chosen so that the mean is vbar/2.  For the uniform case, a and b 
%are ignored and the parameters become 0 and vbar.  The power distribution is 
%just a uniform CDF raised to an exponent of a>=1.  If no value for a is 
%specified, the default is 1.1.
%
%The final input argument, bidstop, is optional.  It is included to allow 
%the user to specify an early cutoff point for computation.  If the user 
%does not need all functional values between 0 and the image of vbar, he/she 
%may specify that only the domain points with functional values up to 
%bidstop need be evaluated.  The program ensures that the final 
%gridpoint is either vbar or the inverse of bidstop, which ever applies.
%
%The syntax is: [pv,bid]=FPAuction(N,vbar,'dist',a,b,bidstop) or 
%               [pv,bid]=FPAuction(N,vbar,'dist',a,b) or
%               [pv,bid]=FPAuction(N,vbar,'dist',a) or
%               [pv,bid]=FPAuction(N,vbar,'dist'), depending on
%               whether the user wishes to specify distributional
%               parameters or an early termination value.
%
%%%%%%%%%%%%%%%%%%%%%%%   OUTPUTS:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The mandatory output arguments are pv, the domain grid (pv stands for 
%private values) and bid, the functional values at the grid points.  The
%optional output arguments are opponent_dens, the density of the second-
%highest order statistic; and opponent_dist, the distribution of the 
%second-highest order statistic.
%
%The syntax with all arguments is:
%   [pv,bid,opponent_dens,opponent_dist] =
%                                       FPauction(N,vbar,dist,a,b,bidstop)
%
%In solving for the function, I use a Simpson's Rule algorithm. The
%stepsize h between gridpoints is chosen on lines 134-138 below and the
%algorithm has an O(h^5) approximation error.
%
%See also EAAUCTION

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  DISTRIBUTIONAL SPECIFICATION:  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 3
    error('At a minimum, you must specify N, vbar and dist');
elseif nargout < 2
    error('Too few output arguments');
end

if ischar(dist)==0
    display('You must specify a distribution over private values.  Type either normal, ');
    display('weibull, or uniform IN SINGLE QUOTES.');
    error();

elseif strcmpi(dist,'lognormal')
    F = inline('tlogncdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(tlogncdf(pv,l,u,a,b).^(N-2)).*tlognpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 5
        b = vbar/6;
    end
    if nargin < 4
        a = log(vbar/2) - (b^2)/2;
    end
    
elseif strcmpi(dist,'normal')
    F = inline('tnormcdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(tnormcdf(pv,l,u,a,b).^(N-2)).*tnormpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 5
        b = vbar/6;
    end
    if nargin < 4
        a = vbar/2;
    end
    
elseif strcmpi(dist,'rayleigh')
    F = inline('traylcdf(pv,l,u,a).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(traylcdf(pv,l,u,a).^(N-2)).*traylpdf(pv,l,u,a)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 4
        a = vbar*sqrt(2)/(2*sqrt(pi));
        b = 1;
    end
    
elseif strcmpi(dist,'weibull')
    F = inline('twblcdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(twblcdf(pv,l,u,a,b).^(N-2)).*twblpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 5
        b = vbar/6;
    end
    if nargin < 4
        a = vbar/2;
    end
    
elseif strcmpi(dist,'uniform')
    F = inline('unifcdf(pv,l,u).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(unifcdf(pv,l,u).^(N-2)).*unifpdf(pv,l,u)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    a = 0;
    b = vbar;
    
elseif strcmpi(dist,'power')
    if nargin < 4
        a = 1.1;
    elseif a<1
        display('When using a power distribution, a must be weakly greater than 1.')
        error()
    end
    F = inline('unifcdf(pv,l,u).^(a*(N-1))','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('a*(N-1)*(unifcdf(pv,l,u).^(a*(N-1)-1)).*unifpdf(pv,l,u)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    b = 1;
end

if nargin < 6
    bidstop = inf;
end

if isfinite(bidstop)
    h = min((vbar)/1000,.001); %stepsize is chosen so that there will be no fewer than 300 grid points
else
    h = min([bidstop/100, (vbar)/1000, .001]);
end

pv  = 0;  %Here, I initialize the private value vector and the bid vector
bid = 0; 
gridsize = 151;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  SIMPSON'S RULE:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while ( (pv(end)< vbar) && (bid(end)<bidstop) )
        x = pv(end) + h;
        k = x/(gridsize-1);
        while k>h  % This loop insures that the Simpson step is never larger than h
            gridsize = gridsize+50;
            k = x/(gridsize-1);
        end
        endpts     = [0:2*k:x]';
        midpts     = [0+k:2*k:x]';
        bidnew     = (k/3)*( F(endpts(1),N,0,vbar,a,b) + F(endpts(end),N,0,vbar,a,b) ) +...
                     (2*k/3)*sum( F(endpts(2:end-1),N,0,vbar,a,b) ) +...
                     (4*k/3)*sum( F(midpts,N,0,vbar,a,b) )   ;
        pv  = [pv; x];
        bid = [bid; x - bidnew*( 1/F(x,N,0,vbar,a,b) )];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  FINAL FUNCTIONAL VALUE:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This final segment ensures that the final domain point in the vector pv is
%exactly either vbar or the inverse of bidstop, whichever applies.

if nargin < 6 %i.e.; if no premature cutoff was specified
    h = vbar - pv(end-1); 
    x = pv(end-1);
    x = pv(end) + h;
    k = x/(gridsize-1);
    while  k > min([bidstop/50, (vbar)/300, .025]) % This loop insures that the Simpson step is never larger than the original choice of h
        gridsize = gridsize+50;
        k = x/(gridsize-1);
    end
    endpts     = [0:2*k:x]';
    midpts     = [0+k:2*k:x]';
    bidnew     = (k/3)*( F(endpts(1),N,0,vbar,a,b) + F(endpts(end),N,0,vbar,a,b) ) +...
                 (2*k/3)*sum( F(endpts(2:end-1),N,0,vbar,a,b) ) +...
                 (4*k/3)*sum( F(midpts,N,0,vbar,a,b) )   ;
    pv(end) = x;
    bid(end) = x - bidnew*( 1/F(x,N,0,vbar,a,b) );
else
    if bid(end)<bidstop  %i.e. if bidstop was never reached before the entire domain was traversed
        h = vbar - pv(end-1)
        x = pv(end-1);
        x = pv(end) + h;
        k = x/(gridsize-1);
        while k > min([bidstop/50, (vbar)/300, .025]) % This loop insures that the Simpson step is never larger than the original choice of h
            gridsize = gridsize+50;
            k = x/(gridsize-1);
        end
        endpts     = [0:2*k:x]';
        midpts     = [0+k:2*k:x]';
        bidnew     = (k/3)*( F(endpts(1),N,0,vbar,a,b) + F(endpts(end),N,0,vbar,a,b) ) +...
                     (2*k/3)*sum( F(endpts(2:end-1),N,0,vbar,a,b) ) +...
                     (4*k/3)*sum( F(midpts,N,0,vbar,a,b) )   ;
        pv(end) = x;
        bid(end) = x - bidnew*( 1/F(x,N,0,vbar,a,b) );
        display('');
        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        display('BIDSTOP not reached before traversing domain [0,vbar].');
        display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%');
        Max_func_val = bid(end)
    else
       y = bidstop;
       x = spline(bid,pv,y); %this finds the inverse of bidstop using MATLAB's built-in cubic spline interpolant utility
       bid = [bid(find( bid<bidstop )); y];
       pv  = [pv(find( bid<bidstop )); x];
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  OTHER OUTPUT ARGUMENTS:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
opponent_dens = 0;
opponent_dist = 0;

if nargout > 2
    opponent_dens = f(pv,N,0,vbar,a,b);
end

if nargout == 4
    opponent_dist = F(pv,N,0,vbar,a,b);
end