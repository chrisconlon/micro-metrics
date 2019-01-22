function [pv,bid,paybid,opponent_dens,opponent_dist] = EAauction(N,vbar,delta,dist,a,b,algorithm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Author: Brent Hickman 03/2008                 %%
%%%          Department of Economics               %%
%%%          University of Iowa                    %%
%%%          brent-hickman@uiowa.edu               %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%EAauction(N,vbar,delta,'dist',a,b,algorithm) calculates the equilibrium bid function 
%of an electronic auction for the case of independent private values 
%distributed on the interval [0,vbar], with N bidders and bid increment  
%delta.  See Hickman (2009) for the underlying theory.
%
%%%%%%%%%%%%%%%%%%%%%%%   INPUTS:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The mandatory inputs are N, vbar, delta and 'dist'. The other three are
%optional.  Although This program explicitly computes the equilibrium of a
%model where private values belong to a compact set, the unbounded case can
%be computed by choosing a large vbar where the difference between the CDF 
%F(vbar) and 1 is near machine precision.
%
%Dist, a and b determine the distribution of private values.  The six 
%options for 'dist' are 'Lognormal', 'Normal', 'Rayleigh', 'Uniform', 
%'Weibull', and 'Power' (note that 'dist' is a string argument).  Arguments
%a and b are distributional parameters.  
%
%For the lognormal case, a is the mean of the underlying normal and b is
%the standard deviation.  If a/b are left unspecified, then b is set
%equal to vbar/6 and a is chosen so that the mean is at vbar/2.  For the 
%normal case, a is the mean and b is the standard deviation.  If a
%is left unspecified, then it is set equal to vbar/2.  If b is left
%unspecified, it is set to vbar/6.  For the Rayleigh case, a is the parameter
%and b is ignored.  If a is left unspecified, then it is chosen so that the 
% mean is vbar/2.  For the Weibull case, a is the scale 
%parameter and b is the shape parameter.  When b=1, the Weibull becomes an 
%Exponential distribution with parameter a^(-1): f(x)=a^(-1)*exp(-x/a).  When a 
%and b are left unspecified, b is set equal to 1 and a is chosen so that 
%the mean is vbar/2.  For the uniform case, a and b are ignored and the
%parameters become 0 and vbar.  The power distribution is just a uniform CDF 
%raised to an exponent of a>=1.  If no value for a is specified, the
%default is 1.1.
%
%The final input parameter, algorithm, has two permissable values: 4 and 6.
%This determines the accuracy of the numerical approximation of the EA
%equilibrium .  A value of 4 means that a fourth-order Runge-Kutta
%algorithm will be used, and a value of 6 means a sixth-order
%RKA.  Either version uses a stepsize, h, chosen on line 146 below.  
%The 4th-order version has an O(h^4) error term and the 
%error on the 6th-order version is O(h^5).  On the other hand, run-time for 
%the latter is typically twice as long, so the default value of algorithm
%is 4.  
%
%The syntax is: 
%               [pv,bid]=EAauction(N,vbar,delta,'dist',a,b,algorithm) or 
%               [pv,bid]=EAauction(N,vbar,delta,'dist',a,b) or 
%               [pv,bid]=EAauction(N,vbar,delta,'dist',a) or
%               [pv,bid]=EAauction(N,vbar,delta,'dist'), depending on
%               whether the user wishes to specify distributional
%               parameters.
%
%%%%%%%%%%%%%%%%%%%%%%%   OUTPUTS:   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%The mandatory output arguments are pv, the domain grid (pv stands for 
%private values) and bid, the functional values at the grid points.  The
%optional output arguments are paybid, which gives the probability of
%paying ones bid conditional on winning; opponent_dens, the density of the 
%second-highest order statistic; and opponent_dist, the distribution of 
%the second-highest order statistic. The syntax with all arguments is:
%   [pv,bid,paybid,opponent_dens,opponent_dist] =
%                                       EAauction(N,vbar,delta,'dist',a,b,algorithm)
%%
%See also FPAUCTION 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  DISTRIBUTIONAL SPECIFICATION:  %%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7
    algorithm = 4;
end

if nargin < 4
    error('At a minimum, you must specify N, vbar, delta and dist');
elseif nargout < 2
    error('Too few output arguments');
end

if ischar(dist)==0
    display('You must specify a distribution over private values.  Type either normal, ');
    display('weibull, uniform, rayleigh, or power IN SINGLE QUOTES.');
    error();

elseif strcmpi(dist,'lognormal')
    F = inline('tlogncdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(tlogncdf(pv,l,u,a,b).^(N-2)).*tlognpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 6
        b = vbar/6;
    end
    if nargin < 5
        a = log(vbar/2) - (b^2)/2;
    end
    
elseif strcmpi(dist,'normal')
    F = inline('tnormcdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(tnormcdf(pv,l,u,a,b).^(N-2)).*tnormpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 6
        b = vbar/6;
    end
    if nargin < 5
        a = vbar/2;
    end
    
elseif strcmpi(dist,'rayleigh')
    F = inline('traylcdf(pv,l,u,a).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(traylcdf(pv,l,u,a).^(N-2)).*traylpdf(pv,l,u,a)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 5
        a = vbar*sqrt(2)/(2*sqrt(pi));
        b = 1;
    end
    
elseif strcmpi(dist,'weibull')
    F = inline('twblcdf(pv,l,u,a,b).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(twblcdf(pv,l,u,a,b).^(N-2)).*twblpdf(pv,l,u,a,b)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    if nargin < 6
        b = vbar/6;
    end
    if nargin < 5
        a = vbar/2;
    end
    
elseif strcmpi(dist,'uniform')
    F = inline('unifcdf(pv,l,u).^(N-1)','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('(N-1)*(unifcdf(pv,l,u).^(N-2)).*unifpdf(pv,l,u)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    a = 0;
    b = vbar;
    
elseif strcmpi(dist,'power')
    if nargin < 5
        a = 1.1;
    elseif a<1
        display('When using a power distribution, a must be weakly greater than 1.')
        error()
    end
    F = inline('unifcdf(pv,l,u).^(a*(N-1))','pv','N','l','u','a','b');  %CDF of the second highest order statistic
    f = inline('a*(N-1)*(unifcdf(pv,l,u).^(a*(N-1)-1)).*unifpdf(pv,l,u)','pv','N','l','u','a','b');  %PDF of the second highest order statistic
    b = 1;
end

% tic;
h = min((vbar-(delta))/1000, .001); %stepsize is chosen so that there will be no fewer than 300 grid points
%Retrieve the initial section of the function from FPauction.m
[pv,bid] = FPauction(N,vbar,dist,a,b,delta);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  RUNGE-KUTTA:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%The following is a Runge-Kutta algorithm for solving the differential 
%equation which defines equilibrium bidding in the EA auction.
%See equation (10) of Hickman 2009.
if algorithm == 4
    while pv(end) < vbar
        x = pv(end);
        y = bid(end);
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%@#@#@#@#@#@#@#@#@#@#@# FOURTH-ORDER @#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
            k1  = h*( x - y )*f(x,N,0,vbar,a,b)/...                             %#@
                    ( F(x,N,0,vbar,a,b) -...                                    %@#
                    F(spline(bid,pv,y-delta),N,0,vbar,a,b) );                   %#@
                                                                                %@#
            k2  = h*( x+h/2 - (y+k1/2) )*f(x+h/2,N,0,vbar,a,b)/...              %#@
                    ( F(x+h/2,N,0,vbar,a,b) -...                                %@#
                    F(spline(bid,pv,(y+k1/2)-delta),N,0,vbar,a,b) );            %#@
                                                                                %@#
            k3  = h*( x+h/2 - (y+k2/2) )*f(x+h/2,N,0,vbar,a,b)/...              %#@
                    ( F(x+h/2,N,0,vbar,a,b) -...                                %@#
                    F(spline(bid,pv,(y+k2/2)-delta),N,0,vbar,a,b) );            %#@
                                                                                %@#
            k4  = h*( x+h - (y+k3) )*f(x+h,N,0,vbar,a,b)/...                    %#@
                    ( F(x+h,N,0,vbar,a,b) -...                                  %@#
                    F(spline(bid,pv,(y+k3)-delta),N,0,vbar,a,b) );              %#@
                                                                                %@#
            pvnew  = x + h;                                                     %#@
            bidnew = y + k1/6 + k2/3 + k3/3 + k4/6;                             %@#
                                                                                %#@
            pv  = [pv; pvnew];                                                  %@#
            bid = [bid; bidnew];                                                %#@
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%@#@#@#@#@#@#@#@#@# END FOURTH-ORDER @#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#    
    end
elseif algorithm ==6    
    while pv(end) < vbar
        x = pv(end);
        y = bid(end);
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%@#@#@#@#@#@#@#@#@#@#@# SIXTH-ORDER #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %The following is a sixth-order Runge-Kutta algorithm for solving the
        %differential equation.    
            k1 = h*( x - y )...
                   *f(x,N,0,vbar,a,b)/...
                   ( F(x,N,0,vbar,a,b) -...
                   F(spline(bid,pv,y-delta),N,0,vbar,a,b) );

            k2 = h*( x + h -...
                   (y+k1) )...
                   *f(x+h,N,0,vbar,a,b)/...
                   ( F(x+h,N,0,vbar,a,b) -...
                   F(spline(bid,pv,(y+k1)-delta),N,0,vbar,a,b) );

            k3 = h*( x + h/2 -...
                   ( y+(3*k1 + k2)/8 ) )...
                   *f(x+h/2,N,0,vbar,a,b)/...
                   ( F(x+h/2,N,0,vbar,a,b) -...
                   F(spline(bid,pv,( y+(3*k1 + k2)/8 )-delta),N,0,vbar,a,b) );

            k4 = h*( x + 2*h/3 -...
                   ( y + (8*k1 + 2*k2 + 8*k3)/27 ) )...
                   *f(x+2*h/3,N,0,vbar,a,b)/...
                   ( F(x+2*h/3,N,0,vbar,a,b) -...
                   F(spline(bid,pv,( y + (8*k1 + 2*k2 + 8*k3)/27 )-delta),N,0,vbar,a,b) );

            k5 = h*( x + (7-21^.5)*h/14 -...
                   ( y + (3*(3*21^.5-7)*k1 - 8*(7-21^.5)*k2 + 48*(7-21^.5)*k3 - 3*(21-21^.5)*k4)/392 ) )...
                   *f(x+(7-21^.5)*h/14,N,0,vbar,a,b)/...
                   ( F(x+(7-21^.5)*h/14,N,0,vbar,a,b) -...
                   F(spline(bid,pv,( y + (3*(3*21^.5-7)*k1 - 8*(7-21^.5)*k2 + 48*(7-21^.5)*k3 - 3*(21-21^.5)*k4)/392 )-delta),N,0,vbar,a,b) );

            k6 = h*( x + (7+21^.5)*h/14 -...
                   ( y + (-5*(231+51*21^.5)*k1 - 40*(7+21^.5)*k2 - 320*(21^.5)*k3 + 3*(21+121*21^.5)*k4 + 392*(6+21^.5)*k5)/1960 ) )...
                   *f(x+(7+21^.5)*h/14,N,0,vbar,a,b)/...
                   ( F(x+(7+21^.5)*h/14,N,0,vbar,a,b) -...
                   F(spline(bid,pv,( y + (-5*(231+51*21^.5)*k1 - 40*(7+21^.5)*k2 - 320*(21^.5)*k3 + 3*(21+121*21^.5)*k4 + 392*(6+21^.5)*k5)/1960 )-delta),N,0,vbar,a,b) );

            k7 = h*( x + h -...
                   ( y + (15*(22+7*21^.5)*k1 + 120*k2 + 40*(7*21^.5-5)*k3 - 63*(3*21^.5-2)*k4 - 14*(49+9*21^.5)*k5 + 70*(7-21^.5)*k6)/180 ) )...
                   *f(x+h,N,0,vbar,a,b)/...
                   ( F(x+h,N,0,vbar,a,b) -...
                   F(spline(bid,pv,( y + (15*(22+7*21^.5)*k1 + 120*k2 + 40*(7*21^.5-5)*k3 - 63*(3*21^.5-2)*k4 - 14*(49+9*21^.5)*k5 + 70*(7-21^.5)*k6)/180 )-delta),N,0,vbar,a,b) );

            pvnew  = x + h;
            bidnew = y + ( 9*k1 + 64*k3 + 49*k5 + 49*k6 + 9*k7 )/180;

            pv  = [pv; pvnew];
            bid = [bid; bidnew];
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
        %%%%@#@#@#@#@#@#@#@#@# END SIXTH-ORDER #@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@
        %%%%#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#@#
    end
else
    display('Input variable algorithm must be either 4 or 6, depending on')
    display('whether a Runge-Kutta algorithm of 4th or 6th order is desired.')
    error()
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  FINAL FUNCTIONAL VALUE:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This final segment ensures that the final domain point in the vector pv is
%exactly vbar.

if pv(end) ~= vbar
    h = vbar - pv(end-1); 
    x = pv(end-1);
    y = bid(end-1);
    k1  = h*( x - y )*f(x,N,0,vbar,a,b)/F(x,N,0,vbar,a,b);
    k2  = h*( x+h/2 - (y+k1/2) )*f(x+h/2,N,0,vbar,a,b)/F(x+h/2,N,0,vbar,a,b);
    k3  = h*( x+h/2 - (y+k2/2) )*f(x+h/2,N,0,vbar,a,b)/F(x+h/2,N,0,vbar,a,b);
    k4  = h*( x+h - (y+k3) )*f(x+h,N,0,vbar,a,b)/F(x+h,N,0,vbar,a,b);
    pv(end)  = x + h;   %This forces the final domain point to be exactly vbar.
    bid(end) = y + k1/6 + k2/3 + k3/3 + k4/6;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%  OTHER OUTPUT ARGUMENTS:  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
paybid = 0;
opponent_dens = 0;
opponent_dist = 0;

if nargout > 2
    paybid=ones(length(find(bid<=delta)),1); 
    temp=pv(find(bid>delta)); 
    paybid = [paybid; (F(temp,N,0,vbar,a,b)-F(spline(bid,pv,bid(find(bid>delta))-delta),N,0,vbar,a,b))./F(temp,N,0,vbar,a,b)];
end

if nargout > 3
    opponent_dens = f(pv,N,0,vbar,a,b);
end

if nargout > 4
    opponent_dist = F(pv,N,0,vbar,a,b);
end


