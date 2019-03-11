%% Root Finding and Fixed Points
%
%
% There are many many examples from economics we seek to either find the root or the 
% fixed point to a non-linear of (often many) equations, which cannot be computed analytically.
%
% * Many estimation algorithms for equilibrium problems involve a nested
% structure where there is some root-finding problem in the inner nest:
% BLP, dynamic discrete choice (rust, labor models), prodction function, trade models. 
% * You may want to simulate a model and find equilibria using neccessary
% conditions. 
% * Model where equilibria may be defined by a simple threshold of a
% non-linear equation. 

%% Roots
% 
% A function $f$ from $R^n$ to $R^n$ is given and one must compute a vector
% $x$  that satifies: $f(x)=0$.

%% Fixed Points
% 
% A function $g$ from $R^n$ to $R^n$ is given and one must compute a vector
% $x$  that satifies: $x=g(x)$.

%%
% Notice that these two forms are equivilent:
%
% $g(x) = x - f(x)$
%
% *OR*
%
% $f(x) = x - g(x)$
%
% so we will use the same methods to solve both types of problems. 

%%
% <html>
% <br><br><br><br>
% <img src="graph1.png" width="25%"/>
% <br><br><br><br>
% </html>

%% Examples
%
% * Demand and supply market clearing
% * FOCs from an optimization problem (more on this later)

%% Iterative Methods
%
% We will consider methods that *systematically* look over the range of $x$
% until $f(x) = 0$

%% Bisection Method
%
% Intermediate Value Theorem: If a continuous real-valued function assumes two 
% distinct values, then it must assume all values in between. 
%
% If $f$ is continuous and $f(a)$ and $f(b)$ have different signs, then there 
% must be at least one root $x$ in $[a,b]$.
%
% Evaluate $f$ at the bisection of $a$ and $b$. Take the new interval to be the 
% bisected interval with endpoints of different signs. Repeat. 

%%
% <html>
% <br><br><br><br>
% <img src="graph1.png" width="25%"/>
% <br><br><br><br>
% </html>

clear
close all
addpath('~/Dropbox/MATLAB/compecon/CEtools')

%% Bisection Example

f = @(x) x.^3;
a = -6;
b = 12;

tol = 1e-4;
s = sign(f(a)); % sign if the left boundary 
x = (a+b)/2; % inital midpoint
d = (b-a)/2; 
xsave=[];

while d>tol
	d=d/2; % length to cut the next interval
	xsave = [xsave x];
	if s == sign(f(x))
		x = x+d;
	else
		x = x-d;
	end
end
xsave = [xsave x];

fprintf('Solution using user written code: %3.9f\n\n\n',x)

disp(xsave)

%%
% Alternatively, we can use the |bisect| function from the COMPECON toolbox.

x2 = bisect (f,-6,12);  % (function name, a, b)
fprintf('Solution using user bisect code: %3.9f\n',x2)

%%
% Pros and cons of bisection
%
% * Pro: Guaranteed to find a root.
% * Con: Slow (no gradient information).
% * Con: Will only find one root.
% * Con: Only good for single variable functions
% * Con: Can be very slow b/c it does not use info on shape of function

%% Function iteration
%
% * Supply a guess $x^0$
% * Use the updating rule $x^{(t+1)} \leftarrow g(x^(t))$.

%% 
% The starting guess must be close to the fixed point where $||g'(x*)||<1$

%%
% <html>
% <br><br><br><br>
% <img src="graph2.png" width="25%"/>
% <br><br><br><br>
% </html>

%% Function Iteration Example
%
% From the COMPECON Toolbox:

g = @(x) x.^0.5;

xFP1 = fixpoint(g,0.1) % (function name, starting vlaue)
xFP2 = fixpoint(g,1.8) % no start from above the FP

%% Write our own Function Routine
% and capture the results graphically 

vidfile = VideoWriter('testmovie-below.mp4','MPEG-4');
vidfile.FrameRate = 1;
open(vidfile);

xvalues = 0:.1:2;
x_init = 0.1;
xtol = 0.001;
error = 100;
niter=0;
while error>xtol && niter<20
    niter=niter+1;
    x_new = g(x_init(end));        
    error = (x_new - x_init(end)).^2;
    x_init = [x_init x_new];
    
    f = figure('visible','off');
    plot(xvalues,g(xvalues),'LineWidth',2);
    hold on
    plot(xvalues,xvalues,'LineWidth',2);
    hold on    
    plot(x_init,x_init,'*','LineWidth',2)
    hold off
    
%     F(niter) = getframe(gcf); 
    F(niter) = getframe(f);
    
    writeVideo(vidfile, F(niter));
    
    
end
close(vidfile)


%%
% <html>
% <br><br><br><br>
% <video width="640" height="480" controls>
% <source src="testmovie-below.mp4" type="video/mp4">
% </video>
% <br><br><br><br>
% </html>


%% Function Interation Example (from above)
vidfile = VideoWriter('testmovie-above.mp4','MPEG-4');
vidfile.FrameRate = 1;
open(vidfile);

xvalues = 0:.1:4;
x_init = 3.8;
xtol = 0.001;
error = 100;
niter=0;
while error>xtol && niter<20
    niter=niter+1;
    x_new = g(x_init(end));        
    error = (x_new - x_init(end)).^2;
    x_init = [x_init x_new];
    
    f = figure('visible','off');
    plot(xvalues,g(xvalues),'LineWidth',2);
    hold on 
    plot(xvalues,xvalues,'LineWidth',2);
    hold on
    plot(x_init,x_init,'*','LineWidth',2)
    hold off
    
%     F(niter) = getframe(gcf); 
    F(niter) = getframe(f);
    
    writeVideo(vidfile, F(niter));
    
    
end
close(vidfile)


%%
% <html>
% <br><br><br><br>
% <video width="640" height="480" controls>
% <source src="testmovie-above.mp4" type="video/mp4">
% </video>
% <br><br><br><br>
% </html>


%% Newton's Method
%
% * Use derivative information
% * Probably most common method.
% * Sometimes we know the derivative (pen a paper).
% * Sometimes we need to approximate the derivative.
% * Same thing goes with second derivatives. 
%
%

%%
% *The idea:*

%%
% # guess a point
% # linearize the function around that point 
% # find the root of the linear function using taylor expansion
% # use that point as your new guess and repeat

%%
% <html>
% <br><br><br><br>
% <img src="graph3.png" width="25%"/>
% <br><br><br><br>
% </html>


%%
% First-order Taylor approximation: $f(x)\approx f(x^t) + f'(x^t)(x - x^t) = 0$
%
% which yields the following iteration rule: $x^{t+1}\leftarrow x^t - [f'(x^t)]^{-1}f(x^t)$
%
% * *What do you notice about this iterative method?*
% * We need to know the derivative!
% * We will discuss this in detail later.

%%
% *Convergence*: Judd Theorem 2.1 (page 130) -- If $x^1$ is "sufficiently" close
% to $x^*$, $f'(x^*)\ne0$ and $\mid \frac{f''(x^*)}{f'(x^*)}<\infty$, then the Newton
% sequence will converge to $x^*$. Also, $f$ needs to be "smooth."
%
% * Warning: if $f'(x^t)$ is close to zero, then it can overshoot and cause
% problems



%% Newton Example
%
% Simple demand function in a separate file:

%%
%
% <include>simpleFunc.m</include>
% 

%%
xvals = .3:.05:.7;
plot(xvals,simpleFunc(xvals))
hold on;
plot(xvals,zeros(size(xvals)))
hold off;

%% 
% Set options for the COMPECON function 'newton':

optset('newton','maxit',20);optset('newton','showiters',1);

%%
% Now, call the 'newton' routine to find the root of simpleFunc.

[fstar,fval,flag] = newton('simpleFunc',.1)

%%
% What happens if we pick a _weird_ starting vlaue?
xvals = 0:.05:.7;
plot(xvals,simpleFunc(xvals))
hold on;
plot(xvals,zeros(size(xvals)))
hold off;


%%
[fstar,fval,flag] = newton('simpleFunc',0)


%% Quasi-Newton Methods
%
% Many times we do not have an analytical derivative:
%
% * It is difficult to compute analytically.
% * Potentially make mistakes in coding.
% * In general humans make mistakes.

%%
% Quasi-Newton methods are the same as the Newton method, except with an 
% approximation of the jacobian.

%%
% _Secant Method_
%
% Univariate Newton method with Jacobian approximation.
%
% Replace $f'$ with an approximation from the last two function values:
%
% $f'(x^t) \approx \frac{f(x^t) - f(x^{t-1})}{x^t - x^{t-1}}$
%
% which yields the following update rule:
%
% $x^{t+1} \leftarrow x^t - \frac{ x^t - x^{t-1} }{ f(x^t) - f(x^{t-1}) } f(x^t)$
%
% You are constructing the approximating line through the two points $(x^t,f(x^t))$ 
% and $(x^{t-1},f(x^{t-1}))$.

%%
% <html>
% <br><br><br><br><br><br>
% </html>
% 

%%
%  [GRAPH]

%%
% <html>
% <br><br><br><br><br><br>
% </html>


%%
% _Broyden's Method_
%
% Multivariate version of the secant method. 
%
% * Generate a sequence of vectors $x^t$ and matrices $A^t$
% * These approximate the root and Jacobian of $f$
% * Guess $x^0$ and $A^0$.
% * $A^0$ is often set to the numerical jacobian at x^0.

%%
% $f(x) \approx f(x^t) + A^t(x-x^t) = 0$
% which yields the following rule
% $x^{t+1} \leftarrow x^t - (A^t)^{-1}f(x^t)$
%

%%
% The Jacobian is also updated iteratively: 
% $A^{t+1} \leftarrow A^t + [f(x^{t+1}) - f(x^t) -
% A^td^t]\frac{d^t}{d^td^t}$
% where $d^t = x^{t+1} - x^t$
%
% In priactice we will update the inverse of the Jacboian to save an inversion step. 
%
% _NOTE:_ The sequence of approximations of the Jacobian DO NOT neccessarily 
% converge to the true Jacobian. 

%%
% This method will work if you start sufficiently close, and $f$ is well behaved...duh!
% 
% In priactice, I have used this method and it has worked very well for problems
% where the Jacobian diagonally dominant. 

%% Newton Methods in Practice.
% See Gravity Example

!rm *_eq*.png







