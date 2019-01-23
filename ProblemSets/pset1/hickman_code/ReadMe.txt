These files are what you need to compute the equilibrium of the first-price auction and the static EA auction as in Hickman (2009).  

The two main files are:
FPauction.m:		MATLAB function file that computes the first-price equilibrium by Simpson's rule.
EAauction.m: 		MATLAB function file that computes the static EA auction equilibrium given by equations (2), (10), and (11) 
				in Hickman (2009).  The program uses a 4th-order Runge-Kutta algorithm.

Both of these were designed as MATLAB function files so that they could be easily embedded into other programs.  For usage instructions and syntax, use the MATLAB explorer to navigate to the folder where the two files are located, and type "help FPauction" or "help EAauction" in the MATLAB command window.  There are explanatory annotations at the beginning of each file and this command will display them in the command window.

There are 8 supporting files which compute truncated versions of the pdf and cdf for the lognormal, normal, Rayleigh, and Weibull distributions.  These files are, respectively:
tlognpdf.m
tlogncdf.m
tnormpdf.m
tnormcdf.m
traylpdf.m
traylcdf.m
twblpdf.m
twblcdf.m

These files must be in the same directory as EAauction.m and FPauction.m in order for the latter two to run.  As before, by typing "help file_name" in the command window a breif explanation of the file will be output to the screen.

