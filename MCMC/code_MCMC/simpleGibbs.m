%A Gibbs sampler of a Bivariate Normal Distribution
%This is just a toy example, it would be much more efficient to just draw
%from the multivariate normal directly.
clear;
close all;

%Set the Seed for the Random Number generator..
% '0' sets the seed of the clock, any other integer will give you the same
% random number sequence each time.
randn('seed', 0);

%Set parameters
mu = [0 2];
sig = [3 3];
rho = .8;

%Setup space to store draws...
%T = 100;
animate = 1;
if animate
    T = 500;
    afterburn = 50;
else
    T = 100000;
    afterburn = 1000;
end

samp = zeros(T, 2);

%Set initial value (in practice, need a burn in period)
%start = [0 0]
start = [40 40]

figure(1)
myaxis = [-15 50 -15 50];

samp(1,:) = start;
for t=2:T
    %Draw first element conditional on the previous draw...
    samp(t,1) = normrnd(mu(1) + rho*(sig(1)/sig(2))*(samp(t-1,2) - mu(2)), ...
        sqrt(sig(1)^2*(1 - rho^2)));
    
    %Draw second element conditional on first element...
    samp(t,2) = normrnd(mu(2) + rho*(sig(2)/sig(1))*(samp(t,1) - mu(1)), ...
        sqrt(sig(2)^2*(1 - rho^2)));
    if animate
        scatter(samp(1:t, 1), samp(1:t,2));
        axis(myaxis);
        pause(.05);
    end
end

%Compute moments, 
m_samp = mean(samp(afterburn:T,:))
sd_samp = std(samp(afterburn:T,:))
corr_samp = corr(samp(afterburn:T,:))

%Show a picture
figure(1);
scatter(samp(:,1), samp(:,2),'.');
title('Draws from the Gibbs Sampler');

%Compute autocorollation of First Parameter
dsamp = samp(:,1) - m_samp(1);
denom = dsamp'*dsamp;

for k = 1:25
   cov = dsamp(1+k:end)'*dsamp(1:end-k);
   acor(k) = cov/denom;
end
figure(2);
bar(acor)
title('Autocorollation Function of First Parameter');

figure(3)
plot(samp(1:min(T,2000),1))
title('Time Series Trace of First Parameter');

save simpGibbsSamp samp
