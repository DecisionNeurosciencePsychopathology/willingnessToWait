pd1 = makedist('Uniform','lower',0,'upper',12);

% Compute the cdfs
x = 0:1:20;
cdf1 = cdf(pd1,x);

% Plot the cdfs
figure;
plot(x,cdf1,'r','LineWidth',2);