function [delay,quantile] = drawSample(distrib,quantile)
% generates a sample from next quantile of the designated distribution
% timing parameters are specified within this function
% distrib may be 'gp' or 'unif'
% quantile is a vector (may be empty)

try quantile;
catch 
    quantile=[];
end

nQuants = 4; % number of partitions to sample w/o replacement
if isempty(quantile), quantile = randperm(nQuants); end
q = quantile(1);
quantile(1) = [];

switch distrib
    case 'gp'
        
        params = {4, 5.75, 0}; 
        maxDelay = 20; % longest possible delay (in s) before a trial wins
        maxX = cdf('gp',maxDelay,params{:}); % greatest x value to submit to icdf
        x = (q-1+rand)*maxX/nQuants;
        delay = icdf('gp',x,params{:});
        
    case 'unif'
        
        params = {0, 12};
        maxX = 1; % no values are out of range
        x = (q-1+rand)*maxX/nQuants;
        delay = icdf('unif',x,params{:});
        
    case 'beta'
        
        params = {0.25, 0.25};
        scaling = 12;
        maxX = 1; % no values are out of range
        x = (q-1+rand)*maxX/nQuants;
        delay = scaling * icdf('beta',x,params{:});
        
    case 'discrete1'
        
        vals = [1, 2, 3, 20]; % equally likely delay values (sec)
        delay = vals(q);
        
end


