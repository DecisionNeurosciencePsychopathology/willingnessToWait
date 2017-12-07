%test null agent that chooses a time point randomly (uniformly distributed)

function [wtw]=null_policy(ntimesteps)

wtw=randi([1,ntimesteps],1);

end