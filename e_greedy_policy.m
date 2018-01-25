%epsilon greedy policy

function wtw = e_greedy_policy(ntimesteps, return_on_policy)

epsilon = 0.5; %hardcoded

num=rand; %random number between 0 and 1

if num < epsilon 
    wtw  = randi([1,ntimesteps],1); %pick random time point
else
    [val, indx] = max(return_on_policy);
    wtw = indx; %pick best opinion
end

end