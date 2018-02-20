function data = simulate_data(params, distrib_num, agent)

%this function creates test data for the wtw_sceptic_vba 

[cost,constr,v_func,value_hist, wtw, mov, ret] = wtw_prototype(params,distrib_num, agent);

data.latency = wtw; %willingness to wait for a reward
data.trialResult = cell(1,1100);
for i=1:length(ret.rew_i) % win or quit
    if ret.rew_i == 10
        data.trialResult(i) = {'win'};
    else 
        data.trialResult(i) = {'quit'};
    end
end
data.trialResult=data.trialResult';
data.initialTime = wtw;
data.payoff = ret.rew_i; %reward received at each time point

end