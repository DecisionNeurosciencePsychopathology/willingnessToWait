function  [ gx ] = g_wtwsceptic_logistic(x_t,P,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - u : [useless]
% - inG :
% OUTPUT
% - gx : p(chosen|x_t) or RT

beta = exp(P(1));
%Need to add discrim as an additional param
discrim = 1./(1+exp(-P(2)));

gamma = exp(P(3));
% gamma=0;
%gamma = exp(phi(2));

iti = 20; 
opp_cost_scaling = 1/iti;


gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
u=x_t(nbasis+1:nbasis*2)*ones(1,ntimesteps) .* gaussmat; %Uncertainty is a function of Kalman uncertainties.

v_func = cumsum(sum(v)); %subjective value by timestep as a sum of all basis functions
u_func = cumsum(sum(u)); %vecotr of uncertainty by timestep

opp_cost = x_t(end).*(1:ntimesteps);
cumulative_reward = v_func;
return_on_policy = cumulative_reward-(gamma*opp_cost_scaling*opp_cost);

%In order to run the multinomial version we still need to compute the
%softmax for both the rt_explore and rt_exploit, then let the sigmoid chose
%an action by giving more weight to a specific softmax.
p_rt_exploit = (exp((return_on_policy-max(return_on_policy))/beta)) / (sum(exp((return_on_policy-max(return_on_policy))/beta))); %Divide by temperature
p_rt_explore = (exp((u_func-max(u_func))/beta)) / (sum(exp((u_func-max(u_func))/beta))); %Divide by temperature

%compared to other models that use a curve over which to choose,
%kalman_uv_logistic computes explore and exploit choices and chooses according to a logistic.
u_final = sum(u_func)/length(u_func);

sigmoid = 1/(1+exp(-discrim.*(u_final - inG.u_threshold))); %Rasch model with tradeoff as difficulty (location) parameter

p_choice_final = (((1 - sigmoid).*p_rt_explore) +  (sigmoid.*p_rt_exploit));



if inG.multinomial
    gx = p_choice_final';
else
    best_rts = find(p_choice_final==max(p_choice_final));
    gx = mean(best_rts);
end