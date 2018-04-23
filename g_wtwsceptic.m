function  [ gx ] = g_wtwsceptic(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT

%temperature value
beta = exp(phi(1));

% gamma = phi(2)/10;
gamma = 10^(phi(2) + 1);
%gamma = exp(phi(2));

%pull variables to set up gaussians
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector


v_func = sum(v); %subjective value by timestep as a sum of all basis functions

% add cumsum
v_func = cumsum(v_func);
 

% add opportunity cost = rr*trial length (max 200)
opp_cost = x_t(end).*(1:ntimesteps);

cumulative_reward_fx = v_func;

return_on_policy = cumulative_reward_fx-(gamma*opp_cost);


%choice rule
%p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

p_choice = (exp((return_on_policy-max(return_on_policy))/beta)) / (sum(exp((return_on_policy-max(return_on_policy))/beta))); 


% % h=figure(99);
% % subplot(2,2,2)
% % plot(p_choice,'r','LineWidth',4)
% % title(sprintf('P Choice AKA gx with beta: %d',beta))
% % subplot(2,2,3)
% % plot(return_on_policy,'k','LineWidth',4)
% % title('Return on policy')
% % subplot(2,2,4)
% % plot(opp_cost,'g','LineWidth',4)
% % title(sprintf('Opportunity cost with gamma: %d',gamma))


%% OUTPUT
gx = p_choice';
end
