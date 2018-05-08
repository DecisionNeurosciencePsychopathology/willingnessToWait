function  [ gx ] = g_wtwsceptic(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT

%temperature (not inverse!)
beta = exp(phi(1));

% opportunity cost sensitivity
gamma = exp(phi(2));

% uncertainty sensitivity first fixed at neutrality
tau = 1;
% gamma=0;
%gamma = exp(phi(2));

iti = 20; 
opp_cost_scaling = 1/iti;
%pull variables to set up gaussians
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;
u = zeros(1,ntimesteps);

%Define value
w = x_t(1:nbasis);
v = w*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector

if inG.kalman.kalman_uv_sum || inG.kalman.fixed_uv || inG.kalman.fixed_decay_uv
%    tau = 1./(1+exp(-phi(3))); %Uncertainty mixing: 0..1
% uncertainty sensitivity
% tau = 1 + phi(3)/1000;
tau = 1./(1+exp(-phi(3)-10)); 

sigma=x_t(nbasis+1:nbasis*2);
   
   %perform tau mixing here
   u = sigma*ones(1,ntimesteps) .* gaussmat;
end
u_func = sum(u);


v_func = sum(v); %subjective value by timestep as a sum of all basis functions

% add cumsum
v_func = cumsum(v_func);

%Scale v_func down to the max of the basis
v_func = v_func./max(v_func)*max(w);

%If we are dividing by 0
if isnan(v_func)
    v_func = zeros(1,ntimesteps);
end

cumulative_reward_fx = v_func;



% incorporate opportunity cost = rr*trial length (max 200)
opp_cost = x_t(end).*(1:ntimesteps);

return_on_policy = cumulative_reward_fx-(gamma.*opp_cost.*opp_cost_scaling);

% uncertainty-seeking only matters past the value bump -- nothing is
% learned by quitting earlier
u_func(1:return_on_policy==max(return_on_policy)) = 0;

%Update with uncertainty if appliciable
return_on_policy = tau .* return_on_policy + (1-tau).*u_func;




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
