function  [ gx ] = g_sceptic(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT


beta = exp(phi(1)); %Temperature
omicron = 1./(1+exp(-phi(2)))*10; %Decay Opportunity cost 
%omicron = phi(2)/100; %Decay Opportunity cost 


gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;


%Get an approximate integral of the local value basis fx
v_integral=cumtrapz(x_t(1:nbasis));

%Multiply hidden state basis by guassmat to tranform into time domain
%v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
v=v_integral*ones(1,ntimesteps).* gaussmat; %use vector outer product to replicate weight vector
choice = x_t(nbasis+1:nbasis*2)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector for choice

%Sum by timestep
v_func = sum(v); %subjective value by timestep as a sum of all basis functions
choice_func = sum(choice); %subjective choice value by timestep as a sum of all basis functions

%Take mean of each hidden state
mu_v = mean(v_func);
mu_choice = mean(choice_func);

%Calc opportunity cost
% reward_rate = mu_v./mu_choice;
% if isnan(reward_rate) || isinf(reward_rate)
%     reward_rate = 0;
% end
% opportunity_cost = reward_rate.*(1:ntimesteps);
%Calc final value 
%v_func = v_func - 1./(1+opportunity_cost.*exp(omicron));




opportunity_cost = mu_v./mu_choice;
if isnan(opportunity_cost) || isinf(opportunity_cost)
    opportunity_cost = 0;
end

v_func = v_func - opportunity_cost.*omicron;

%Non-sutocorrelated softmax
p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

rt_prev = u(1); %% retrieve previous RT

if strcmp(inG.autocorrelation,'exponential') || strcmp(inG.autocorrelation,'softmax_multitrial')
    pi =  1./(1+exp(-phi(2))); %% introduce a choice autocorrelation parameter lambda
    chi =  1./(1+exp(-phi(3))); %% control the extent of choice autocorrelation
    %% try a Gaussian chi to enable it to go negative (choice anticorrelation as in Lau & Glimcher 2005 at t-1)
    %chi =  phi(3)./100;
end

if strcmp(inG.autocorrelation,'exponential')
    p_choice = p_choice + chi.*(pi.^(abs((1:ntimesteps) - rt_prev)));  %% incorporate an exponential choice autocorrelation function
    p_choice = p_choice./(sum(p_choice));  %% re-normalize choice probability so that it adds up to 1
elseif strcmp(inG.autocorrelation,'softmax_multitrial') || strcmp(inG.autocorrelation,'softmax_multitrial_smooth')
    %% build a matrix of past rts
 if  u(3)>0
    lambda = pi; %When writing the equations we decided lambda is fine for schonberg but we should change it for AR1 so we did...to pi.
    trial = u(3);
    choice_history = inG.rts(1:trial);
    discounted_choice_history = zeros(size(1:ntimesteps));

    for bin = 1:ntimesteps
        if sum(choice_history(1:trial-1)==bin)>0
            when_occurred = find(choice_history(1:trial-1)==bin);
            last_occurred = when_occurred(end);
            trials_ago = trial - last_occurred;
            discounted_choice_history(bin) = lambda^trials_ago; %% it goes to 41 instead of 40  WHY?
        else
            discounted_choice_history(bin) = 0;
        end
    end
    if strcmp(inG.autocorrelation,'softmax_multitrial_smooth')
         iota =  1./(1+exp(-phi(4))); %% control the extent of choice autocorrelation
        discounted_choice_history = smooth(discounted_choice_history,(10*ntimesteps*iota))';
    end
    p_choice = (exp((v_func-max(v_func)+chi.*max(v_func).*discounted_choice_history)/beta)) / (sum(exp((v_func-max(v_func)+chi.*max(v_func).*discounted_choice_history)/beta))); %Divide by temperature
 end
 
 elseif strcmp(inG.autocorrelation,'softmax_tdf')
    %% Just update p_choice
    chi =  1./(1+exp(-phi(3))); %% control the extent of choice autocorrelation
    p_choice = (exp((v_func-max(v_func)+chi.*max(v_func).*choice_func)/beta)) / (sum(exp((v_func-max(v_func)+chi.*max(v_func).*choice_func)/beta))); %Divide by temperature
end
    gx = p_choice';
end



