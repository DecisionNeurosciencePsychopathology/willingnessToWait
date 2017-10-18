%Diagnose the discounted and undiscounted value basis from wtw

gaussmat=options.inG.gaussmat;
ntimesteps = options.inG.ntimesteps;
nbasis = options.inG.nbasis;
omicron = posterior.muPhi(2); %Decay Opportunity cost 

%This will plot the value basis functions per trial
value_history=out.suffStat.muX(1:n_basis,:);
figure(99)
clf;
surf(value_history);

%This will plot the value basis cumulative integral -- irrelevent
value_history_integral=cumtrapz(out.suffStat.muX(1:n_basis,:));
figure(102)
clf;
surf(value_history_integral);

x_t = out.suffStat.muX;
v_integral=cumtrapz(x_t(1:nbasis,:));

%Overwrite omicron
%omicron = .01;

%Start trial loop
for i = 1:length(x_t);
    
    %Multiply hidden state basis by guassmat to tranform into time domain
    %v=x_t(1:n_basis,i)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    v=v_integral(:,i)*ones(1,ntimesteps).* gaussmat; %use vector outer product to replicate weight vector
    choice = x_t(n_basis+1:n_basis*2,i)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector for choice
    
    %Sum by timestep
    v_func(i,:) = sum(v); %subjective value by timestep as a sum of all basis functions
    choice_func(i,:) = sum(choice); %subjective choice value by timestep as a sum of all basis functions
    v_func_undiscounted(i,:) = sum(v);
    
    %Take mean of each hidden state
    mu_v(i,:) = mean(v_func(i,:));
    mu_choice(i,:) = mean(choice_func(i,:));
    
%Calc opportunity cost
% reward_rate(i) = mu_v(i,:)./mu_choice(i,:);
% if isnan(reward_rate(i)) || isinf(reward_rate(i))
%     reward_rate(i) = 0;
% end
% opportunity_cost(i,:) = reward_rate(i).*(1:ntimesteps);
%Calc final value 
%v_func(i,:) = v_func(i,:) - 1./(opportunity_cost(i,:).*exp(omicron));





opportunity_cost(i) = mu_v(i,:)./mu_choice(i,:);
if isnan(opportunity_cost(i)) || isinf(opportunity_cost(i))
    opportunity_cost(i) = 0;
end
v_func(i,:) = v_func(i,:) - opportunity_cost(i).*omicron;
end


figure(100)
clf;
surf(v_func);
title('Discounted value per trial by timestep')

figure(101)
clf;
surf(v_func_undiscounted);
title('Undiscounted value per trial by timestep')

% figure(301)
% clf;
% surf(choice_func)
