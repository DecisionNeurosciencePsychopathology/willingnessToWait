function  [ gx ] = g_wtwsceptic_null(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT

% %temperature value
% beta = exp(phi(1));
% 
% gamma = phi(2)/10;

%pull variables to set up gaussians
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector


v_func = sum(v); %subjective value by timestep as a sum of all basis functions

% add cumsum
v_func = cumsum(v_func);
 

%choice rule
% p_choice = (exp((v_func-max(v_func)))) / (sum(exp((v_func-max(v_func))))); %Divide by temperature
%p_choice = v_func;
p_choice = repmat(1/ntimesteps,1,ntimesteps);

%Make the model pick one radom time point
% % x=round(1+rand(1)*(200-1));
% % gx = zeros(1,length(v_func));
% % gx(x)=1;

rt_prev = u(1); %% retrieve previous RT

%% OUTPUT
gx = p_choice';
end
