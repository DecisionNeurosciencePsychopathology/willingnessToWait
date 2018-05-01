function  [ gx ] = g_sceptic(x_t,phi,u,inG)
% INPUT
% - x : Q-values (2x1)
% - beta : temperature (1x1)
% - OR K: mean response tendency
% - inG : multinomial
% OUTPUT
% - gx : p(chosen|x_t) or RT

%temperature value
beta = exp(phi(1));

%pull variables to set up gaussians
gaussmat=inG.gaussmat;
ntimesteps = inG.ntimesteps;
nbasis = inG.nbasis;

v=x_t(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector


v_func = sum(v); %subjective value by timestep as a sum of all basis functions


% add cumsum
v_func = cumsum(v_func);

p_choice = (exp((v_func-max(v_func))/beta)) / (sum(exp((v_func-max(v_func))/beta))); %Divide by temperature

% add opportunity cost 

rt_prev = u(1); %% retrieve previous RT

%% OUTPUT
gx = p_choice';
end
