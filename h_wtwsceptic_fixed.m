function  [fx] = h_wtwsceptic_fixed(x_t, theta, u, inF)
% evolution function of q-values of a RL agent (2-armed bandit problem)
% [fx,dfdx,dfdP] = f_Qlearn2(x,P,u,in)
% Here, there are only two q-values to evolve, i.e. there are only two
% actions to reinforce (2-armed bandit problem).
% IN:
%   - x_t : basis values/heights (nbasis x 1)
%   - theta : theta(1) = prop_spread; theta(2) = alpha;
%   - u : u(1) = rt; u(2) = reward
%   - inF : struct of input options (has nbasis and ntimesteps)
% OUT:
%   - fx: evolved basis values/heights (nbasis x 1)

alpha = 1./(1+exp(-theta(1)));
%alpha = sig(theta(1));

% alpha2 = 1./(1+exp(-theta(2)));

alpha2 = alpha;


if inF.fit_propspread
    prop_spread = 1./(1+exp(-theta(3))); %0..1 SD of Gaussian eligibility as proportion of interval
    sig_spread=prop_spread*range(inF.tvec); %determine SD of spread function in time units (not proportion)
    
    %if prop_spread is free, then refspread must be recomputed to get AUC of eligibility correct
    refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));
else
    sig_spread = inF.sig_spread; %precomputed sig_spread based on fixed prop_spread (see setup_rbf.m)
    refspread = inF.refspread; %precomputed refspread based on fixed prop_spread
end
rt = u(1);
reward = u(2);
time = u(3);
rt_time_bin = u(4);

iti = 20; 
tau = iti + rt;

if inF.fit_nbasis
    %% convert normally distributed theta(2) to discrete uniform number of bases
    nbasis_cdf = cdf('Normal',theta(3),inF.muTheta2, inF.SigmaTheta2);
    nbasis = unidinv(nbasis_cdf,inF.maxbasis);
else
    nbasis = inF.nbasis;
end
% ntimesteps = inF.ntimesteps;

%refspread = sum(gaussmf(min(inF.tvec)-range(inF.tvec):max(inF.tvec)+range(inF.tvec), [sig_spread, median(inF.tvec)]));

%compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
% elig = gaussmf(inF.tvec, [inF.sig_spread, rt]);
elig = gaussmf(inF.tvec, [sig_spread, rt]);

%compute sum of area under the curve of the gaussian function
auc=sum(elig);

%divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
%this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
%note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
%will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
elig=elig/auc*refspread;

%If r is a 'win' i.e. the subject quit and recived a 1 as reward elig is 0
%after the rt else we have a boxcar from the peak of the elig back to 0
 peak_elig = max(elig);

if inF.use_boxcar_elig
   if reward==1
       elig(rt_time_bin+1:end)=0;
       elig(1:rt_time_bin)=peak_elig;

   else
       %peak_elig = elig(rt_time_bin);
       elig(rt_time_bin+1:end)=0;
   end
end


%If we want to replace the elig trace with a dirac impulse function
%I don't think this is correct
if inF.use_dirac
    elig = zeros(size(elig));
    elig(rt_time_bin)=reward;
end


% figure(88)
% plot(elig)

%compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
%this is essentially summinF the area under the curve of each truncated RBF weighted by the truncated
%Gaussian spread function.
e = sum(repmat(elig,nbasis,1).*inF.gaussmat_trunc, 2);


%value 
value = x_t(1:nbasis);

%initalize fx
fx = zeros(length(x_t),1);

%1) compute prediction error, scaled by eligibility trace
delta = e.*(reward - value);

fx(1:nbasis,:) = value + alpha.*delta;


%if inF.diagnos_model
% % h=figure(99);
% % subplot(2,2,1)
% % plot(fx(1:end-1), 'Linewidth',5)
% % title(sprintf('value of basis function with RT: %.2f', rt))
%end


%add in reward rate as hidden state
if inF.tau_rr 
    fx(end) = ((1-alpha2).^tau)*x_t(end)+(1-(1-alpha2).^tau)*(reward/tau);
else 
    fx(end) = x_t(end) + alpha2*(reward/rt - x_t(end));
end


% % %For diagnosis purposes only
% % gamma = 10^(0.1688 + 1);
% % gaussmat=inF.gaussmat;
% % ntimesteps = inF.ntimesteps;
% % nbasis = inF.nbasis;
% % 
% % v=fx(1:nbasis)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
% % 
% % 
% % v_func = sum(v); %subjective value by timestep as a sum of all basis functions
% % 
% % % add cumsum
% % v_func = cumsum(v_func);
% %  
% % 
% % % add opportunity cost = rr*trial length (max 200)
% % opp_cost = fx(end).*(1:ntimesteps);
% % 
% % cumulative_reward_fx = v_func;
% % 
% % return_on_policy = cumulative_reward_fx-(gamma*opp_cost);
% % fx(nbasis+1:end-1) = return_on_policy;


end



