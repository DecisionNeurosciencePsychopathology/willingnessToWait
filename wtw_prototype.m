%Author: Jonathan Wilson & Alex Dombrovski
%Date last modified: 10/8/14
%Matlab Version: R2012a
%Model of RT on the clock task using Gaussian stimuli for actions or "CSs"
%free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
%learning rule - modified RW with TD(lambda)-type credit assignment
%choice rule: stochastic continuous RT weighted by value
%calls on RewFunction.m (Frank lab code)

%inputs
%params = [alpha, lambda, gamma]
%distrib_num = 1 - uniform, 2 - gen. Pareto, 3 - beta, 4 - discrete1
%agent = null, softmax, e_greedy, uncertainty, logistic

function [cost,constr,v_func,value_hist, wtw, mov, ret] = wtw_prototype(params,distrib_num, agent, trial_length, min_trials, timeout, large_rew, small_rew,sampling_reward)
if (~exist('prr', 'var')), load prr; end
%if (~exist('pev', 'var')), load pev; end;

number_of_stimuli = 12;
trial_plots = 1;
%agent = 'fixedLR_softmax';
beta=0.5; %Softmax temperture



if nargin < 4, trial_length = 200; end %200 100ms bins = 20s
if nargin < 5, min_trials = 100; end % how many trials the agent will sample if he always waits
if nargin < 6, timeout = 20; end % timeout after an immediate quit (2s - 20 10ms bins)

task_length = min_trials*(trial_length+timeout);
maxtrials = ceil(task_length./timeout);          %maximum number of trials in a run
ntimesteps = trial_length; %number of timesteps in action space (set based on pseudorandom sampling of contingency above)

if nargin < 7, large_rew = 10; end % the larger reward for waiting
if nargin < 8, small_rew = 1; end % the smaller reward for quitting
if nargin < 9, sampling_reward = 0.01; end % scaling constant representing the update to the
%uncertainty function in arbitrary units of uncertainty.  Scaled to make
%the sigmoid choice rule work.


%% initialize the movie
mov=repmat(struct('cdata', [], 'colormap', []), maxtrials,1);

%%
pars = {{0, 1},{.4, .575, 0},{0.25, 0.25}, {5, 1}, {maxtrials,300}};
conds = {'unif', 'gp', 'beta','discrete uniform'};
output_names = {'unif', 'gp', 'beta', 'discrete1', 'camel'};
%distrib_num = 4; % 1 - uniform, 2 - gen. Pareto, 3 - early beta, 4 - late
%beta, 5 - piecewise normal

%assign each based on which distribution
cond = output_names{distrib_num};
distrib_pars = pars{distrib_num};
distrib = conds{distrib_num};

%set up seed
softmax_seed=21; %hardcoded
softmax_stream = RandStream('mt19937ar','Seed',softmax_seed);


constr = [];

% if strcmpi(cond, 'unif') %compare string (case insensitive)
%     reward_times = prr.unif*ntimesteps;
%     %     ev_times = pev.unif*ntimesteps;
% elseif strcmpi(cond, 'gp')
%     reward_times = prr.gp*ntimesteps;
%     %     ev_times = pev.gp*ntimesteps;
% elseif strcmpi(cond, 'beta')
%     reward_times = prr.beta*ntimesteps;
%     %     ev_times = pev.beta*ntimesteps;
% elseif strcmpi(cond, 'late_beta')
%     reward_times = prr.late_beta*ntimesteps;
%     %     ev_times = pev.late_beta*ntimesteps;
% elseif strcmpi(cond, 'camel')
%     reward_times = prr.camel*ntimesteps;
%     distrib = fitdist(prr.camel', 'kernel'); %Create PD object for EV computation
% end

%find rewardtimes using the drawSample function
for i = 1 : maxtrials
   reward_times(i) =  drawSample(cond,[])*10; %multiply by 10 to convert to ms
end

% %create histogram of the reward times
% figure(99)
% histogram(reward_times)


%--program start

%Initialize time step vector and allocate for memory
t=1:ntimesteps;


%Selected time points
step_num = number_of_stimuli - 1;
step = (ntimesteps-1)/step_num;
c = -step:step:ntimesteps+step;
nbasis=length(c);


%% define radial basis
[~, ~, tvec, sig_spread, gaussmat, gaussmat_trunc, refspread] = setup_rbf(ntimesteps, nbasis);

%% free parameters: learning rate (alpha), temporal decay (lambda, as in TD(lambda))
alpha = params(1);
lambda = params(2);
gamma = params(3);
%epsilon = params(3);
% log_k = params(4);
% k = exp(log_k);
% k is the hyperbolic discounting parameter
% % % sigma of Gaussian hard coded for now
% % sig = (150.*ntimesteps)./trial_length;


%% initialize RTs as a vector of zeros with a random initial RT;
wtw = zeros(1,maxtrials);

%Have initial rts be constant
%Willingness to wait (WTW) is the quit time, to which the model pre-commits
%at the start of the trial
%NB - times here are in 100ms (centisec) bins
%wtw(1) = ceil(.5*trial_length/10); <- trial_length/10 = deciseconds
wtw(1) = ceil(.5*trial_length);

% i = trial
% t = time step within trial, in centiseconds (1-500, representing 0-5 seconds)

% get a trial loop going
vh =            zeros(maxtrials, nbasis);     % initialize basis height values
deltah =        zeros(maxtrials, nbasis);     % prediction error assigned to each microstimulus
rh =            zeros(maxtrials, nbasis);     % reward assigned to each microstimulus
eh =            zeros(maxtrials, nbasis);     % eligibility traces for each microstimulus in relation to RT (US)
value_by_h =    zeros(nbasis, ntimesteps);   % value by microstimulus (rows for every microstimulus and columns for time points within trial)
value_hist =    zeros(maxtrials, ntimesteps); % history of value by trial
rew_i =           zeros(1, maxtrials);          % actual reward for each trial
%disc_rew =      zeros(1, maxtrials);          % discounted reward for each trial
uh =            zeros(maxtrials, nbasis);     % uncertainty at each microstimulus
udeltah =       zeros(maxtrials, nbasis);     % uncertainty prediction error for each microstumulus (hour)
sampling_h =    ones(maxtrials, nbasis);     % the amount of sampling assigned to each microstimulus at each trial
u_by_h =        zeros(nbasis, ntimesteps);   % uncertainty by microstimulus (rows for every microstimulus and columns for time points within trial)
u_hist =        zeros(maxtrials,ntimesteps);  % history of uncertainty by trial
ev     =        nan(1,maxtrials);             % expected rewards

%choice rule set up
mu_ij =         nan(maxtrials, nbasis);       %means of Gaussians for Kalman
mu_ij(1,:) =    0; %expected reward on first trial is initialized to 0 for all Gaussians.


reward_rate = zeros(1,maxtrials);
opportunity_cost = zeros(trial_length,maxtrials)';
cumulative_reward_fx = zeros(trial_length,maxtrials)';
return_on_policy = zeros(trial_length,maxtrials)';



    %%%%%%%%%% uncertainty driven %%%%%%%%%%%%%
if strcmpi(agent, 'logistic') || strcmpi(agent, 'uncertainty')
        %assign reward
        switch distrib
            case 'gp'
            case 'unif'
            case 'beta'
            case 'discrete uniform'        
                %assign reward at all times as 1
                reward = ones(1,200);
                %assign 10, 20, 30, 200 sec as a reward of 10
                reward(10) = 10; 
                reward(20) = 10;
                reward(30) = 10;
                reward(200) = 10;    
        end

        sigma_noise = repmat(std(reward)^2, 1, nbasis);

        %As in Frank, initialize estimate of std of each Gaussian to the noise of returns on a sample of the whole contingency.
        %This leads to an effective learning rate of 0.5 since k = sigma_ij / sigma_ij + sigma_noise
        sigma_ij(1,:) = sigma_noise;
    
end

%Set up to run multiple runs for multiple ntrials
% for i = 1:ntrials
i=0;
task_time = 0;
while task_time<task_length
    i=i+1; %while loop indexer starts at 0
    
    %         disp(i)
    if wtw(i)> reward_times(i)
        %% rew = real reward; disc_rew = reward discounted for the wait time
        rew_i(i) = large_rew;
        %disc_rew(i) = large_rew.*(1./(1+reward_times(i).*k));
        % get eligibility traces for each stimulus (h)
        % let's assume eligibility decays in inverse proporation to time
        % elapsed from or remaining until the peak of a stimulus
        %% NB - not sure if temporal generalization can distort the ...
        %contingency when the contingency is discontinuous, e.g. two neighboring bumps
        %eh(i,:) = lambda.^(abs(reward_times(i)-c)+1);
        task_time = task_time+reward_times(i)+timeout;
        sampling_h(i,:) =  sampling_h(i,:)+sampling_reward; %sampling_h = vector of 0.01
        update_time(i) = ceil(reward_times(i));
    else
        rew_i(i) = small_rew;
        %disc_rew(i) = small_rew.*(1./(1+reward_times(i).*k));
        %eh(i,:) = lambda.^(abs(wtw(i)-c)+1);
        task_time = task_time+wtw(i)+timeout;
        idx=find(eh(i,:)==max(eh(i,:)));
        sampling_h(i,idx:end) = (eh(i,idx:end)).*sampling_reward;
        sampling_h(i,1:idx) = max(eh(i,:)).*sampling_reward;
        update_time(i) = wtw(i);
    end
    
    %% eligiblity trace
    % get symmetric eligibility traces for each basis function (temporal generalization)
    % generate a truncated Gaussian basis function centered at the RT and with sigma equal to the free parameter.
    
    %compute gaussian spread function with mu = rts(i) and sigma based on free param prop_spread
    elig = gaussmf(tvec, [sig_spread, update_time(i)]);
    
    %compute sum of area under the curve of the gaussian function
    auc=sum(elig);
    
    %divide gaussian update function by its sum so that AUC=1.0, then rescale to have AUC of a non-truncated basis
    %this ensures that eligibility is 0-1.0 for non-truncated update function, and can exceed 1.0 at the edge.
    %note: this leads to a truncated gaussian update function defined on the interval of interest because AUC
    %will be 1.0 even for a partial Gaussian where part of the distribution falls outside of the interval.
    elig=elig/auc*refspread;
    
    %truncated gaussian eligibility
    %figure(7); plot(tvec, elig);
    
    %compute the intersection of the Gaussian spread function with the truncated Gaussian basis.
    %this is essentially summing the area under the curve of each truncated RBF weighted by the truncated
    %Gaussian spread function.
    e_ij(i,:) = sum(repmat(elig,nbasis,1).*gaussmat_trunc, 2);
    
    
    %% update value
    %1) compute prediction error, scaled by eligibility trace
    %this is the basis-wise update, which does not converge to underlying EV
    %delta_ij(i,:) = e_ij(i,:).*(rew_i(i) - mu_ij(i,:));
    
    v_i=sum(mu_ij(i,:)'*ones(1,ntimesteps) .* gaussmat);
    curv = v_i(ceil(update_time(i))); %scalar value estimate at chosen response
    
    %distribute estimated value at RT according to eligibility (works correctly)
    delta_ij(i,:) = e_ij(i,:).*(rew_i(i) - curv);
    
    
    %Variants of learning rule
    if ismember(agent, {'softmax','e_greedy'})
      mu_ij(i+1,:) = mu_ij(i,:) + alpha.*delta_ij(i,:);
         % update mu for uncertainty
    else %this is for uncertainty
        %Kalman gain
        Q_ij(i,:)=0; %for now may remove in future
        k_ij(i,:) = (sigma_ij(i,:) + Q_ij(i,:))./(sigma_ij(i,:) + Q_ij(i,:) + sigma_noise); 

        % update mu for uncertainty
        mu_ij(i+1,:) = mu_ij(i,:) + k_ij(i,:).*delta_ij(i,:);
        
    end
    
    %compute summed/evaluated value function across all timesteps
    v_jt=mu_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector
    v_func = sum(v_jt); %subjective value by timestep as a sum of all basis functions
    v_it(i+1,:) = v_func; %return value vector
    
    %If we run into indexing issues?
    %if i == ntrials, break; end %do not compute i+1 choice on the final trial (invalid indexing problem)
    
    rh(i,:) = rew_i(i); %document reward history
    
    % find the RT corresponding to exploitative choice (choose randomly if value unknown)
    % NB: we added just a little bit of noise
    
    %compute final value function to use for choice
    v_final = v_func; % just use value curve for choice
    
    
%     %%%%%%%%%% uncertainty driven %%%%%%%%%%%%%
   
    %     %include Q???
    %     Q_ij(i+1,:) = p.omega.*abs(delta_ij(i,:)); %use abs of PE so that any large surprise enhances effective gain.
    %     %Compute the Kalman gains for the current trial (potentially adding process noise)

if strcmpi(agent, 'logistic') || strcmpi(agent, 'uncertainty')
        
        %Update posterior variances on the basis of Kalman gains
        sigma_ij(i+1,:) = (1 - e_ij(i,:).*k_ij(i,:)).*(sigma_ij(i,:));

        %Uncertainty is a function of Kalman uncertainties.
        u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;  
        u_func = sum(u_jt); %vector of uncertainties by timestep


        p.beta = params(2);
        p.tau  = .0001; % tau from Greek ???? -- value, price, cf ???? as trophys in Homer
end

    if strcmpi(agent,'uncertainty')   
        uv_func=p.tau*v_func + (1-p.tau)*u_func; %mix together value and uncertainty according to tau
        v_final = uv_func;
        uv_it(i+1,:) = uv_func;
        
    end 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if strcmpi(agent, 'logistic')
        %compared to other models that use a curve over which to choose (either by softmax or egreedy selection),
        %kalman_uv_logistic computes explore and exploit choices and chooses according to a logistic.
        u = sum(u_func)/length(u_func);

        if u == 0
            rt_explore = ceil(.5*ntimesteps);
        else
            rt_explore = find(u_func==max(u_func), 1); %return position of first max (and add gaussian noise?)
        end

        %what do to about p.discrim and u_threshold??
        p.discrim=0;
        u_threshold=0;
        sigmoid = 1/(1+exp(-p.discrim.*(u - u_threshold))); %Rasch model with tradeoff as difficulty (location) parameter
        
        choice_rand=rand(ntimesteps,1);
        
        %compute hard max of value function alone 
        if i == 1
                rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
        else
            rt_exploit = find(v_func==max(v_func), 1); %only take the first max if there are two identical peaks.
            if rt_exploit > max(tvec), rt_exploit = max(tvec); end
        end
        
        if choice_rand(i) < sigmoid
            %explore according to hardmax u
            wtw(i+1) = rt_explore;
        else
            wtw(i+1) = rt_exploit;
        end
        
        v_final = v_func; %no alterations of value function for logistic
    end
    
    
%     softmax_seed=21; %hardcoded random number, used by softmax agent for weighted sampling of softmax function using randsample
% 
%     
%     %setup random number generator for softmax function if relevant
%     if ismember(agent, {'fixedLR_softmax'})
%         %setup a random number stream to be used by the softmax choice rule to make choices consistent during optimization
%         softmax_stream = RandStream('mt19937ar','Seed',softmax_seed);
%     end
    
    
    if i == 1   %indexing contingency
        reward_rate(i) = 0;
    else
        %update reward rate via delta rule TODO add time out to reward rate
        %calc
        reward_rate(i) = reward_rate(i-1) + alpha*(rew_i(i)/update_time(i) - reward_rate(i-1));
    end
    
    %Opportunity cost = reward rate per trial * trial length
    opportunity_cost(i,:) = reward_rate(i).*(1:trial_length);
    %cumulative_reward_fx(i,:) = cumtrapz(1:trial_length, v_final);
    cumulative_reward_fx(i,:) =  v_final;
   
    
     
    %% normalize
%     
%     opportunity_cost(i,:) = (opportunity_cost(i,:)-min(opportunity_cost(i,:))) ./ (max(opportunity_cost(i,:)-min(opportunity_cost(i,:))));
%     if isnan(opportunity_cost(i,:))==1
%         opportunity_cost(i,:)=0;
%     end
%     cumulative_reward_fx(i,:) = (cumulative_reward_fx(i,:)-min(cumulative_reward_fx(i,:))) ./ (max(cumulative_reward_fx(i,:)-min(cumulative_reward_fx(i,:))));
%     %return_on_policy(i,:) = (return_on_policy(i,:)-min(return_on_policy(i,:))) ./ (max(return_on_policy(i,:)-min(return_on_policy(i,:))));
%     
    
    %%
    
    
    return_on_policy(i,:) = cumulative_reward_fx(i,:)-(gamma*opportunity_cost(i,:));

    
    %if statement if return_on_policy is inf
    if any(return_on_policy(i,:)==inf)
        return_on_policy = (0)*ones(1,200);
    end
    

    
    %% choice rule
    
    %compute choice rule according to agent
    %NB: all other models use a softmax choice rule over the v_final curve.
    %p_choice(i,:) = (exp((v_final-max(v_final))/beta)) / (sum(exp((v_final-max(v_final))/beta))); %Divide by temperature
    p_choice(i,:) = (exp((return_on_policy(i,:)-max(return_on_policy(i,:)))/beta)) / (sum(exp((return_on_policy(i,:)-max(return_on_policy(i,:)))./beta))); 
%     if p_choice(i,:)<0
%         p_choice(i,:)=0;
%     end
    %if (all(v_final==0)), v_final=rand(1, length(v_final)).*1e-6; end; %need small non-zero values to unstick softmax on first trial
    
    
    
    % different exploration policies
    if strcmpi(agent,'null')
        wtw(i+1) = null_policy(ntimesteps);
    elseif strcmpi(agent,'e_greedy')
        wtw(i+1) = e_greedy_policy(ntimesteps, return_on_policy(i,:));
    elseif strcmpi(agent, 'softmax') || strcmpi(agent, 'uncertainty')
        %wtw(i+1) = softmax_policy(tvec, p_choice(i,:));
        wtw(i+1) = randsample(softmax_stream, tvec, 1, true, p_choice(i,:));
    end
    
    %populate v_it for tracking final value function
    vfinal_it(i+1,:) = v_final; %store choice function for return according to model
    
    verbose=0;
    if verbose == 1
        fprintf('Trial: %d, Rew(i): %.2f, Rt(i): %.2f\n', i, rew_i(i), wts(i));
        %fprintf('w_i,k:    '); fprintf('%.2f ', mu_ij(i,:)); fprintf('\n');
        %fprintf('delta_ij:   '); fprintf('%.2f ', delta_ij(i,:)); fprintf('\n');
        %fprintf('w_i+1,k:  '); fprintf('%.2f ', mu_ij(i+1,:)); fprintf('\n');
        fprintf('\n');
    end
    
    
    %Oppertunity cost code integrate after initial proto is up and running 
% %     % 10/27/14 Updated choice Rule intergration/analytical method
% %     if i == 1   %indexing contingency
% %         reward_rate(i) = 0;
% %     else
% %         %update reward rate via delta rule
% %         reward_rate(i) = reward_rate(i-1) + alpha*(rew_i(i)/reward_times(i) - reward_rate(i-1));
% %     end
% %     
% %     %Opportunity cost = reward rate per trial * trial length
% %     opportunity_cost(i,:) = reward_rate(i).*(1:trial_length);
% %     cumulative_reward_fx(i,:) = cumtrapz(1:trial_length, value_all);
% %     
% %     
% %     return_on_policy(i,:) = cumulative_reward_fx(i,:)./opportunity_cost(i,:);
% %     
% %     if sum(value_all) == 0 || any(return_on_policy(i,:)==inf)
% %         %rt_exploit = ceil(rand(1)*ntrials); %random number within space
% %         rt_exploit = ceil(.5*trial_length); %default to mid-point of time domain
% %     else
% %         %rt_exploit = max(round(find(value_all==max(value_all))));
% %         rt_exploit = find(return_on_policy(i,:)==max(return_on_policy(i,:)));
% %         if rt_exploit > trial_length
% %             rt_exploit = trial_length;
% %             
% %         elseif rt_exploit < 0 %changed from if to elseif
% %             rt_exploit = 0;
% %         end
% %     end
% %     
% %   
% %     % find the RT corresponding to uncertainty-driven exploration (try random exploration if uncertainty is uniform)
% %     
% %     % u -- total amount of uncertainty on this trial (starts at 0 and decreases)
% %     u = mean(u_all);
% %     if u == 0
% %         rt_explore = ceil(.5*trial_length); %Changed!!! from 5000 previously
% %         
% %     else
% %         rt_explore = max(round(find(u_all(1:trial_length)==max(u_all))));
% %         
% %     end
% %     
% %     
% %     discrim = 100; %need to increase steepness of logistic given the tiny values we have here. Could free later
% %     sigmoid = 1/(1+exp(-discrim*(u - epsilon))); %Rasch model with epsilon as difficulty (location) parameter
% %     
% %     %hard classification of exploration for now at 0.5
% %     if i < maxtrials %do not populate rt on final trial
% %         if sigmoid > 0.5
% %             wtw(i+1) = rt_explore;
% %         else
% %             wtw(i+1) = rt_exploit;
% %         end
% %     end

    %% Compute the expected value of choice for the cost function
    if ~strcmpi(cond, 'camel')
        ev(i+1) = cdf(distrib, wtw(i+1)./trial_length, distrib_pars{:});
    else
        ev(i+1) = cdf(distrib, wtw(i+1)./trial_length);
    end
    
    if trial_plots == 1
        figure(1); %clf;
        
        %% figures for the movie
        subplot(4,2,1:4)
        title('black: wtw blue: RT(reward) red: RT(quit)');  hold on; ...
            plot(find(rew_i(1:maxtrials)==large_rew),update_time(rew_i(1:maxtrials)==large_rew),'bo','LineWidth',2);
        plot(find(rew_i(1:maxtrials)==small_rew),wtw(rew_i(1:maxtrials)==small_rew),'ro','LineWidth',2); hold off;
        
        axis([1 150 1 200])
        
        subplot(4,2,5)
        plot(t,v_func);
        ylabel('value')
        subplot(4,2,6)
        %barh(sigmoid); axis([-.1 1.1 0 2]);
        plot(t(1:ntimesteps),v_jt);
        hold on
        plot(update_time(i),mean(v_jt),'r*')
        ylabel('value temporal basis')
        hold off
        %title(sprintf('trial # = %i', h)); %
        %         xlabel('time(ms)')
        %         ylabel('reward value')
        
        subplot(4,2,7)
        plot(cumulative_reward_fx(i,:));
        ylabel('cumlulative reward function')
        
        subplot(4,2,8)
        plot(return_on_policy(i,:));
        ylabel('return on policy')
        
        drawnow update;
        mov(i) = getframe(gcf);
        
        
        
        %% movie 2
        figure(2)
        plot(1:200,cumulative_reward_fx(i,:),1:200,return_on_policy(i,:),1:200,opportunity_cost(i,:),1:200, v_final,'LineWidth',2.5)
        legend('cdf','return on policy', 'opp cost','v_final')
        %         plot(cumulative_reward_fx(i,:),'r','LineWidth',2); hold on;
%         plot(return_on_policy(i,:),'b','LineWidth',2); hold on;
%         plot(opportunity_cost(i,:),'g','LineWidth',2); hold on;
%         plot(v_final,'m','LineWidth',2);
%         
        drawnow update;
        mov(i)=getframe(gcf);
        
        if strcmpi(agent, 'logistic') || strcmpi(agent, 'uncertainty')
        %% movie 3
        figure(3)
        plot(1:200, u_func, 1:200, v_func)
        legend('u_func','v_func')
        
        drawnow update;
        mov(i)=getframe(gcf);
        end
    end
    %     disp([i rts(i) rew(i) sum(value_all)])
end

%Cost function
cost = -sum(ev(~isnan(ev)));

ret.agent=agent;
ret.nbasis = nbasis;
ret.ntimesteps = ntimesteps;
ret.mu_ij = mu_ij;
ret.delta_ij = delta_ij;
ret.e_ij = e_ij;
ret.v_it = v_it;
ret.vfinal_it = vfinal_it;
ret.rew_i = rew_i;
ret.p_choice=p_choice;
ret.wtw = wtw;
ret.reward_rate = reward_rate;

%plot_wtw_model_data(ret,wtw)