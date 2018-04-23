function [posterior,out] = wtw_sceptic_vba(data1,id,sample_to_use,model,n_basis, multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion,tau_rr, saveresults, graphics)



%working models = fixed, null, uncertainty (uv sum), logistic


%% fits SCEPTIC model to Clock Task subject data using VBA toolbox
% example call:
% [posterior,out]=clock_sceptic_vba(10638,'modelname',nbasis,multinomial,multisession,fixed_params_across_runs,fit_propsrpead)
% id:           5-digit subject id in Michael Hallquist's BPD study
% only works with 'fixed' (fixed learning rate SCEPTIC) so far
% n_basis:      8 works well, 4 also OK
% multinomial:  if 1 fits p_chosen from the softmax; continuous RT (multinomial=0) works less well
% multisession: treats runs/conditions as separate, helps fit (do not allow X0 to vary though)
% fixed_params_across_runs -- self-explanatory
% fit_propspread -- makes temporal generalization within the eligibility trace a free parameter
% n_steps:      number of time bins
% u_aversion:   allow for uncertainty (ambiguity) aversion for UV_sum
%%
close all

%% uncertainty aversion for UV_sum
if nargin<11
    u_aversion = 0;
    saveresults = 0; %change later
    graphics = 0;
elseif nargin<12
    saveresults = 0; %change later
    graphics = 0;
elseif nargin<13
    graphics = 0;
end


global rew_rng_state no_gamma
rew_rng_seed = 99;
options.inG.autocorrelation = 'none';

%Do not grpah Volterra kernals
options.plotKernals = 0;

if graphics
    options.DisplayWin = 1;
else 
    options.DisplayWin = 0;
end
%% set up dim defaults
n_theta = 1;
n_phi = 1;

%% fit as multiple runs
% multisession = 1;
% fix parameters across runs
% fixed_params_across_runs = 1;

%% Load in the data struct for now, in theory we should just paa the subject's data to the function at this point instead of constantly loading it from file.
if tau_rr    
    if ~exist('willingnessToWait/tau_rr','dir')
        mkdir('tau_rr')
    end
    
    rr_type = 'tau_rr';
else
    if ~exist('willingnessToWait/normal_rr','dir')
        mkdir('normal_rr')
    end
    
    rr_type = 'normal_rr';
end

results_dir = sprintf('ss_test/%s/%s', rr_type,num2str(id));

if ~exist(results_dir,'dir')
    mkdir(results_dir)
end


%remove row from data if there is nothing in it
i = 1;
for ii = 1:length(data1)
   if isempty(data1(ii).latency) == 0
       data(i) = data1(ii);
       i= i + 1;
   end
end

%Options changing the mechanics of the model evolution function
options.inF.tau_rr = tau_rr; %RossOtto Reward Rate
options.inF.use_boxcar_elig = 1; %If we want to use the boxcar idea to update the eligibility traces
options.inF.use_dirac = 0; %Use the dirac impulse fuciton (will overwrite the boxcar method)


options.inF.fit_nbasis = 0;
range_RT = 200; %The max time until new trial is 20 seconds
n_t = size(~isnan([data.latency]),2);
iti=20; %2 second inter-trial interval
n_runs = n_t/50;
trialsToFit = 1:n_t;
% fit_propspread = 1;
options.inF.fit_propspread = fit_propspread;


%% set up models within evolution/observation Fx
%Note: we might need to add option.inF.model to make the kalman models
%easier to deal with...
options.inF.nbasis = n_basis;
options.inF.ntimesteps = n_steps;
options.inG.ntimesteps = n_steps;
options.inG.multinomial = multinomial;
options.inG.nbasis = n_basis;
options.inG.maxRT = range_RT;
%%
options.TolFun = 1e-6;
options.GnTolFun = 1e-6;
options.verbose=1;
% options.DisplayWin=1;

%% set up kalman defaults
options.inF.kalman.processnoise = 0;
options.inF.kalman.kalman_sigmavolatility  = 0;
options.inF.kalman.kalman_softmax = 0;
options.inF.kalman.kalman_logistic = 0;
options.inF.kalman.kalman_uv_logistic = 0;
options.inF.kalman.kalman_uv_sum = 0;
options.inF.kalman.kalman_uv_sum_sig_vol = 0;
options.inF.kalman.fixed_uv = 0;
options.inF.kalman.kalman_sigmavolatility_local =0;
options.inF.kalman.kalman_sigmavolatility_precision=0;


%% set up basis
[~, ~, options.inF.tvec, options.inF.sig_spread, options.inG.gaussmat, options.inF.gaussmat_trunc, options.inF.refspread] = setup_rbf(options.inF.ntimesteps, options.inF.nbasis, .08);


%Set up sigma noise for every point in u or hidden state?
rng(rew_rng_seed); %inside trial loop, use random number generator to draw probabilistic outcomes using RewFunction
rew_rng_state=rng;

%From evalDisctrete we get te hiRes policy payoff for quitting at time t in
%.01 bins. We'll use this as the variance of returns from typical runs for
%the Kalman gain.
load('policyPay_hiRes.mat')
sigma_noise = policyPay_hiRes;
sigma_noise = mean(sigma_noise);
options.inF.sigma_noise = sigma_noise;
options.inF.gaussmat = options.inG.gaussmat;

%% split into conditions/runs
if multisession %improves fits moderately
    options.multisession.split = repmat(n_t/n_runs,1,n_runs); % two sessions of 120 datapoints each
    %% fix parameters
    if fixed_params_across_runs
        options.multisession.fixed.theta = 'all';
        options.multisession.fixed.phi = 'all';
        %
        % allow unique initial values for each run?x
        options.multisession.fixed.X0 = 'all';
    end
    
end



%Determine which evolution funciton to use
options.inF.kalman.(model)=1; %Declare which model to use if kalman

switch model
    %fixed learning rate (alpha) for PE+ and PE-; softmax choice rule
    case 'fixed'
        h_name = @h_wtwsceptic_fixed;
        hidden_variables = 1; %tracks only value
        priors.muX0 = zeros(hidden_variables*n_basis+1,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis+1); % added 1 for rr
        n_phi = 2; % one for return on policy, one for choice rule
        n_theta = 1;


    case 'null' %no updates
        h_name = @h_wtwsceptic_null;
        hidden_variables = 1;
        n_phi = 0;
        n_theta = 0;
        priors.muX0 = zeros(hidden_variables*n_basis+1,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis+1); % added 1 for rr
%     case 'fixed_opportunity_cost'
%         h_name = @h_sceptic_fixed_opportunity_cost;
%         hidden_variables = 2; %tracks only value
%         priors.muX0 = zeros(hidden_variables*n_basis,1);
%         priors.SigmaX0 = zeros(hidden_variables*n_basis);
%         n_theta = 2;
%         n_phi = 3- strcmp(options.inG.autocorrelation, 'none');
% %         priors.SigmaX0 = 10*ones(hidden_variables*n_basis);
    case 'fixed_decay'
        h_name = @h_wtwsceptic_fixed_decay;
        hidden_variables = 1; %tracks only value
        priors.muX0 = zeros(hidden_variables*n_basis+1,1);
        priors.SigmaX0 = zeros(hidden_variables*n_basis+1);
        n_theta = 2; %learning rate and decay outside of the eligibility trace
        n_phi = 2; % one for return on policy, one for choice rule

        
        %kalman learning rule (no free parameter); softmax choice over value curve
    case 'kalman_softmax'
        %Prop_spread is the only variable in this model
        if fit_propspread
            n_theta = 0;
        end
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis+1,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis +1);
        h_name = @h_sceptic_kalman;
        
        %old kalman with explore/exploit hardmax selection according to logistic function
    case 'kalman_logistic'
        %Prop_spread is the only variable currently in this model
        if fit_propspread
            n_theta = 0;
        end
        n_phi  = 3;  %Beta and discrim and gamma (for return on policy)
        n_theta = 2;
        %Different observation function than other kalman models
%         g_name = @g_sceptic_logistic;
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis+1,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis+1); %This is Discrim not Beta for this model
        h_name = @h_wtwsceptic_kalman;
        %Define indifference point between explore and exploit (p = 0.5) as proportion reduction in variance from initial value
        tradeoff = 0.1209529; %From what was the optimized overall
        options.inG.u_threshold = (1 - tradeoff * sigma_noise);
        %Predetermined random trials
        %options.inG.choice_rand=rand(n_steps,1);
        
        %kalman learning rule and uncertainty update; V and U are mixed by tau; softmax choice over U+V
    case 'kalman_uv_sum'
        hidden_variables = 2; %tracks value and uncertainty
        priors.muX0 = [zeros(n_basis,1); sigma_noise*ones(n_basis,1);zeros(1,1)];
        priors.SigmaX0 = zeros(hidden_variables*n_basis+1);
        h_name = @h_wtwsceptic_kalman;
        n_theta = 2;
        n_phi = 2; % one for return on policy, one for choice rule
        options.inF.u_aversion = u_aversion;
        options.inG.u_aversion = u_aversion;
    otherwise
        disp('The model you have entered does not match any of the default names, check spelling!');
        return
        
end

options.inF.hidden_state = hidden_variables;

%Map the necessary options from F to G
options.inG.hidden_state = options.inF.hidden_state;
options.inG.kalman = options.inF.kalman;

%Was the trial a win or not
wins = cellfun(@(x) strcmp(x,'win'), {data.trialResult});
%Was the trial an immediate quit or not
imm_quit = cellfun(@(x) strcmp(x,'optSmall'), {data.initialPos}) &  [data.latency]<0.1;

%Save for later
options.inF.wins = wins;
options.inF.imm_quit = imm_quit;


if multinomial
    %rtrnd = round([data.latency]'*(n_steps-iti)/range_RT)' + iti; %Add in ITI of 2 seconds
    rtrnd = round([data.latency]'*10);
    rtrnd(rtrnd==0)=1;
    %rtrnd(rtrnd>200)=200;
    max_rt = max(rtrnd);
    
    %add if statement for rr (+1 to n)
    dim = struct('n',hidden_variables*n_basis+1,'n_theta',n_theta+fit_propspread,'n_phi',n_phi,'p',n_steps);
    options.sources(1) = struct('out',1:n_steps,'type',2);
    
    
   %% compute multinomial response -- renamed 'y' here instead of 'rtbin'
    y = zeros(n_steps, length(trialsToFit));
    
    for i = 1:length(trialsToFit)
        if wins(i)~=1 %Only include non-wins
            y(rtrnd(i),i) =1;
        end
    end
    
    skip_mat = y; %This keeps track of all 'win' trials ie the trials the subjects waited we don't want to fit these.
    win_or_quit = {data.trialResult}; %Did the subject win or quit
    win_or_quit_idx = []; %For graphing purposes
    
    %For debug purposes save a pre-transformed rtrnd
    pre_trsnaformed_rtrnd = rtrnd;
    
    priors.a_alpha = Inf;   % infinite precision prior
    priors.b_alpha = 0;
    priors.a_sigma = 1;     % Jeffrey's prior
    priors.b_sigma = 1;     % Jeffrey's prior
    options.binomial = 1;
%     priors.muPhi = zeros(dim.n_phi,1); % exp tranform
%     priors.SigmaPhi = 1e1*eye(dim.n_phi);
    % Inputs
    time  = [data.initialTime];
    u = [[data.latency]*10; [data.payoff]; time(trialsToFit); rtrnd']; %data.latency * 10 to get from sec to ms
    u = [nan(size(u,1),1) u(:,1:end-1)];
    
    % Observation function
    switch model
        case 'kalman_logistic'
            g_name = @g_wtwsceptic_logistic;
        case 'win_stay_lose_switch'
            g_name = @g_WSLS;
        case 'stay'
            g_name = @g_WSLS;
        case 'null'
            g_name = @g_wtwsceptic_null;
        otherwise
            g_name = @g_wtwsceptic;
    end    
end
%% skip first trial and all trials with 'wins' ie they waited
options.skipf = zeros(1,n_t);
options.skipf(1) = 1;
%Do not use imm_quits in evolution function
options.skipf(imm_quit) = 1;
%options.skipf(wins) = 1;
%options.skipf(logical(sum(isnan(skip_mat)))) = 1;

%% priors
%dim.n_theta = 1;
if strcmp(model,'kalman_logistic')
    priors.muTheta = zeros(dim.n_theta,1);
    %priors.muTheta = 0.5;
    priors.SigmaTheta = 1e1*eye(dim.n_theta); % lower the learning rate variance -- it tends to be low in the posterior
    priors.muPhi = zeros(dim.n_phi,1); %Beta, gamma
    priors.SigmaPhi = 1e2*eye(dim.n_phi);
    
else
    priors.muTheta = zeros(dim.n_theta,1);
    %priors.muTheta = 0.5;
    priors.SigmaTheta = 1e1*eye(dim.n_theta); % lower the learning rate variance -- it tends to be low in the posterior
    priors.muPhi = [0, 0]; %Beta, gamma
    priors.SigmaPhi = 1e2*eye(dim.n_phi);
end

options.priors = priors;
options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)


%Remove the trials in which there was a large reward & the latency and between .1<20s
%3 conditions:
%all - all trials except wins
%qdf - all + no immediate quits
%wdf - qdf + no waits over 19.9

if strcmp(sample_to_use,'all')
    trial_index_to_remove = wins;
elseif strcmp(sample_to_use,'qdf')
    trial_index_to_remove = wins | imm_quit;
elseif strcmp(sample_to_use,'wdf')
    trial_index_to_remove = wins | imm_quit | [data.latency]>19.9;
else
    error('Which sample mto use all, qdf or wdf?')
end

options.isYout = zeros(size(y));
one_after_rt = 1;
if one_after_rt
    %Fill times after rt to nan
    for k = 1:length(wins)
        if wins(k)==1
            rt_idx = find(y(:,k)==1);
            try
                options.isYout(rt_idx+1:end,k) = 1;
            catch
            end
            options.isYout(rt_idx,k)=0;
        end
    end
end

if strcmp(sample_to_use,'qdf') || strcmp(sample_to_use,'wdf')
    options.isYout(:,imm_quit)=1;
end

figure(500)
subplot(1,2,1)
[time_bin,~]=find(y==1);
plot(time_bin, 'bo','LineWidth',3,'MarkerSize',10)
axis([0,length(time_bin),0,max(time_bin)])
title('Y Original')
subplot(1,2,2)
tmp_idx = sum(options.isYout)~=200;
surviving_data_points=time_bin;
surviving_data_points(~tmp_idx)=nan;
% [surviving_data_points,~]=find(y(:,tmp_idx)==1);
plot(surviving_data_points, 'ro','LineWidth',3,'MarkerSize',10)
axis([0,length(time_bin),0,max(time_bin)])
title(sprintf('Y of sub-sample %s',sample_to_use))

% options.priors.iQy = cell(1, length(u));
% for jj = 1 : length(trial_index_to_remove)
%     if trial_index_to_remove(jj)
%         options.priors.iQy(jj) = {repmat(0.00001,n_t,n_t)};
%     else
%         options.priors.iQy(jj) = cell(1,1);
%     end
% end

%%u(:,trial_index_to_remove)=[];

%IF we want to return return on policy for plotting
plot_rop = 0;
if plot_rop
    options.inF.plot_rop = 1;
    priors.muX0 = zeros(hidden_variables*n_basis+201,1);
    priors.SigmaX0 = zeros(hidden_variables*n_basis+201); %rop then rr
    dim.n=hidden_variables*n_basis+201;
    options.priors = priors;
    options.inG.priors = priors; %copy priors into inG for parameter transformation (e.g., Gaussian -> uniform)
end


 [posterior,out] = VBA_NLStateSpaceModel(y,u,h_name,g_name,dim,options);


% % if graphics==1
% %     figure(2); plot(pre_trsnaformed_rtrnd,'o')
% %     hold on
% %     pre_trsnaformed_rtrnd(~win_or_quit_idx)=nan;
% %     plot(pre_trsnaformed_rtrnd,'ro')
% %     title('Red is waits Blue is quits')
% %     
% %     %diagnose_wtw_sceptic()
% % end

if saveresults
%     cd(results_dir);
    %% save output figure
%     h = figure(1);
%     savefig(h,sprintf('%d_%s_multinomial%d_multisession%d_fixedParams%d_tau_rr%d',id,model,multinomial,multisession,fixed_params_across_runs,tau_rr))
     save([results_dir filesep sprintf('sample_%s_%d_%s_multinomial%d_multisession%d_fixedParams%d_uaversion%d_tau_rr%d_sceptic_vba_fit',sample_to_use, id, model, multinomial,multisession,fixed_params_across_runs, u_aversion,tau_rr)], 'posterior', 'out');
%      cd ..
%      cd ..
end