%Pass the wtw data to the vba function
clear;
clc;

%Load in main data matrix
load('C:\kod\Neuropsych_preproc\matlab\analysis\willingness to wait\data\wtw_data.mat');

%clean struct (ie remove empties)
empty_elms = cellfun(@(s) isempty(s), {wtw_struct.id});
wtw_struct(empty_elms) = [];

%Set up variables
modelnames = {'fixed_opportunity_cost' 'fixed' 'fixed_uv' 'fixed_decay' 'kalman_softmax' 'kalman_processnoise' 'kalman_uv_sum' 'kalman_sigmavolatility' 'kalman_logistic'};
%model = 'fixed_decay';
nbasis = 16;
multinomial = 1;
multisession = 0;
fixed_params_across_runs = 1;
fit_propspread = 0;
iti=20;
n_steps = 200+iti;
u_aversion = 1; % allow for uncertainty aversion in UV_sum
saveresults = 0; %don't save to prevent script from freezing on Thorndike
graphics = 1;

%Grab the ids
id = 46069;
ids = [wtw_struct.id];







for m=1:length(modelnames)
    model = char(modelnames(m));
    for i = 1:length(ids)
        %Grab subj id
        id = ids(i);
        
        %Find where to pull the subj data from the main struct
        subj_idx=find([wtw_struct.id] == id);
        data = wtw_struct(subj_idx).trialData;
        
        %Run the main vba function
        [posterior,out] = wtw_sceptic_vba(data,id,model,nbasis, multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion, saveresults, graphics);
        
        %Update user
        fprintf('Fitting %s subject %d \r',model,id)
        
        %Save the log evidence
        L(m,i) = out.F;
    end
end