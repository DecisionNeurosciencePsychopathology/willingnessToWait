
function wtw_analysis

%% default variables
    n_basis = 16;
    multinomial = 1;
    multisession = 0;
    fixed_params_across_runs = 0;
    fit_propspread = 0;
    n_steps = 200;
    u_aversion = 0;
    saveresults = 0; %change later
    graphics = 0;

%models

model = {'fixed', 'null','kalman_uv_sum', 'kalman_logistic'};

%% start 
load('wtw_data.mat')

%collect all ids
for j = 1:length(wtw_struct)
    if isempty(wtw_struct(j).id) == 0
        ids(l) = wtw_struct(j).id;
        l= l+1;
    else  l = j;
    end
end

for h = 0:1
    tau_rr = h;
    
    disp('Tau_RR = ')
    disp(h)
    
    for k = 1:length(model)
        model = char(model(k));
        disp('Model: ')
        disp(model)

        for i = 1:length(ids)
            for p = 1:length(wtw_struct)
               if wtw_struct(p).id == ids(i)
                    data = wtw_struct(p).trialData;
               end
            end

            
            id = ids(i);

            disp('ID: ')
            disp(id)

            wtw_sceptic_vba(data,id,model,n_basis, multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion, tau_rr, saveresults, graphics);

        end
    end
end

end