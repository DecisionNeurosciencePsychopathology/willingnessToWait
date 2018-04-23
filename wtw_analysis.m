
function L = wtw_analysis

%% default variables
    n_basis = 16;
    multinomial = 1;
    multisession = 0;
    fixed_params_across_runs = 0;
    fit_propspread = 0;
    n_steps = 200;
    u_aversion = 0;
    saveresults = 0; %change later
    graphics = 1;

%models

%model_list = {'fixed', 'fixed_decay', 'null','kalman_uv_sum', 'kalman_logistic'};
samples_to_use = {'all', 'qdf', 'wdf'};
model_list = {'fixed'};

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

%initialize
L = zeros((length(model_list)*2),(length(ids)));
L_id = zeros(1,length(ids));

%Override ids so it is only a good subject
%ids=214710;
ids = 216174;

for h = 0%:1 %Tau_rr or no
    tau_rr = h;
    
    disp('Tau_RR = ')
    disp(h)
    
    for k = 1:length(model_list)
        model = char(model_list(k));
        disp('Model: ')
        disp(model)
        
        if tau_rr
            m = 4+k; %row for L
        else 
            m = k;
        end

        for i = 1:length(ids)
            for p = 1:length(wtw_struct)
               if wtw_struct(p).id == ids(i)
                    data = wtw_struct(p).trialData;
               end
            end

            
            id = ids(i);
            
            disp('ID: ')
            disp(id)
% %             try                
                [post, out] = wtw_sceptic_vba(data,id,samples_to_use{2},model,n_basis, multinomial,multisession,fixed_params_across_runs,fit_propspread,n_steps,u_aversion, tau_rr, saveresults, graphics);
                diagnose_mux
                L(m, i) = out.F;
                L_id(1,i) = id;
% %             catch
% %                 L(m, i) = nan;
% %                 L_id(1,i) = nan;
% %                 fprintf('Subject %d did not complete the State space model\n\n',id)
% %                 uncompleted_ids{m,i} = id;
% %             end
        end
    end
end


save(sprintf('all_outF_%s', date),'L');

end