%For diagnosis purposes only
gamma = exp(post.muPhi(2));
fx=out.suffStat.muX;
gx = out.suffStat.gx;
iti = 20; 
opp_cost_scaling = 1/iti;
tau = 1;
for i=1:out.dim.n_t
    
gaussmat = out.options.inF.gaussmat;
ntimesteps = out.options.inF.ntimesteps;
nbasis = out.options.inF.nbasis;

v=fx(1:nbasis,i)*ones(1,ntimesteps) .* gaussmat; %use vector outer product to replicate weight vector



% add cumsum
w = fx(1:nbasis,i);
u = zeros(1,ntimesteps);
if out.options.inG.kalman.kalman_uv_sum || out.options.inG.kalman.fixed_uv
   tau = 1./(1+exp(-post.muPhi(3)-10)); %Uncertainty mixing: 0..1
   sigma=fx(nbasis+1:nbasis*2,i);
   
   %perform tau mixing here
   u = sigma*ones(1,ntimesteps) .* gaussmat;
end
u_func = sum(u);
v_func = sum(v); %subjective value by timestep as a sum of all basis functions

% add cumsum
v_func = cumsum(v_func);

%Scale v_func down to the max of the basis
v_func = v_func./max(v_func)*max(w);

%If we are dividing by 0
if isnan(v_func)
    v_func = zeros(1,ntimesteps);
end



% add opportunity cost = rr*trial length (max 200)
opp_cost(i,:) = fx(end,i).*(1:ntimesteps);

cumulative_reward_fx(i,:) = v_func;

return_on_policy(:,i) = cumulative_reward_fx(i,:)-(gamma*opp_cost_scaling*opp_cost(i,:));

%Update with uncertainty if appliciable
return_on_policy(:,i) = tau .* return_on_policy(:,i) + (1-tau).*u_func;

end

%Define the wins
rtrnd = round([data.latency]'*10);
win_data = nan(out.dim.n_t,1);
win_data(out.options.inF.wins) = rtrnd(out.options.inF.wins);

%Define the plot data for the y input
plot_y_data = nan(out.dim.n_t,1);

%Find where the y was actually being used 
[y_data,idx]=find(out.y==1);
plot_y_data(idx) = y_data;

%Remove the immediate quits
plot_y_data(logical(out.options.skipf))=nan;

figure(999)
%Value - MUX
subplot(6,1,1)
%surf(out.suffStat.muX(1:nbasis,:))
imagesc(out.suffStat.muX(1:nbasis,:))
title('nbasis functions')
hold on
plot(1:out.dim.n_t,plot_y_data*nbasis/200,'r*')
plot(1:out.dim.n_t,win_data*nbasis/200,'b*')
colorbar

%RR - reward rate
subplot(6,1,2)
plot(out.suffStat.muX(end,:),'Linewidth',3)
colorbar
title('Reward Rate')

%Cumulative reward funciton
subplot(6,1,3)
imagesc(cumulative_reward_fx')
colorbar
title('cumulative reward fx ')

%Opportunity cost
subplot(6,1,4)
%[data,idx]=find(out.options.skipf~=1);
imagesc(opp_cost')
colorbar
title(sprintf('Opportunity cost with gamma: %.2f', gamma))

%Rop - return on policy
subplot(6,1,5)
%[data,idx]=find(out.options.skipf~=1);
imagesc(return_on_policy)
title('Return on policy')
colormap bone
colorbar
hold on
plot(1:out.dim.n_t,plot_y_data,'r*')

subplot(6,1,6)
%[data,idx]=find(out.options.skipf~=1);
imagesc(gx)
title('Model response')
colormap bone
colorbar
hold on
plot(1:out.dim.n_t,plot_y_data,'r*')









