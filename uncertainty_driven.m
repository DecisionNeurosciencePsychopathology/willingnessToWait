function uncertainty_driven(i,distrib,nbasis,ntimesteps)

%assign reward
switch distrib
    case 'gp'
        
        
    case 'unif'
        

    case 'beta'
        
 
        
    case 'discrete1'        
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
    
    
%     %include Q???
%     Q_ij(i+1,:) = p.omega.*abs(delta_ij(i,:)); %use abs of PE so that any large surprise enhances effective gain.
%     %Compute the Kalman gains for the current trial (potentially adding process noise)
%     k_ij(i,:) = (sigma_ij(i,:) + Q_ij(i,:))./(sigma_ij(i,:) + Q_ij(i,:) + sigma_noise); 

    %Kalman gain
    k_ij(i,:)=0.5;
    
    %Uncertainty is a function of Kalman uncertainties.
    u_jt=sigma_ij(i+1,:)'*ones(1,ntimesteps) .* gaussmat;        
    u_func = sum(u_jt); %vector of uncertainties by timestep
    
    %compute a "tilt" vector that linearly scales with timepoint according to PE size * parameter
    %use different parameters, kappa and lambda, to handle PE+ and PE-, respectively
    if max(delta_ij(i,:)) > 0
        discount = p.kappa*(max(delta_ij(i,:)))*tvec; %% also add valence-dependent parameters: kappa for PE+, lambda for PE-
    else
        discount = p.lambda*(min(delta_ij(i,:)))*tvec;
    end
    
    
    v_disc = v_func + discount;
    %v_it_undisc(i+1,:) = v_func;
    uv_func=p.tau*v_disc + (1-p.tau)*u_func;
    v_final = uv_func;

end
    
    
    
        