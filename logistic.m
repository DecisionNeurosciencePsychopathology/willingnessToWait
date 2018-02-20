function [wtw,shape,sigmoid]=logistic(a,choice_rand,u_func, ntimesteps, return_on_policy, tvec, discrim, u_threshold)

        %compared to other models that use a curve over which to choose (either by softmax or egreedy selection),
        %kalman_uv_logistic computes explore and exploit choices and chooses according to a logistic.
        u = sum(u_func)/length(u_func);

        if u == 0
            rt_explore = ceil(.5*ntimesteps);
        else
            rt_explore = find(u_func==max(u_func), 1); %return position of first max (and add gaussian noise?)
        end
    
        %Rasch model with tradeoff as difficulty (location) parameter
        sigmoid = 1/(1+exp(-discrim.*(u - u_threshold))); 
        
        %compute hard max of value function alone 
        if a == 1
                rt_exploit = ceil(.5*ntimesteps); %default to mid-point of time domain
        else
            rt_exploit = find(return_on_policy==max(return_on_policy), 1); %only take the first max if there are two identical peaks.
            if rt_exploit > max(tvec), rt_exploit = max(tvec); end
        end
        
        if choice_rand(a) < sigmoid
            %explore according to hardmax u
            wtw = rt_explore;
            shape = 150;
        else
            wtw = rt_exploit;
            shape = 100;
        end
        
       
        
end