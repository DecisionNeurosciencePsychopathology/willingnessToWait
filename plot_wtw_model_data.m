function plot_wtw_model_data(ret,wtw)

%plot RTs

%output total reward

rew=ret.rew_i;

tot_rew=sum(rew);

i=1;
while wtw(i)~=0
    wtw_(i)=wtw(i);
    i=i+1;
end

plot(wtw_)

end