%test softmax policy
% no need for this 

function wtw=softmax_policy(tvec, p_choice,softmax_stream)

wtw = randsample(softmax_stream, tvec, 1, true, p_choice);

end