%test softmax policy

function wtw=softmax_policy(tvec, p_choice)

softmax_seed=21; %hardcoded

softmax_stream = RandStream('mt19937ar','Seed',softmax_seed);

wtw = randsample(softmax_stream, tvec, 1, true, p_choice);

end