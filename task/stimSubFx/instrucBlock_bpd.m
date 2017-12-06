function [] = instrucBlock_bpd(params,rects,display,blockMins,pointConversion,nBks)
% presents instructions screens for qTask

% unpack parameters
wid = params.wid;
bkgd = params.bkgd;
payoff = params.payoff;

% first instruction
msg = {};
msg{1} = 'You will see a yellow light on the screen. The light represents a chance for you to earn money.';
msg{2} = 'You can use the touchscreen to choose between two options.';
msg{3} = sprintf('One option is, "Wait for $0.%02d." If you choose this, the yellow light will eventually go out and pay you %d cents.',payoff.win,payoff.win);
msg{4} = sprintf('The other option is, "Take $0.01." If you choose this, you can put the light out immediately and get paid 1 cent.');
msg{5} = 'Click to see some practice.';

txt = sprintf('%s\n\n',msg{:});
showMsg(wid,bkgd,txt);

% present 2 practice trials
SetMouse(rects.optLarge.center(1),rects.optLarge.center(2),wid); % initially position in the 'waiting' box
trialPhase = 'iti';
phaseOnset = GetSecs;
exptOnset = GetSecs;
sessSecs = 60;
trialNum = 1;
durations = [10, 10, 10, 10]; % determines the number of practice trials
itiLengths = display.iti*ones(size(durations)); % for the 2 practice trials (first is long)
%while trialNum<=length(durations)
while (GetSecs - exptOnset) < sessSecs


    % update the option selected
    [x,y] = GetMouse(wid);
    for b = 1:2
        if ~strcmp(display.choice,display.optNames{b}) && IsInRect(x,y,rects.(display.optNames{b}).rect) % if selection has changed
            display.choice = display.optNames{b};
            changeDisp = 1;
        end
    end

    % determine how far we are into the current trial phase
    eventLatency = GetSecs-phaseOnset;

    % if the ITI is ending, begin a trial
    if strcmp(trialPhase,'iti') && eventLatency>itiLengths(1) % SPECIAL duration of instruction ITI
        trialPhase = 'waiting';
        phaseOnset = GetSecs;
        waitDuration = durations(1); % SPECIAL for the example trial
        display.lightState = 'waiting';
        %Reset the choice to wait after each trial
        display.choice = display.optNames{2}; 
        SetMouse(rects.optLarge.center(1),rects.optLarge.center(2),wid); 
        changeDisp = 1;

    % if the small/soon option is being selected, end the trial
    elseif strcmp(trialPhase,'waiting') && strcmp(display.choice,'optSmall')
        trialPhase = 'feedback';
        phaseOnset = GetSecs;
        display.lightState = 'small';
        display.lightText = sprintf('%d',payoff.quit);
        display.totalWon = display.totalWon+payoff.quit;
        changeDisp = 1;

    % if a win has arrived, end the trial
    elseif strcmp(trialPhase,'waiting') && eventLatency>waitDuration
        trialPhase = 'feedback';
        phaseOnset = GetSecs;
        display.lightState = 'large';
        display.lightText = sprintf('%d',payoff.win);
        display.totalWon = display.totalWon+payoff.win;
        changeDisp = 1;

    % if the feedback interval is over
    elseif strcmp(trialPhase,'feedback') && eventLatency>display.fbackDur
        trialNum = trialNum + 1; % special for tracking practice trials
        trialPhase = 'iti';
        phaseOnset = GetSecs;
        display.lightState = 'off';
        display.lightText = '';
        changeDisp = 1;
    end

    % update the time remaining
    timeLeft = ceil((exptOnset+sessSecs)-GetSecs);
    if display.timeLeft~=timeLeft
        display.timeLeft = timeLeft;
        changeDisp = 1;
    end

    % display any updates to the display on this cycle
    if changeDisp==1
        changeDisp = 0;
        refreshDisplay(wid,bkgd,rects,display);
    end

    WaitSecs(.001);

end

% 2nd instruction
msg = {};
msg{1} = sprintf('The yellow light always wins %d cents eventually, but sometimes it may take a long time. You can change your choice anytime.',payoff.win);
if nBks>1 % instrucs differ for 1 block or multiple blocks.
    bkDescrip = sprintf('You will play %d blocks of the task, each lasting %d minutes.',nBks,blockMins);
else
    bkDescrip = sprintf('You will have %d minutes to play.',blockMins);
end
msg{2} = sprintf('%s You should try to earn as much as you can. As an incentive, you will be paid whatever you win.',...
    bkDescrip); %,pointConversion);
msg{3} = 'Click to start.';

txt = sprintf('%s\n\n',msg{:});
showMsg(wid,bkgd,txt);




