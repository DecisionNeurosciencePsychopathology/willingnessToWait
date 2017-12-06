function [display,trialData,datarow] = showTrials(params,rects,display,trialData,datarow)
% show a block of choice trials

% unpack parameters
wid = params.wid;
bkgd = params.bkgd;
blockSecs = params.blockSecs;
distrib = params.distrib;
payoff = params.payoff;
bkIdx = params.bkIdx;

% display initialization: prior to every block
display.lightState = 'off'; % may take on 'waiting', 'small', 'large'
display.lightText = '';
display.choice = 'null'; % may take on 'optSmall', 'optLarge'
display.timeLeft = blockSecs;
quantile = []; % this variable will control draws from the timing distribution

% show the initial screen
refreshDisplay(wid,bkgd,rects,display);
changeDisp = false; % register whether anything has changed on each loop

% show trials
% each trial consists of 3 phases: feedback, iti, waiting
SetMouse(rects.optLarge.center(1),rects.optLarge.center(2),wid); % initially position in the 'waiting' box
blockOnset = GetSecs;
trialPhase = 'iti';
phaseOnset = GetSecs;
while (GetSecs-blockOnset) < blockSecs % proceed continuously until time is up
%     reset=1;
%     if reset
%         changeDisp = false;
%         display.choice = display.optNames{2};
%         refreshDisplay(wid,bkgd,rects,display);
%     end
    
    
    % update the option selected
    [x,y] = GetMouse(wid);
    
    for b = 1:2
        if ~strcmp(display.choice,display.optNames{b}) && IsInRect(x,y,rects.(display.optNames{b}).rect) % if selection has changed
            display.choice = display.optNames{b};
            changeDisp = true;
        end
    end

    % determine how far we are into the current trial phase
    eventLatency = GetSecs-phaseOnset;

    % if the ITI is ending, begin a trial
    if strcmp(trialPhase,'iti') && eventLatency>display.iti
        trialPhase = 'waiting';
        phaseOnset = GetSecs;
        [waitDuration,quantile] = drawSample(distrib,quantile); % waiting time on this trial
        datarow = datarow+1; % number for the new trial
        display.lightState = 'waiting';
        %Reset the choice to wait after each trial
        display.choice = display.optNames{2}; 
        SetMouse(rects.optLarge.center(1),rects.optLarge.center(2),wid); 
        changeDisp = true;

        % record data pertaining to the new trial
        trialData(datarow).bkIdx = bkIdx; % block number
        trialData(datarow).initialTime = GetSecs-blockOnset; % trial onset time
        trialData(datarow).initialPos = display.choice; % initial cursor position

    % if the small/soon option is being selected, end the trial
    elseif strcmp(trialPhase,'waiting') && strcmp(display.choice,'optSmall')
        trialPhase = 'feedback';
        phaseOnset = GetSecs;
        display.lightState = 'small';
        display.lightText = sprintf('%d',payoff.quit);
        display.totalWon = display.totalWon+payoff.quit;
        changeDisp = true;

        % log data
        trialData(datarow).designatedWait = waitDuration;
        trialData(datarow).trialResult = 'quit';
        trialData(datarow).latency = eventLatency;
        trialData(datarow).payoff = payoff.quit;
        trialData(datarow).totalEarned = display.totalWon;
        trialData(datarow).timeLeft = display.timeLeft; % the displayed time left
        trialData(datarow).outcomeTime = GetSecs-blockOnset; % exact time elapsed

    % if a win has arrived, end the trial
    elseif strcmp(trialPhase,'waiting') && eventLatency>waitDuration
        trialPhase = 'feedback';
        phaseOnset = GetSecs;
        display.lightState = 'large';
        display.lightText = sprintf('%d',payoff.win);
        display.totalWon = display.totalWon+payoff.win;
        changeDisp = 1;

        % log data
        trialData(datarow).designatedWait = waitDuration;
        trialData(datarow).trialResult = 'win';
        trialData(datarow).latency = eventLatency;
        trialData(datarow).payoff = payoff.win;
        trialData(datarow).totalEarned = display.totalWon;
        trialData(datarow).timeLeft = display.timeLeft; % the displayed time left
        trialData(datarow).outcomeTime = GetSecs-blockOnset; % exact time elapsed

    % if the feedback interval is over, enter the ITI
    elseif strcmp(trialPhase,'feedback') && eventLatency>display.fbackDur
        trialPhase = 'iti';
        phaseOnset = GetSecs;
        display.lightState = 'off';
        display.lightText = '';
        changeDisp = true;
    end

    % update the time remaining
    timeLeft = ceil((blockOnset+blockSecs)-GetSecs);
    if display.timeLeft~=timeLeft
        display.timeLeft = timeLeft;
        changeDisp = true;
    end

    % display any updates to the display on this cycle
    if changeDisp
        changeDisp = false;
        refreshDisplay(wid,bkgd,rects,display);
    end

    WaitSecs(.001);

end


