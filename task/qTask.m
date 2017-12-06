function [] = qTask()
% presents the qTask experiment

% skip screen test
Screen('Preference', 'SkipSyncTests', 1)
Screen('Preference', 'Verbosity', 1)
try
    
    %%% modifiable parameters
    % timing
    blockMins = 5; % block duration in minutes
    itiSecs = 2; % was previously 0.8
    timingDistribs = {'discrete1'}; 
        % timing distributions: will be presented in order
        % these are defined in the function DRAWSAMPLE
    % payoff contingencies
    payoff.win = 10; % points
    payoff.quit = 1; % points
    %pointConversion = 400; % points per dollar
	pointConversion = 100; % 100 cents per dollar
    
    %%% display preferences (set differently for windows)
    if IsWin
        % modify text display settings
        % set font to helvetica (default is courier new)
        Screen('Preference', 'DefaultFontName', 'Helvetica');
        % set style to normal (=0) (default is bold [=1])
        Screen('Preference', 'DefaultFontStyle', 0);
    end
    
    % set path
    path(path,'stimSubFx');
    
    % name datafile
    [dataFileName,dataHeader] = gatherSubInfo('qtask_upmc');
    
    % standard tasks
    bkgd = 80; % set shade of gray for screen background
    [wid,origin,dataHeader] = manageExpt('open',dataHeader,bkgd); % standard initial tasks
    dataHeader.distribs = timingDistribs; % log the sequence of distributions
    
    % set detailed timing parameters
    nBks = length(timingDistribs); % number of blocks to present
    blockSecs = blockMins*60; % convert minutes to seconds
    display.fbackDur = itiSecs/2; % first half of ITI is feedback
    display.iti = itiSecs/2; % second half of ITI is blank
    
    % initialize data logging structure
    trialData = struct([]);
    datarow = 0;
    
    % set screen locations for stimuli and buttons
    % rects has fields targ, button{1}, and button{2}
    rects = setRects(origin);

    % display initialization: prior to first block only
    display.totalWon = 0; % earnings initialized at zero
    if payoff.quit==1, pluralStr = ''; else pluralStr = 's'; end % buttons say either "point" or "points"
    %display.buttonText = {sprintf('Take %d cent%s',payoff.quit,pluralStr) sprintf('Wait for %d cents',payoff.win)};
	display.buttonText = {sprintf('Take $0.%02d',payoff.quit) sprintf('Wait for $0.%02d',payoff.win)};
    display.optNames = {'optSmall','optLarge'}; % to identify buttons for data logging
    
    % display initialization: prior to every block
    display.lightState = 'off'; % may take on 'waiting', 'small', 'large'
    display.lightText = '';
    display.choice = 'null'; % may take on 'optSmall', 'optLarge'
    display.timeLeft = blockSecs;
    
    % pack params into a struct
    params.wid = wid;
    params.bkgd = bkgd;
    params.blockSecs = blockSecs;
    params.payoff = payoff;
    
    % instructions and practice
    instrucBlock(params,rects,display,blockMins,pointConversion,nBks);
    
    % present individual blocks
    for bkIdx = 1:nBks
    
        % set block-specific parameters
        params.bkIdx = bkIdx;
        params.distrib = timingDistribs{bkIdx};
        
        % present a block of choice trials
        [display,trialData,datarow] = showTrials(params,rects,display,trialData,datarow);
        
        % save data
        save(dataFileName,'dataHeader','trialData');
        
        % intermediate instructions screen
        % (shown after all except the last block)
        if bkIdx<nBks
            msg = sprintf('Block %d complete.\n\nClick to begin block %d.',...
                bkIdx,bkIdx+1);
            showMsg(wid,bkgd,msg);
        end
        
    end % loop over blocks
    
    % show the provisional final screen
    %msg = sprintf('Complete!\n\nTotal money: %d ($%2.2f).\n\nPlease see the experimenter.',...
    %    display.totalWon,display.totalWon/pointConversion);
    msg = sprintf('Complete!\n\nTotal money: $%2.2f.\n\nPlease see the experimenter.',...
        display.totalWon/pointConversion);
	showMsg(wid,bkgd,msg,'q');
    
    % very last screen
    msg = 'Task complete.';
    showMsg(wid,bkgd,msg,'q');
    
    % close the expt
    manageExpt('close');
    
    % print some information
    fprintf('\n\nParticipant: %s\n',dataHeader.dfname);
    %fprintf('Winnings: %d points ($%2.2f)\n\n',display.totalWon,display.totalWon/pointConversion);
	fprintf('Winnings: $%2.2f\n\n',display.totalWon/pointConversion);
    
catch ME
    
    % close the expt
    disp(getReport(ME));
    manageExpt('close');
    
end % try/catch loop

    

    


