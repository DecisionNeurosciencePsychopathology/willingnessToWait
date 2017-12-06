function [] = refreshDisplay(wid,bkgd,rects,display)
% places stimuli on the screen
% trial-specific parameters are in the struct display

% clear the screen
Screen('FillRect',wid,bkgd);

% reward indicator
lightColor.off = bkgd;
lightColor.waiting = [200 200 0];
lightColor.small = [250 80 80];
lightColor.large = [80 80 250];
Screen('FillOval',wid,lightColor.(display.lightState),rects.light.rect);
if ~isempty(display.lightText)
    Screen('TextSize',wid,30);
    bounds = Screen('TextBounds',wid,display.lightText);
    DrawFormattedText(wid,display.lightText,'center',rects.light.center(2)-RectHeight(bounds)*.4);
end

% buttons
buttonNames = {'optSmall','optLarge'};
for b = 1:2
    if strcmp(display.choice,buttonNames{b}) % if this button is now selected
        Screen('FillRect',wid,255,rects.(buttonNames{b}).border); % highlight the border
    end
    Screen('FillRect',wid,130,rects.(buttonNames{b}).rect); % show the button
    Screen('TextSize',wid,20);
    bounds = Screen('TextBounds',wid,display.buttonText{b});
    x = rects.(buttonNames{b}).center(1)-RectWidth(bounds)/2;
    y = rects.(buttonNames{b}).center(2)-RectHeight(bounds)/2;
    DrawFormattedText(wid,display.buttonText{b},x,y); % mark it
    Screen('TextSize',wid,18);
end

% earnings
%earningsStr = sprintf('Points earned:  %d',display.totalWon);
earningsStr = sprintf('Money earned: $%2.2f',display.totalWon/100);
DrawFormattedText(wid,earningsStr,rects.earnings(1),rects.earnings(2));

% time
timeStr = sprintf('Time left:  %02d:%02d',floor(display.timeLeft/60),floor(mod(display.timeLeft,60)));
DrawFormattedText(wid,timeStr,rects.time(1),rects.time(2));

% show the screen
Screen('Flip',wid);


