function [] = quickLook()
% single-subject analysis

% identify the data file
[fname,pathname] = uigetfile('data/*');
d = load(fullfile(pathname,fname));

% format data
[subInfo, trials] = formatData(d);

% individual timecourse plot
for b = 1:length(trials)
    subplot(1,length(trials),b);
    ssPlot(b,subInfo,trials(b));
end

end



% SUBFUNCTIONS

%%%%%
% subfunction to format one subject's data
function [subInfo, trials] = formatData(d)

% assess the number of blocks
bkIdx = [d.trialData.bkIdx]';
nBks = max(bkIdx);

trials = struct([]);
for b = 1:nBks

    idx = bkIdx==b;
    
    % add data fields for trials
    trials(b).trialNums = (1:sum(idx))';
    trials(b).designatedWait = [d.trialData(idx).designatedWait]';
    trials(b).outcomeWin = strcmp({d.trialData(idx).trialResult}','win');
    trials(b).outcomeQuit = strcmp({d.trialData(idx).trialResult}','quit');
    trials(b).latency = [d.trialData(idx).latency]';
    trials(b).startTime = [d.trialData(idx).initialTime]';
    trials(b).outcomeTime = [d.trialData(idx).outcomeTime]';
    trials(b).totalEarned = [d.trialData(idx).totalEarned]';
    trials(b).initialPosLarge = strcmp({d.trialData(idx).initialPos}','optLarge');
    trials(b).initialPosSmall = strcmp({d.trialData(idx).initialPos}','optSmall');
    
end

% display some info
subInfo.distribs = d.dataHeader.distribs;
subInfo.id = d.dataHeader.id;
subInfo.points = trials(b).totalEarned(end);
subInfo.money = subInfo.points/400;
fprintf('id: %s\npoints: %d\nmoney: %2.2f\n',...
    subInfo.id,subInfo.points,subInfo.money);

end % end of subfunction



%%%%%
% subfunction to create individual timecourse plot
function [] = ssPlot(b,sub,d)

titleText = sprintf('%s, cond = %s',sub.id,sub.distribs{b});

isWin = d.outcomeWin;
isQuit = d.outcomeQuit;

hold on;
plot(d.trialNums(isWin),d.designatedWait(isWin),'bo-','LineWidth',1);
plot(d.trialNums(isQuit),d.latency(isQuit),'ro-','LineWidth',1);
ceilingWait = d.designatedWait; 
ceilingWait(ceilingWait>20) = 20;
plot(d.trialNums(isQuit),ceilingWait(isQuit),'ko','LineWidth',1);
hold off;
h=get(gca,'Children');
for i=1:length(h), set(h(i),'MarkerSize',3); end
set(gca,'YLim',[0 20],'FontSize',8);
ylabel('time (s)');
xlabel('trial number');
%legend('Success','Quit');
title(titleText);

end % end of subfunction




