function [dataFileName,dataHeader] = gatherSubInfo(tag)
% GATHERSUBINFO requests the subject's ID code.
% Creates a unique filename for saving the data.  
% Returns some relevant info in dataHeader.
% Input:  
%  TAG is the experiment name to be used in the data file. (a text string)

id = [];
while isempty(id)
    id = input('Subject ID:  ','s');
end

session = 1;
nameVetted = false;
while ~nameVetted
    dataFileName = fullfile('data',sprintf('%s_%s_%d',tag,id,session));
    if exist(sprintf('%s.mat',dataFileName),'file')==2
        session = session+1;
    elseif exist(sprintf('%s.txt',dataFileName),'file')==2
        session = session+1;
    else
        nameVetted = true;
    end
end

dataHeader.id = id;
dataHeader.dfname = dataFileName;

input(['Data will be saved in ',dataFileName,' (ENTER to continue)']);