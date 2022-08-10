function [SCORE, trials] = wMazeScoreTrials(basePath, sr, saveMat)
% [SCORE, trials] = wMazeScoreTrials(basePath, sr, saveMat)
% calculate performance W-maze

cd(basePath)
if nargin == 1
    sr = 20000;
end
if nargin < 3
   saveMat = 1;
end

%%
data = LoadBinary('digitalin.dat','nChannels',1);
% digitalin channel code:
%   1- optitrack
%   2- right arm 
%   3- center arm 
%   4- left arm 
%   5- reward deliverly

%Separating the data to the channels they belong.
DigIn = zeros(length(data),5);

aux = data-32;
auxIdx = find(aux>=0);
DigIn(auxIdx,5) = 1;
data(auxIdx) = data(auxIdx)-32;

aux = data-16;
auxIdx = find(aux>=0);
DigIn(auxIdx,4) = 1;
data(auxIdx) = data(auxIdx)-16;

aux = data-8;
auxIdx = find(aux>=0);
DigIn(auxIdx,3) = 1;
data(auxIdx) = data(auxIdx)-8;

aux = data-4;
auxIdx = find(aux>=0);
DigIn(auxIdx,2) = 1;
data(auxIdx) = data(auxIdx)-4;

% auxIdx = find(data>0);
% DigIn(auxIdx,1) = 1;  

%% Get sensor crossing times
auxDI = diff(DigIn);
auxLowSiz = find(auxDI(:,3) == 1);
auxSupSiz = find(auxDI(:,3) == -1);
auxSiz = find((auxSupSiz - auxLowSiz)>(sr/2)-1);
auxTS = find(auxDI(:,3) == 1);
Marm = auxTS(auxSiz) + 1; % middle arm sensor crossing

auxLowSiz = find(auxDI(:,2) == 1);
auxSupSiz = find(auxDI(:,2) == -1);
auxSiz = find((auxSupSiz - auxLowSiz)>(sr/2)-1);
auxTS = find(auxDI(:,2) == 1);
Rarm = auxTS(auxSiz) + 1; % right arm sensor crossing

auxLowSiz = find(auxDI(:,4) == 1);
auxSupSiz = find(auxDI(:,4) == -1);
auxSiz = find((auxSupSiz - auxLowSiz)>(sr/2)-1);
auxTS = find(auxDI(:,4) == 1);
Larm = auxTS(auxSiz) + 1; % left arm sensor crossing

auxLowSiz = find(sum(auxDI(:,2:4),2) == 1);
auxSupSiz = find(sum(auxDI(:,2:4),2) == -1);
if( length(auxSupSiz) ~= length(auxLowSiz) )
    auxSupSiz = auxSupSiz( (length(auxSupSiz) - length(auxLowSiz)) + 1 : end );
end
auxSiz = find((auxSupSiz - auxLowSiz)>(sr/2)-1);
auxTS = find(sum(auxDI(:,2:4),2) == 1);
Trials = auxTS(auxSiz) + 1; % all sensor crossing

auxLowSiz = find(auxDI(:,5) == 1);
auxSupSiz = find(auxDI(:,5) == -1);
auxSiz = find((auxSupSiz - auxLowSiz)>(sr/2)-1);
auxTS = find(auxDI(:,5) == 1);
Reward = auxTS(auxSiz) + 1; % all sensor crossing with reward delivery 

%% Scoring Behavior
    % Asuming rat always starts in right arm    
iBerror = 0; oBerror = 0;
iBtrial = 0; oBtrial = 0;
% arm code: 0=left; 1=middle; 2=right;

for idx = 1:length(Trials)
    if sum(Trials(idx) >= Larm-3 & Trials(idx) <= Larm+3) % L arm sensor
        arm(idx) = 0;
        timestamps(idx) = Trials(idx);
        if sum(Trials(idx) >= Reward-5 & Trials(idx) <= Reward+5) % if there was a reward
            trialsS(idx) = 1; % correct trial
        else
            trialsS(idx) = 0; % error trial
        end
        %checking wrong trials
        if idx>1
            if (arm(idx - 1) ~= 1)
                iBerror = iBerror + 1;
                iBtrial = iBtrial + 1;
            elseif (idx>2) && (arm(idx-2) == 0 && arm(idx - 1) == 1)
                oBerror = oBerror + 1;
                oBtrial = oBtrial + 1;
            elseif (idx>2) &&(arm(idx-2) == 2 && arm(idx - 1) == 1)
                oBtrial = oBtrial + 1;
            end
        end
    elseif sum(Trials(idx) >= Marm-3 & Trials(idx) <= Marm+3) % M arm sensor
        arm(idx) = 1;
        timestamps(idx) = Trials(idx);
        if sum(Trials(idx) >= Reward-5 & Trials(idx) <= Reward+5)
            trialsS(idx) = 1;
        else
            trialsS(idx) = 0;
        end
        iBtrial = iBtrial + 1;
    elseif sum(Trials(idx) >= Rarm-5 & Trials(idx) <= Rarm+5)  % R arm sensor
        arm(idx) = 2;
        timestamps(idx) = Trials(idx);
        if sum(Trials(idx) >= Reward-5 & Trials(idx) <= Reward+5)
            trialsS(idx) = 1;
        else
            trialsS(idx) = 0;
        end
        if idx>2
            if (arm(idx - 1) ~= 1)
                iBerror = iBerror + 1;
                iBtrial = iBtrial + 1;
            elseif (arm(idx-2) == 2 && arm(idx-1) == 1)
                oBerror = oBerror + 1;
                oBtrial = oBtrial + 1;
            elseif (arm(idx-2) == 0 && arm(idx - 1) == 1)
                oBtrial = oBtrial + 1;
            end
        end
    end
end
if arm(1) ~= 2
    warning('rat did not start in correct arm');
end

% output behav score
SCORE.trials = trialsS;
SCORE.arm = arm;
SCORE.timestamps = timestamps/sr;
SCORE.performance = (-oBerror-iBerror+iBtrial+oBtrial)/(oBtrial+iBtrial);
SCORE.inbound = (iBtrial - iBerror) / iBtrial;
SCORE.outbound = (oBtrial - oBerror) / oBtrial;
SCORE.iBtrial = iBtrial;
SCORE.iBerror = iBerror;
SCORE.oBtrial = oBtrial;
SCORE.oBerror = oBerror;

if saveMat
    save([basePath '.performance2.mat'],'SCORE');
end

%% get trial structure

for i = 2:length(Trials)
    % trials start / stop
        trials.int(i-1,1) = Trials(i-1)/sr;
        trials.int(i-1,2) = Trials(i)/sr;
    
    % trial type
    if arm(i) == 1 % inbound 
       trials.inbound(i-1,1) = 1;
       trials.outbound(i-1,1) = 0;
    else
       trials.outbound(i-1,1) = 1;
    end
    
    % trial direction
    if arm(i) == 1
       if arm(i-1) == 0 % inbound coming from left
          trials.left(i-1,1) = 1; 
          trials.right(i-1,1) = 0; 
       elseif arm(i-1) == 2
          trials.left(i-1,1) = 0;
          trials.right(i-1,1) = 1; 
       end
    else
       if arm(i) == 0 % inbound coming from left
          trials.left(i-1,1) = 1; 
          trials.right(i-1,1) = 0; 
       elseif arm(i) == 2
          trials.left(i-1,1) = 0;
          trials.right(i-1,1) = 1; 
       end        
    end
    
    % correct / error trial
    trials.correct = trialsS(2:end)';
    
end

if saveMat
   save([basePath '.WmazeTrials.mat'],'trials');
end

