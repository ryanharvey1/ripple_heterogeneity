% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'OR'),:);
% df.basepath = {'Z:\Data\FujisawaS\EE\EE0622fm'};
% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\tmaze_sessions.csv');
% opts = delimitedTextImportOptions("NumVariables", 1);
% opts.DataLines = [2, Inf];
% opts.Delimiter = ",";
% opts.VariableNames = "basepath";
% opts.VariableTypes = "string";
% opts.ExtraColumnsRule = "ignore";
% opts.EmptyLineRule = "read";
% opts = setvaropts(opts, "basepath", "WhitespaceRule", "preserve");
% opts = setvaropts(opts, "basepath", "EmptyFieldRule", "auto");
% df = readtable("Z:\home\ryanh\projects\ripple_heterogeneity\tmaze_sessions.csv", opts);

df.basepath = {pwd};
% visual of current edges (states) for linearized pos
%    ___ ___
%   | 2 | 1 |
%   |6  |0  |5
%   |__ | __|
%     4   3
% converted to
%    ___ ___
%   |   |   |
%   |   |0  |
%   |__ | __|
%   2       1
%
% This script will also bring down new state 2 to same level as state 1

for i = 1:length(df.basepath)
    disp(df.basepath{i})
    run(df.basepath{i})
end

function run(basepath)

basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

if ~isfield(behavior,'states')
    return
end
states = unique(behavior.states);
states = states(~isnan(states));
if length(states) == 3
    return
end


% figure to verify states
load(fullfile(basepath,[basename,'.session.mat']))
start = [];
stop = [];
for ep = 1:length(session.epochs)
    if contains(session.epochs{ep}.environment,'tmaze') ||...
            contains(session.epochs{ep}.environment,'Mwheel') ||...
            contains(session.epochs{ep}.environment,'Tmaze')
        
        start = session.epochs{ep}.startTime;
        stop = session.epochs{ep}.stopTime;
        
        [idx,~,~] = InIntervals(behavior.timestamps,...
            [start,stop]);
        
        if ~any(idx)
            continue
        else
            break
        end
    end
end
[epoch_idx,~,~] = InIntervals(behavior.timestamps,...
    [start,stop]);
epoch_idx = epoch_idx';


figure;
states = unique(behavior.states(epoch_idx));
states = states(~isnan(states));
subplot(1,2,1)
for s = states
    scatter(behavior.position.x(epoch_idx & behavior.states == s),...
        behavior.position.y(epoch_idx & behavior.states == s),'DisplayName',num2str(s))
    hold on
end
legend

subplot(1,2,2)
for s = states
    scatter(behavior.timestamps(epoch_idx & behavior.states == s),...
        behavior.position.linearized(epoch_idx & behavior.states == s),'DisplayName',num2str(s))
    hold on
end
legend

% drop states 2 and 4 down to the level of the other two arms (1&3) so
% right and left will be contingous
idx = epoch_idx & (behavior.states == 4 |...
    behavior.states == 6 |...
    behavior.states == 2);

behavior.position.linearized(idx) = behavior.position.linearized(idx) -...
    min(behavior.position.linearized(idx)) +...
    max(behavior.position.linearized(epoch_idx & behavior.states == 0));

% relabel states
states = behavior.states;
idx = epoch_idx & (behavior.states == 4 |...
    behavior.states == 6 |...
    behavior.states == 2);
states(idx) = 2;

idx = epoch_idx & (behavior.states == 3 |...
    behavior.states == 5 |...
    behavior.states == 1);
states(idx) = 1;
behavior.states(epoch_idx) = states(epoch_idx);

% update statenames
behavior.stateNames = {'middle','right','left'};

colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
figure;
subplot(1,2,1)
states = unique(behavior.states(epoch_idx));
states = states(~isnan(states));
for s = states
    scatter(behavior.position.x(epoch_idx & behavior.states == s),...
        behavior.position.y(epoch_idx & behavior.states == s),[],colors{s+1})
    hold on
end
title('whole first epoch')

subplot(1,2,2)
states = unique(behavior.states(epoch_idx));
states = states(~isnan(states));
for s = states
    scatter(behavior.timestamps(epoch_idx & behavior.states == s),...
        behavior.position.linearized(epoch_idx & behavior.states == s),[],colors{s+1})
    hold on
end
% xlim([start,start+60*5])
% title('first 5 minutes linear')


prompt = "Does this look good? y/n: ";
txt = input(prompt,"s");
if contains(txt,'y')
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end

end