% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'OR'),:);
df.basepath = {pwd}
% visual of current edges (states) for linearized pos
%
%   |   |   |
%   |4  |0  |3
%   |__ | __|
%     2   1
% converted to 
%   |   |   |
%   |   |0  |
%   |__ | __|
%   2       1
%
for i = 1:length(df.basepath)
    disp(df.basepath{i})
   run(df.basepath{i}) 
end

function run(basepath)

basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.animal.behavior.mat']))

states = unique(behavior.states);
states = states(~isnan(states));
if length(states) == 3
   return 
end
% drop states 2 and 4 down to the level of the other two arms (1&3) so 
% right and left will be contingous
idx = behavior.states == 2 | behavior.states == 4;

behavior.position.linearized(idx) = behavior.position.linearized(idx) -...
    min(behavior.position.linearized(idx)) +...
    max(behavior.position.linearized(behavior.states == 0));

% relabel states
states = behavior.states;
idx = behavior.states == 2 | behavior.states == 4;
states(idx) = 2;

idx = behavior.states == 1 | behavior.states == 3;
states(idx) = 1;
behavior.states = states;

% update statenames
behavior.stateNames = {'middle','right','left'};

% figure to verify states
load(fullfile(basepath,[basename,'.session.mat']))
start = [];
stop = [];
for ep = 1:length(session.epochs)
    if contains(session.epochs{ep}.environment,'wmaze')
        start = [start,session.epochs{ep}.startTime];
        stop = [stop,session.epochs{ep}.stopTime];
    break
    end
end
[idx,~,~] = InIntervals(behavior.timestamps,...
    [start,stop]);

colors = {[0 0.4470 0.7410],[0.8500 0.3250 0.0980],[0.9290 0.6940 0.1250]};
figure;
subplot(1,2,1)
states = unique(behavior.states);
states = states(~isnan(states));
for s = states
   scatter(behavior.position.x(idx' &behavior.states == s),...
       behavior.position.y(idx' & behavior.states == s),[],colors{s+1})
   hold on
end
title('whole first epoch')

subplot(1,2,2)
states = unique(behavior.states);
states = states(~isnan(states));
for s = states
   scatter(behavior.timestamps(idx' & behavior.states == s),...
       behavior.position.linearized(idx' & behavior.states == s),[],colors{s+1})
   hold on
end
xlim([start,start+60*5])
title('first 5 minutes linear')


prompt = "Does this look good? y/n: ";
txt = input(prompt,"s");
if contains(txt,'y')
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end

end

% figure;
% nansum(behavior.position.x - save_test.behavior.position.x)
% 
% behavior.states = ones(1,length(behavior.timestamps));
% 
% hold on
% scatter(behavior.position.projected_x,behavior.position.projected_y)
% Xedges = floor(min(behavior.position.projected_x)):3:round(max(behavior.position.projected_x))+3
% Yedges = floor(min(behavior.position.projected_y)):3:round(max(behavior.position.projected_y))+3
% 
% [N,Xedges,Yedges] = histcounts2(behavior.position.projected_x,behavior.position.projected_y,Xedges,Yedges)
% 
% figure;
% imagesc(N,[0 1000]);
% hold on
% zlim([0 .01])
%%
% figure;
% hold on
% states = unique(behavior.states);
% states = states(~isnan(states));
% for s = states
%    scatter(behavior.position.x(behavior.states == s),...
%        behavior.position.y(behavior.states == s))
% end

% figure;
% hold on
% states = unique(behavior.states);
% states = states(~isnan(states));
% for s = states
%    scatter(behavior.timestamps(behavior.states == s),...
%        behavior.position.linearized(behavior.states == s))
% end
% 
% 
% idx = behavior.states == 2 | behavior.states == 4;
% behavior.position.linearized(idx) = behavior.position.linearized(idx) -...
%     min(behavior.position.linearized(idx)) +...
%     max(behavior.position.linearized(behavior.states == 0))

% idx = behavior.states == 0;
% 
% figure;
% scatter(behavior.timestamps(idx),behavior.position.linearized(idx))
% 
% figure;
% scatter(behavior.position.x(idx),behavior.position.y(idx))
% 
% figure;
% scatter(behavior.position.x,behavior.position.y,[],behavior.states )

% figure;
% hold on
% states = unique(behavior.states);
% states = states(~isnan(states));
% for s = states
%    disp(s) 
%    scatter(behavior.position.x(behavior.states == s),...
%        behavior.position.y(behavior.states == s))
% end

% 
% figure;
% hold on
% states = unique(behavior.states);
% states = states(~isnan(states));
% for s = states
%    scatter(behavior.timestamps(behavior.states == s),...
%        behavior.position.linearized(behavior.states == s))
% end
% 
% figure;
% hold on
% states = unique(behavior.states);
% states = states(~isnan(states));
% for s = states
%     subplot(1,5,s+1)
%    scatter(behavior.position.x(behavior.states == s),...
%        behavior.position.y(behavior.states == s))
%    title(s)
% end


% 
% myListFile = pyrunfile("python D:\github\ripple_heterogeneity\ripple_heterogeneity\utils\linearization_pipeline.py 'Z:\Data\ORproject\OR15\day1'")
% 
% command = 'python D:\github\ripple_heterogeneity\ripple_heterogeneity\utils\linearization_pipeline.py ''Z:\Data\ORproject\OR15\day1'''
% status = system(command)
% 
% disp(pyversion)
% pathToCode = fileparts(which('linearization_pipeline.py'));
% if count(py.sys.path,pathToCode)==0
%     % adds the code path to the python path
%     insert(py.sys.path,int32(0),pathToCode) 
% end
% pyOut = py.linearization_pipeline.run('Z:\Data\ORproject\OR15\day1')
% 
% 
% 
% pathToCode = fileparts(which('find_and_move_videos.py'));
% count(py.sys.path,pathToCode)==0
