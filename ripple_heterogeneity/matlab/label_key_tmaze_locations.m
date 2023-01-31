% label_key_tmaze_locations
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

% locate tmaze sessions
epochs = [];
for i = 1:length(df.basepath)
    basepath = df.basepath{i};
    basename = basenameFromBasepath(basepath);
    load(fullfile(basepath,[basename,'.session.mat']),'session')
    
    epochs = [epochs;load_epoch('basepath', basepath)];
end
% locate all tmaze sessions
basepaths = unique(epochs(ismember(epochs.environment,["tmaze","Mwheel","Tmaze"]),:).basepath);


% 
% start = [];
% stop = [];
% for ep = 1:length(session.epochs)
%     if ~contains(session.epochs{ep}.environment,'sleep')
%         start = [start,session.epochs{ep}.startTime];
%         stop = [stop,session.epochs{ep}.stopTime];
%     end
% end
% epochs = IntervalArray([start',stop']);
% 
% % let the user click around the coordinates
% while true
%     [X,Y]=ginput(1);
%     if isempty(X)
%         break
%     end
%     corners(i,:)=[X,Y];
%     plot(corners(:,1),corners(:,2),'r',X,Y,'*r')
%     i=i+1;
% end