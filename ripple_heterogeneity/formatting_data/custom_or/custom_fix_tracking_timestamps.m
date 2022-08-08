


df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'ORproject'),:);
basepath = pwd
basename = basenameFromBasepath(basepath)
load(fullfile(basepath,[basename,'.tracking.behavior.mat']))
load(fullfile(basepath,[basename,'.MergePoints.events.mat']))
load(fullfile(basepath,[basename,'.animal.behavior.mat']))

digitalIn_ttl = tracking.timestamps;
t = [];
x = [];
y = [];
z = [];
% iter over mergepoint folders to extract tracking
for k = 1:length(MergePoints.foldernames)
    % search for optitrack .tak file as there may be many .csv files
    if ~isempty(dir(fullfile(basepath,MergePoints.foldernames{k},'*.csv')))
        % locate the .tak file in this subfolder
        file = dir(fullfile(basepath,MergePoints.foldernames{k},'*.csv'));
        % use func from cellexplorer to load tracking data
        % here we are using the .csv
        try
            optitrack = optitrack2buzcode('basepath', basepath,...
                'basename', basename,...
                'filenames',fullfile(MergePoints.foldernames{k},[file.name]),...
                'unit_normalization',1,...
                'saveMat',false,...
                'saveFig',false,...
                'plot_on',false);
        catch
            warning(fullfile(MergePoints.foldernames{k},[file.name]),' IS NOT OPTITRACK FILE')
            continue
        end
        % find timestamps within current session
        ts_idx = digitalIn_ttl >= MergePoints.timestamps(k,1) & digitalIn_ttl <= MergePoints.timestamps(k,2);
        ts = digitalIn_ttl(ts_idx);
        
        % cut-to-size method of syncing ttls to frames
        try
            t = [t,ts(1:length(optitrack.position.x))'];
        catch
            t = [t,ts'];
            x = [x,optitrack.position.x(1:length(ts))];
            y = [y,optitrack.position.y(1:length(ts))];
            z = [z,optitrack.position.z(1:length(ts))];
            continue
        end
       
        
        % store xyz
        x = [x,optitrack.position.x];
        y = [y,optitrack.position.y];
        z = [z,optitrack.position.z];
    end
end
behavior.timestamps = t;
behavior.position.x = x;
behavior.position.y = y;
behavior.position.z = z;

behavior.position.units = 'cm';
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
    
% startTime = [];
% stopTime = [];
% for ep = session.epochs
%     if contains(ep{1}.environment,'wmaze')
%         startTime = [startTime;ep{1}.startTime];
%         stopTime = [stopTime;ep{1}.stopTime];
%     end
% end
% maze_epochs = [startTime,stopTime];