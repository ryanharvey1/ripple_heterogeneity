df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'AYA') | contains(df.basepath,'AB'),:);
df = df(contains(df.basepath,'ORproject'),:);

% df = df(contains(df.basepath,'Kenji'),:);
% maze_type = 'tmaze';
% maze_sizes = 160;

% maze_type = 'cheeseboard';
% maze_sizes = 120;

maze_type = 'wmaze';
maze_sizes = 110;


for i = 1:length(df.basepath)
    disp(df.basepath{i})
    convert_it(df.basepath{i},maze_type,maze_sizes)
end

function convert_it(basepath,maze_type,maze_sizes)

basename = basenameFromBasepath(basepath);

load([basepath,filesep,[basename,'.session.mat']]);

has_maze = [];
for ep_i = 1:length(session.epochs)
    if contains(session.epochs{ep_i}.environment,maze_type)
        has_maze(ep_i) = true;
    else
        has_maze(ep_i) = false;
    end
end
if ~any(has_maze)
    return
end

load(fullfile(basepath,[basename,'.animal.behavior.mat']))

if ~isempty(behavior.position.x)
    
    startTime = [];
    stopTime = [];
    for ep = session.epochs
        if contains(ep{1}.environment,maze_type)
            startTime = [startTime;ep{1}.startTime];
            stopTime = [stopTime;ep{1}.stopTime];
        end
    end
    maze_epochs = [startTime,stopTime];
    
    figure;
    for ep = 1:size(maze_epochs,1)
        [idx,~,~] = InIntervals(behavior.timestamps,maze_epochs(ep,:));
        
        subplot(1,size(maze_epochs,1),ep)
        plot(behavior.position.x(idx),behavior.position.y(idx))
        title(session.epochs{ep}.name)
        axis equal
    end
    
    txt = input('does this need reloaded from raw? [y/n]',"s");
    if contains(txt,'y')
        behavior = fix_tracking_timestamps(basepath);
    end
    
    txt = input('does this need manual clean up? [y/n]',"s");
    
    if contains(txt,'y')
        good_idx = manual_trackerjumps(behavior.timestamps,...
            behavior.position.x,...
            behavior.position.y,...
            startTime,...
            stopTime,...
            basepath,'darkmode',false);
        
        behavior.position.x(~good_idx) = NaN;
        behavior.position.y(~good_idx) = NaN;
    end
    txt = input('does this need convert to cm? [y/n]',"s");
    if contains(txt,'y')
        for ep = 1:size(maze_epochs,1)
            [idx,~,~] = InIntervals(behavior.timestamps,maze_epochs(ep,:));
            
            figure;
            plot(behavior.position.x(idx),behavior.position.y(idx))
            
            linpos = behavior.position.x(idx);
            pos_range = max(linpos) - min(linpos);
            convert_pix_to_cm_ratio = (pos_range / maze_sizes);
            
            % convert xy to cm as well
            if ~isempty(behavior.position.x)
                behavior.position.x(idx) = behavior.position.x(idx) /...
                    convert_pix_to_cm_ratio;
                behavior.position.y(idx) = behavior.position.y(idx) /...
                    convert_pix_to_cm_ratio;
            end
        end
    end
    behavior.position.units = 'cm';
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
    close all
end
end

function behavior = fix_tracking_timestamps(basepath)

basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.tracking.behavior.mat']));
load(fullfile(basepath,[basename,'.MergePoints.events.mat']));
load(fullfile(basepath,[basename,'.animal.behavior.mat']));

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
behavior.position.x = x;
behavior.position.y = y;
behavior.position.z = z;
behavior.timestamps = t;

behavior.position.units = 'cm';
save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end
% df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'GirardeauG') ,:);
%
% % maze_size = [];
% figure;
% for i = 1:length(df.basepath)
% basepath = df.basepath{i};
% basename = basenameFromBasepath(basepath);
% load(fullfile(basepath,[basename,'.animal.behavior.mat']))
% %
%
% plot(behavior.position.x*0.43,behavior.position.y*0.43)
% hold on
% maze_size(i) = max(behavior.position.x*0.43) - min(behavior.position.x*0.43);
% end
%
% figure;
% plot(behavior.position.x)
%
% figure;
% plot(behavior.position.linearized)

% linearTrackBehavior('remove_extra_fields',true,'split_linearize',true,just_save_animal_behavior,true)

% linearTrackBehavior('maze_sizes',220,'remove_extra_fields',true,'split_linearize',true)

% close all
