df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
df = df(contains(df.basepath,'AYA') | contains(df.basepath,'AB'),:);
% df = df(contains(df.basepath,'Kenji'),:);
% maze_type = 'tmaze';
% maze_sizes = 160;

maze_type = 'cheeseboard';
maze_sizes = 120;

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
    behavior.position.units = 'cm';
    save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')
end
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
