df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');
% df = df(contains(df.basepath,'AYA') |...
%     contains(df.basepath,'AB') |...
%     contains(df.basepath,'GirardeauG') |...
%     contains(df.basepath,'Kenji'),:);
df = df(contains(df.basepath,'Kenji'),:);

for i = 1:length(df.basepath)
    disp(df.basepath{i})
    convert_it(df.basepath{i})
end

function convert_it(basepath)

basename = basenameFromBasepath(basepath);

use_maze_size = true;
if contains(basepath,'AB1')
    maze_sizes = 120;
elseif contains(basepath,'GirardeauG')
    maze_sizes = 215;
elseif contains(basepath,'Kenji')
    maze_sizes = 250;
    use_maze_size = false;
else
    maze_sizes = 240;
end

load([basepath,filesep,[basename,'.session.mat']]);

has_linear = false;
for ep = session.epochs
    if contains(ep{1}.environment,'linear')
        has_linear = [has_linear;true];
    end
end
if ~any(has_linear)
    return
end

load(fullfile(basepath,[basename,'.animal.behavior.mat']))


if ~isempty(behavior.position.linearized) && use_maze_size
    
    startTime = [];
    stopTime = [];
    for ep = session.epochs
        if contains(ep{1}.environment,'linear')
            startTime = [startTime;ep{1}.startTime];
            stopTime = [stopTime;ep{1}.stopTime];
        end
    end
    linear_epochs = [startTime,stopTime];
    
    for ep = 1:size(linear_epochs,1)
        [idx,~,~] = InIntervals(behavior.timestamps,linear_epochs(ep,:));
        
        linpos = behavior.position.linearized(idx);
        pos_range = max(linpos) - min(linpos);
        convert_pix_to_cm_ratio = (pos_range / maze_sizes);
        linpos = linpos / convert_pix_to_cm_ratio;
        
        behavior.position.linearized(idx) = linpos';
        
        % convert xy to cm as well
        if ~isempty(behavior.position.x)
            behavior.position.x(idx) = behavior.position.x(idx) /...
                convert_pix_to_cm_ratio;
            behavior.position.y(idx) = behavior.position.y(idx) /...
                convert_pix_to_cm_ratio;
        end
    end
else
    if isempty(behavior.position.x) || all(isnan(behavior.position.x))
        return
    end
    n_linear = sum(has_linear);
    if use_maze_size
        behavior = linearTrackBehavior(...
            'basepath',basepath,...
            'remove_extra_fields',true,...
            'split_linearize',true,...
            'just_save_animal_behavior',true,...
            'show_fig',false,...
            'maze_sizes',repmat(maze_sizes,1,n_linear),...
            'savemat',false...
            );
    else
        behavior = linearTrackBehavior(...
            'basepath',basepath,...
            'remove_extra_fields',true,...
            'split_linearize',true,...
            'just_save_animal_behavior',true,...
            'show_fig',false,...
            'savemat',false...
            );
    end
end
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
