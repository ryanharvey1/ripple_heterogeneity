
basepath = 'Z:\Data\OMLproject\OML18\day1'
basename = basenameFromBasepath(basepath);
load([basepath,filesep,[basename,'.session.mat']]);

bad_idx = behavior.position.x > -150 &...
    behavior.position.x < 150 &...
    behavior.position.y > -8 &...
    behavior.position.y < 22;

x = behavior.position.x;
y = behavior.position.y;

x(~bad_idx) = NaN;
y(~bad_idx) = NaN;

figure;plot(x,y)

behavior.position.x = x;
behavior.position.y = y;

start = [];
stop = [];
for ep = 1:length(session.epochs)
    if ~contains(session.epochs{ep}.environment,'sleep')
        start = [start,session.epochs{ep}.startTime];
        stop = [stop,session.epochs{ep}.stopTime];
    end
end
    
good_idx = manual_trackerjumps(behavior.timestamps,...
    behavior.position.x,...
    behavior.position.y,...
    start,...
    stop,...
    basepath,'darkmode',true);

behavior.position.x(~good_idx) = NaN;
behavior.position.y(~good_idx) = NaN;

save(fullfile(basepath,[basename,'.animal.behavior.mat']),'behavior')

behavior = linearTrackBehavior(...
    'basepath',basepath,...
    'remove_extra_fields',true,...
    'split_linearize',false,...
    'just_save_animal_behavior',true,...
    'show_fig',false,...
    'savemat',true...
    );