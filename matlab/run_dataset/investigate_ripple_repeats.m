
% investigate_ripple_repeats

% ts_0_diff_ = diff(ripples.timestamps(:,1)) == 0;
% 
% sum(ts_0_diff_)

% n_repeated_start = length(ripples.timestamps(:,1)) - length(unique(ripples.timestamps(:,1)));
% 
% n_repeated_stop = length(ripples.timestamps(:,2)) - length(unique(ripples.timestamps(:,2)));
% 
% n_repeated_peak = length(ripples.peaks) - length(unique(ripples.peaks));

df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);


for i = 1:length(basepaths)
    basepath = basepaths{i};
    disp(basepath)
    [repeats(i,:),detectorname{i,1}] = main(basepath);
end


sum(repeats(:,1)>0)
min(repeats(repeats(:,1)>0,1))
max(repeats(repeats(:,1)>0,1))

median(repeats(repeats(:,1)>0,1))

unique(detectorname(repeats(:,1)>0))



function [repeats,detectorname] = main(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.ripples.events.mat']))

% dir(fullfile(basepath,[basename,'.ripples.events.mat']))

n_repeated_start = length(ripples.timestamps(:,1)) - length(unique(ripples.timestamps(:,1)));

n_repeated_peak = length(ripples.peaks) - length(unique(ripples.peaks));

n_repeated_stop = length(ripples.timestamps(:,2)) - length(unique(ripples.timestamps(:,2)));

repeats = [n_repeated_start,n_repeated_peak,n_repeated_stop];

try
    detectorname = ripples.detectorinfo.detectorname;  
catch
    detectorname = ripples.detectorName;
end
end