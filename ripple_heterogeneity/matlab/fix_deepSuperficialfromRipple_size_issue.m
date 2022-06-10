% find when cell arrays are not the same size and fix
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

% iter over basepaths
for i = 1:length(df.basepath)
    disp(df.basepath{i})
    run(df.basepath{i})
end

function run(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.session.mat']),'session');

% fields to check
fields = {'ripple_average','ripple_power','ripple_amplitude','SWR_diff','SWR_amplitude'};
% only run if chanCoords is already not updated
if exist(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'file')
    load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'deepSuperficialfromRipple');
    
    for i = 1:length(fields)
        n(i) = length(deepSuperficialfromRipple.(fields{i}));
    end
    if ~all(n == n(1))
        test=true
    end
end
% finally, classify 
% classification_DeepSuperficial(session);
% close all
end