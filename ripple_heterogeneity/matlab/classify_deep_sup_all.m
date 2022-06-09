% classify_deep_sup_all
df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\sessions.csv');

% iter over basepaths
for i = 1:length(df.basepath)
    disp(df.basepath{i})
    run(df.basepath{i})
end

function run(basepath)
basename = basenameFromBasepath(basepath);
load(fullfile(basepath,[basename,'.session.mat']),'session');

% only re-run if chanCoords exists
if ~isfield(session.extracellular,'chanCoords')
    return
end
% only run if chanCoords is already not updated
if exist(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'file')
    load(fullfile(basepath,[basename,'.deepSuperficialfromRipple.channelinfo.mat']),'deepSuperficialfromRipple');
    if isfield(deepSuperficialfromRipple.processinginfo.params,'deepSuperficial_ChDistance_method')
        if contains(deepSuperficialfromRipple.processinginfo.params.deepSuperficial_ChDistance_method,'chanCoords')
            return
        end
    end
end
% finally, classify 
classification_DeepSuperficial(session);
close all
end