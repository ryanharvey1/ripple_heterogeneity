% custom_for_kenji_SWRunitMetrics

load('A:\Data\Kenji\ElePosition.mat')
shank_region = ElePosition(:,6:end);
for i = 1:size(shank_region,1)
    for j = 1:size(shank_region,2)
        shank_region{i,j}=lower(shank_region{i,j});
    end
end

idx = any(strcmp(shank_region,'ca1') |...
    strcmp(shank_region,'ca1c') |...
    strcmp(shank_region,'ca') |...
    strcmp(shank_region,'ca3') |...
    strcmp(shank_region,'ca2') |...
    strcmp(shank_region,'dg') |...
    strcmp(shank_region,'dgca3'),2);

sessions = ElePosition(idx,2);


data_path = 'A:\Data\Kenji\';

% keep this at top of path
addpath('D:\github\ripple_heterogeneity\matlab')
for i = 1:length(sessions)
    basename = sessions{i};
    basepath = [data_path,basename];
    disp(basepath)

    if exist(fullfile(basepath,[basename,'.SWRunitMetrics.mat']),'file')
        continue
    end
    
    load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))

    load(fullfile(basepath,[basename,'.ripples.events.mat']))

    load(fullfile(basepath,[basename '.session.mat']))

    parse_pre_task_post(session,basepath,basename,ripples,spikes)
    
end