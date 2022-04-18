% update_with_bz_getSessionInfo
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

% loop through each session
for i = 1:length(sessions)
    basename = sessions{i};
    if ~exist(basepath,'dir')
        continue
    end
    basepath = [data_path,basename];
    disp(basepath)
    if ~exist(fullfile(basepath,[basename,'.sessionInfo.mat']),'file')
        sessionInfo = bz_getSessionInfo(basepath,'saveMat',true,...
            'noPrompts',true);
    end
end