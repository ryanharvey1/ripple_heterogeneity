

% update_channel_mapping



df = readtable('D:\projects\ripple_heterogeneity\swr_pipe_all.csv');

basepaths = unique(df.basepath);
basepaths = basepaths(contains(basepaths,'Kenji'));
for i = 1:length(basepaths)
    disp(basepaths{i})
%     if exist(fullfile(basepaths{i},'anatomical_map.csv'),'file')
%         continue
%     end
    channel_mapping('basepath',basepaths{i},...
                    'fig',true,...
                    'pull_from_cell_metrics',true,...
                    'force_cell_metric_overwrite',true)
    close all
end

% following some manual editing of the csv files, run this to add to
% .session
basepaths = basepaths(contains(basepaths,'Kenji'));
for i = 1:length(basepaths)
    disp(basepaths{i})
    channel_mapping('basepath',basepaths{i},...
                    'fig',true)
    close all
end







basepaths
'A:\Data\AB3\AB3_38_41'
channel_mapping('fig',true)