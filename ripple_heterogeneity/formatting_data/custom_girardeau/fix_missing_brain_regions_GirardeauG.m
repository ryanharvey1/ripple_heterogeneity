% fix_missing_brain_regions_GirardeauG

%% locate sessions
data_path = 'A:\Data\GirardeauG';
tic
files = dir([data_path,'\**\*.cell_metrics.cellinfo.mat']);
disp(toc)

for i = 1:length(files)
   basepath{i} = files(i).folder;
   basename{i} = basenameFromBasepath(files(i).folder);
end
%% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepath);

%% load shank location info
load(fullfile('A:\Data\GirardeauG\Rat11','RecordingLocations.mat'))

% iter through each rat and label brain regions by RecordingLocations.mat
rats = unique(cell_metrics.animal);
for i = 1:length(rats)
    disp(rats{i})
    for area = 1:11
       disp(shankInfo{1,area})
       shanks = shankInfo{i+1,area};
       if ~isempty(shanks)
           idx = contains(cell_metrics.animal,rats{i}) &...
               ismember(cell_metrics.shankID,unique(shanks(:,2)));
           cell_metrics.brainRegion(idx) = {shankInfo{1,area}};
       end
    end
end
% relabel hpc to ca1 (all hpc in this data is ca1)
cell_metrics.brainRegion(contains(cell_metrics.brainRegion,'hpc')) = {'CA1'};
% print out uniqe regions for a check
unique(cell_metrics.brainRegion)

%% load into CellExplorer and then save all to update CellMetrics
cell_metrics = CellExplorer('metrics',cell_metrics);

%%



% 
% % 
% cell_metrics.brainRegion(cell_metrics.shankID < 5 & contains(cell_metrics.sessionName,'Rat08')) = {'CA1'};
% 
% cell_metrics.brainRegion(cell_metrics.shankID > 16 & ~contains(cell_metrics.sessionName,'Rat08')) = {'CA1'};

% cell_metrics.shankID(contains(cell_metrics.brainRegion,"HIP"));


% unique(cell_metrics.brainRegion(cell_metrics.shankID < 5 & contains(cell_metrics.sessionName,'Rat08')))



