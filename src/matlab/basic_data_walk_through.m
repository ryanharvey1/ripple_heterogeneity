% Basic walk through of the data structure
% By Ryan H - 2022

%% basepath and basename
basepath = 'Z:\Data\OMLproject\OML18\day1';
basename = basenameFromBasepath(basepath);

%% basename.session
load(fullfile(basepath,[basename,'.session.mat']))
gui_session(basepath)

%% channel_mapping
channel_mapping()

%% Cell Explorer
CellExplorer('basepath',basepath)

%% Multi-session CellExplorer
files = dir(['Z:\Data\OMLproject\OML18','\**\*.cell_metrics.cellinfo.mat']);

% pull out basepaths and basenames
for i = 1:length(files)
    basepaths{i} = files(i).folder;
    basenames{i} = basenameFromBasepath(files(i).folder);
end

% load all cell metrics
cell_metrics = loadCellMetricsBatch('basepaths',basepaths,'basenames',basenames);

% pull up gui to inspect all units in your project
cell_metrics = CellExplorer('metrics',cell_metrics);

%% cell_metrics
load(fullfile(basepath,[basename,'.cell_metrics.cellinfo.mat']))
%% spikes
load(fullfile(basepath,[basename,'.spikes.cellinfo.mat']))
%% ripples
load(fullfile(basepath,[basename,'.ripples.events.mat']))
%% animal position
load(fullfile(basepath,[basename,'.animal.behavior.mat']))
%% states
load(fullfile(basepath,[basename,'.SleepState.states.mat']))
