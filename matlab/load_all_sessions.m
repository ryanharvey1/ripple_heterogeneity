% load_all_sessions
df = readtable('D:\projects\ripple_heterogeneity\sessions.csv');
basepaths = unique(df.basepath);
for i = 1:length(basepaths)
   basenames{i} = basenameFromBasepath(basepaths{i});
end
cell_metrics = loadCellMetricsBatch('basepaths',basepaths','basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);