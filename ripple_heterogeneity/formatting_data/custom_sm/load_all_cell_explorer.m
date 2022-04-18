% cell_explorer_all

df = readtable('Z:\home\ryanh\projects\ripple_heterogeneity\mouse_sessions.csv');
basepaths = unique(df.basepath);
for i = 1:length(basepaths)
   basenames{i} = basenameFromBasepath(basepaths{i});
end
cell_metrics = loadCellMetricsBatch('basepaths',basepaths','basenames',basenames);
cell_metrics = CellExplorer('metrics',cell_metrics);